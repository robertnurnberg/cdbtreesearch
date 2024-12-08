#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string_view>
#include <unistd.h>

#include "external/ThreadPool.h"
#include "external/chess.hpp"
#include "external/parallel_hashmap/phmap.h"

#include "cdbdirect.h"

using namespace chess;

using PackedBoard = std::array<std::uint8_t, 24>;

namespace std {
template <> struct hash<PackedBoard> {
  size_t operator()(const PackedBoard pbfen) const {
    std::string_view sv(reinterpret_cast<const char *>(pbfen.data()),
                        pbfen.size());
    return std::hash<std::string_view>{}(sv);
  }
};
} // namespace std

using fen_map_t = phmap::parallel_flat_hash_map<
    PackedBoard, double, std::hash<PackedBoard>, std::equal_to<PackedBoard>,
    std::allocator<std::pair<PackedBoard, double>>, 8, std::mutex>;

// get memory in MB
std::pair<size_t, size_t> get_memory() {
  size_t tSize = 0, resident = 0, share = 0;
  std::ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident;
  buffer.close();

  long page_size = sysconf(_SC_PAGE_SIZE);
  return std::make_pair(tSize * page_size / (1024 * 1024),
                        resident * page_size / (1024 * 1024));
};

// locate command line arguments
inline bool find_argument(const std::vector<std::string> &args,
                          std::vector<std::string>::const_iterator &pos,
                          std::string_view arg,
                          bool without_parameter = false) {
  pos = std::find(args.begin(), args.end(), arg);

  return pos != args.end() &&
         (without_parameter || std::next(pos) != args.end());
}

// local data and time
std::string getCurrentDateTime() {
  // Get current time
  std::time_t now = std::time(nullptr);

  // Convert it to local time structure
  std::tm *localTime = std::localtime(&now);

  // Create a stringstream to format the date and time
  std::ostringstream dateTimeStream;
  dateTimeStream << (1900 + localTime->tm_year) << "-"
                 << (localTime->tm_mon + 1) << "-" << localTime->tm_mday << " "
                 << localTime->tm_hour << ":" << localTime->tm_min << ":"
                 << localTime->tm_sec;

  // Convert to string and return
  return dateTimeStream.str();
};

struct Stats {
  std::atomic<size_t> gets;
  std::atomic<size_t> hits;
  std::atomic<size_t> nodes;
  std::atomic<size_t> max_ply;

  void clear() { gets = hits = nodes = max_ply = 0; };
};

template <typename T>
T atomic_fetch_max_explicit(
    std::atomic<T> &pv, typename std::atomic<T>::value_type v,
    std::memory_order m = std::memory_order_seq_cst) noexcept {
  auto t = pv.load(m);
  while (std::max(v, t) != t) {
    if (pv.compare_exchange_weak(t, v, m, m))
      break;
  }
  return t;
}

using score_t = double;
using depth_t = int;

using search_result_t = std::pair<score_t, depth_t>;

search_result_t gameOver(Board &board, int depth) {

  // game ends, no attempt for now to delay mate as long as possible, i.e. all
  // mates scored the same.
  auto res = board.isGameOver();
  if (res.second == GameResult::LOSE)
    return std::make_pair(double(-30000), depth);

  if (res.second == GameResult::DRAW)
    return std::make_pair(double(0), depth);

  return std::make_pair(double(0), -1);
}

search_result_t softmax(Board &board, int depth, double base, bool deriv,
                        fen_map_t &fen_map, Stats &stats, const double alpha,
                        const std::uintptr_t handle, const bool showTree) {

  stats.nodes++;
  search_result_t ret = gameOver(board, depth);
  if (ret.second >= 0)
    return ret;

  // probe DB
  stats.gets++;
  std::vector<std::pair<std::string, int>> result =
      cdbdirect_get(handle, board.getFen(false));
  size_t n_elements = result.size();
  int ply = result[n_elements - 1].second;

  if (ply == -2 || n_elements <= 1) {

    // non-existent node, let's maintain its deriv info
    if (deriv) {
      PackedBoard pbfen = Board::Compact::encode(board);
      fen_map.lazy_emplace_l(
          std::move(pbfen),
          [&pbfen, &base](fen_map_t::value_type &p) { p.second += base; },
          [&pbfen, &base](const fen_map_t::constructor &ctor) {
            ctor(std::move(pbfen), base);
          });
    }

    return std::make_pair(0, -1);
  }

  stats.hits++;

  // collect for all the cdb moves the softmax values
  std::vector<std::pair<std::string, double>> softmaxes;
  softmaxes.reserve(result.size() - 1);

  double max = -30001;

  // get the correct softmax values for all moves, and the maximum value
  for (auto &pair : result) {
    if (pair.first != "a0a0") {
      if (depth == 0) {
        softmaxes.emplace_back(std::make_pair(pair.first, double(pair.second)));
      } else {
        Move m = uci::uciToMove(board, pair.first);
        board.makeMove<true>(m);
        search_result_t moveresult =
            softmax(board, depth - 1, base, false, fen_map, stats, alpha,
                    handle, showTree);
        board.unmakeMove(m);
        if (moveresult.second < 0)
          softmaxes.emplace_back(
              std::make_pair(pair.first, double(pair.second)));
        else
          softmaxes.emplace_back(std::make_pair(pair.first, -moveresult.first));
      }
      if (softmaxes.back().second > max)
        max = softmaxes.back().second;
    }
  }

  std::sort(softmaxes.begin(), softmaxes.end(),
            [](const std::pair<std::string, double> &x,
               const std::pair<std::string, double> &y) {
              return x.second > y.second;
            });

  // compute softmax
  // for numerical stability, we subtract max from all arguments, and add back
  // as needed (leaving softmax invariant).
  double normalize = 0;
  double average = 0;
  for (auto &pair : softmaxes) {
    double weight = std::exp(alpha * (pair.second - max));
    normalize += weight;
    average += (pair.second - max) * weight;
  }
  double softmaxvalue = max + average / normalize;

  // if requested propagate derived.
  if (deriv) {
    for (auto &pair : softmaxes) {
      double weight = std::exp(alpha * (pair.second - max)) / normalize;
      double deriv = weight * (1 + alpha * (pair.second - softmaxvalue));

      if (showTree) {
        int ply =
            board.fullMoveNumber() * 2 - (board.sideToMove() == Color::WHITE);
        std::cout << std::string(ply * 4, ' ') << std::setw(6) << pair.first
                  << " : " << std::fixed << std::setw(6) << std::setprecision(2)
                  << pair.second << " " << std::scientific << std::setw(12)
                  << std::setprecision(2) << base * deriv
                  << (depth == 0 ? "!" : "") << std::endl;
      }

      if (depth > 0) {
        Move m = uci::uciToMove(board, pair.first);
        board.makeMove<true>(m);
        softmax(board, depth - 1, -base * deriv, true, fen_map, stats, alpha,
                handle, showTree);
        board.unmakeMove(m);
      }
    }
  }

  return std::make_pair(softmaxvalue, depth);
}

// 8 levels of threadpool
std::array<ThreadPool, 8> threadpools{};

//
// approximate derivatives are computed based on cdb evals of the node,
// rather than first computing the exact softmaxes of the moves.
//
search_result_t cdbsoftmax(Board &board, double base, size_t budget,
                           fen_map_t &fen_map, Stats &stats, const double alpha,
                           const std::uintptr_t handle, const bool showTree) {

  stats.nodes++;
  size_t ply =
      board.fullMoveNumber() * 2 - (board.sideToMove() == Color::WHITE);
  atomic_fetch_max_explicit(stats.max_ply, ply);

  search_result_t ret = gameOver(board, 0);
  if (ret.second >= 0)
    return ret;

  // probe DB
  stats.gets++;
  budget--;
  std::vector<std::pair<std::string, int>> result =
      cdbdirect_get(handle, board.getFen(false));
  size_t n_elements = result.size();
  ply = result[n_elements - 1].second;

  if (ply == -2 || n_elements <= 1) {

    // non-existent node, let's maintain its deriv info
    PackedBoard pbfen = Board::Compact::encode(board);
    fen_map.lazy_emplace_l(
        std::move(pbfen),
        [&pbfen, &base](fen_map_t::value_type &p) { p.second += base; },
        [&pbfen, &base](const fen_map_t::constructor &ctor) {
          ctor(std::move(pbfen), base);
        });

    return std::make_pair(0, -1);
  }

  stats.hits++;

  // compute softmax based on cdb values, we'll use that for deriv.
  // for numerical stability, we subtract max from all arguments, and add back
  // as needed (leaving softmax invariant).
  double max = result[0].second;
  double normalize = 0;
  double average = 0;
  for (auto &pair : result) {
    if (pair.first != "a0a0") {
      double weight = std::exp(alpha * (pair.second - max));
      normalize += weight;
      average += (pair.second - max) * weight;
    }
  }
  double softmaxvalue = max + average / normalize;

  // collect for all the cdb moves the real softmax values
  std::vector<std::future<std::pair<std::string, double>>> fsoftmaxes;
  fsoftmaxes.reserve(result.size() - 1);

  double remaining = 1.0;

  // propage the deriv
  for (auto &pair : result) {
    if (pair.first != "a0a0") {
      double weight = std::exp(alpha * (pair.second - max)) / normalize;
      double deriv = weight * (1 + alpha * (pair.second - softmaxvalue));
      remaining -= weight;
      size_t next_budget = budget - std::llrint(budget * remaining);
      budget -= next_budget;
      bool follow = next_budget > 0;

      ply = board.fullMoveNumber() * 2 - (board.sideToMove() == Color::WHITE);

      if (showTree) {
        std::cout << std::string(ply * 4, ' ') << std::setw(6) << pair.first
                  << " : " << std::fixed << std::setw(6) << std::setprecision(2)
                  << pair.second << " " << std::scientific << std::setw(12)
                  << std::setprecision(2) << base * deriv
                  << (!follow ? " ! " : " ") << budget << " " << next_budget
                  << std::endl;
      }

      if (follow) {

        Move m = uci::uciToMove(board, pair.first);
        board.makeMove<true>(m);
        // needs proper game ply.
        std::string fen = board.getFen();

        auto f = [pair, fen, base, deriv, next_budget, &fen_map, &stats, &alpha,
                  &handle, &showTree]() {
          Board board(fen);
          search_result_t moveresult =
              cdbsoftmax(board, -base * deriv, next_budget, fen_map, stats,
                         alpha, handle, showTree);
          return std::make_pair(pair.first, moveresult.second < 0
                                                ? double(pair.second)
                                                : -moveresult.first);
        };

        // if we have a threadpool at this depth, and the cost of the task is
        // sufficient, schedule it on this threadpool
        if (ply - 2 < threadpools.size() && next_budget > 100) {
          fsoftmaxes.emplace_back(threadpools[ply - 2].enqueue(f));
        } else {
          std::promise<std::pair<std::string, double>> res;
          res.set_value(f());
          fsoftmaxes.emplace_back(res.get_future());
        }

        board.unmakeMove(m);

      } else {

        std::promise<std::pair<std::string, double>> res;
        res.set_value(std::make_pair(pair.first, double(pair.second)));
        fsoftmaxes.emplace_back(res.get_future());
      }
    }
  }

  std::vector<std::pair<std::string, double>> softmaxes;
  softmaxes.reserve(result.size() - 1);
  for (auto &res : fsoftmaxes)
    softmaxes.push_back(res.get());

  // find actual maximum of returned values
  max = std::numeric_limits<decltype(max)>::lowest();
  for (auto &[k, v] : softmaxes)
    max = std::max(max, v);

  normalize = 0;
  average = 0;
  for (auto &pair : softmaxes) {
    double weight = std::exp(alpha * (pair.second - max));
    normalize += weight;
    average += (pair.second - max) * weight;
  }
  softmaxvalue = max + average / normalize;

  // debugging, this should never happen.
  if (std::isnan(softmaxvalue)) {
    std::cout << std::endl
              << "Size: " << softmaxes.size() << " max: " << max
              << "normalize: " << normalize << "average: " << average
              << std::endl;
    for (auto &pair : softmaxes)
      std::cout << pair.second << std::endl;
    std::cout << std::endl << std::endl;
    assert(false);
  }

  /* softmaxvalue = std::round(softmaxvalue); */

  return std::make_pair(softmaxvalue, 0);
}

int main(int argc, char const *argv[]) {

  const std::vector<std::string> args(argv + 1, argv + argc);
  std::vector<std::string>::const_iterator pos;

  // starting fen, defaults to 1. g4
  std::string fen = "rnbqkbnr/pppppppp/8/8/6P1/8/PPPPPP1P/RNBQKBNR b KQkq -";

  if (find_argument(args, pos, "--fen"))
    fen = {*std::next(pos)};

  if (fen == "startpos")
    fen = constants::STARTPOS;

  std::vector<std::string> inputFens;

  if (find_argument(args, pos, "--fenFile")) {
    std::string filename = {*std::next(pos)};
    std::ifstream infile(filename);
    while (std::getline(infile, fen))
      inputFens.push_back(fen);
    infile.close();
  } else {
    inputFens.push_back(fen);
  }

  // ensure without move counters
  for (size_t i = 0; i < inputFens.size(); ++i)
    inputFens[i] = Board(inputFens[i]).getFen(false);

  // (initial) depth for searches
  int depth = 0;
  if (find_argument(args, pos, "--depth"))
    depth = std::stoi(*std::next(pos));

  double alpha = 1;
  if (find_argument(args, pos, "--alpha"))
    alpha = std::stod(*std::next(pos));

  if (find_argument(args, pos, "--inverseAlpha"))
    alpha = 1.0 / std::stod(*std::next(pos));

  size_t budget = 1000000;
  if (find_argument(args, pos, "--budget"))
    budget = std::stoll(*std::next(pos));

  bool showTree = find_argument(args, pos, "--showTree", true);
  bool exact = find_argument(args, pos, "--exact", true);
  bool allmoves = find_argument(args, pos, "--moves", true);

  std::cout << "Opening DB" << std::endl;
  std::uintptr_t handle = cdbdirect_initialize(CHESSDB_PATH);

  Stats stats;

  std::cout << "Starting softmax search with alpha " << alpha << std::endl;
  std::cout << "depth " << std::setw(12) << depth << std::endl;
  std::cout << "budget " << std::setw(12) << budget << std::endl;

  std::vector<std::pair<std::string, std::string>> fens;

  for (std::string &fen : inputFens) {
    if (allmoves) {
      Board board(fen);
      std::vector<std::pair<std::string, int>> result =
          cdbdirect_get(handle, fen);
      for (auto &[move, score] : result) {
        if (move != "a0a0") {
          Move m = uci::uciToMove(board, move);
          board.makeMove<true>(m);
          fens.push_back(std::make_pair(board.getFen(false), move));
          board.unmakeMove(m);
        }
      }
    } else {
      fens.push_back(std::make_pair(fen, "none"));
    }
  }

  std::vector<std::pair<PackedBoard, double>> unknown_fens;

  for (auto &[fen, move] : fens) {
    std::cout << "Search fen:" << std::setw(70) << fen << " after move "
              << std::setw(5) << move << std::flush;
    Board board(fen);
    stats.clear();
    double base = -1.0;
    bool deriv = true;
    fen_map_t fen_map;
    search_result_t sr;

    if (exact) {
      sr = softmax(board, depth, base, deriv, fen_map, stats, alpha, handle,
                   showTree);
    } else {
      sr = cdbsoftmax(board, base, budget, fen_map, stats, alpha, handle,
                      showTree);
    };
    std::cout << " score " << std::fixed << std::setw(8) << std::setprecision(2)
              << sr.first;
    std::cout << " unknown " << std::setw(10) << fen_map.size();
    std::cout << " hits    " << std::setw(10) << stats.hits;
    std::cout << " max_ply " << std::setw(10) << stats.max_ply << std::endl;

    for (auto &v : fen_map)
      unknown_fens.push_back(v);
  }

  std::sort(unknown_fens.begin(), unknown_fens.end(),
            [](const std::pair<PackedBoard, double> &x,
               const std::pair<PackedBoard, double> &y) {
              return std::abs(x.second) > std::abs(y.second);
            });

  // store unknown fens
  std::ofstream ufile("unknown.epd");
  assert(ufile.is_open());

  for (auto &[pbfen, deriv] : unknown_fens)
    ufile << std::setw(70) << Board::Compact::decode(pbfen).getFen(false)
          << " : " << std::scientific << std::setw(12) << std::setprecision(2)
          << deriv << std::endl;

  ufile.close();

  std::cout << "Closing DB" << std::endl;
  handle = cdbdirect_finalize(handle);

  return 0;
}
