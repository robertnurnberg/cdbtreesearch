TERARKDBROOT = /home/vondele/chess/noob/terarkdb
CDBDIRECTROOT = /home/vondele/chess/vondele/cdbdirect

LDFLAGS = -L$(TERARKDBROOT)/output/lib -L$(CDBDIRECTROOT)
LIBS = -lcdbdirect -lterarkdb -lterark-zip-r -lboost_fiber -lboost_context -ltcmalloc -pthread -lgcc -lrt -ldl -ltbb -laio -lgomp -lsnappy -llz4 -lz -lbz2

all: cdbtreesearch cdbtreesearchnosoft

CXXFLAGS = -std=c++20 -O3 -g -march=native -fno-omit-frame-pointer -fno-inline
CXXFLAGS = -std=c++20 -O3 -g -march=native

cdbtreesearchnosoft: main.cpp
	g++ $(CXXFLAGS) -I$(CDBDIRECTROOT) -o cdbtreesearchnosoft main.cpp $(LDFLAGS) $(LIBS)

cdbtreesearch: main.cpp
	g++ $(CXXFLAGS) -I$(CDBDIRECTROOT) -o cdbtreesearch softmax.cpp $(LDFLAGS) $(LIBS)

clean:
	rm -f cdbtreesearch cdbtreesearchnosoft

format:
	clang-format -i main.cpp
