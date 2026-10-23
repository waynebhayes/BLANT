CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Werror
SRCS      = options.cpp network.cpp extsort.cpp samples.cpp cluster.cpp
OBJS      = $(SRCS:.cpp=.o)

all: blant-clusters

blant-clusters: blant-clusters.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

test/selftest: test/selftest.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

test: test/selftest
	./test/selftest

clean:
	rm -f *.o test/*.o blant-clusters test/selftest test/check

.PHONY: all test clean
