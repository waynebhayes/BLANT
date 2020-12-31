CC=gcc
mt19937: MTGenerator.hpp c_code.c mt19937.cpp mt19937.h
	$(CC) -c c_code.c
	$(CXX) -std=c++11 -c mt19937.cpp
	$(CXX) -o mt19937 c_code.o mt19937.o

clean:
	@/bin/rm -f *.o mt19937
