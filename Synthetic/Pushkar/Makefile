CXX = g++
CXXFLAGS = -std=c++98 -O3  -fopenmp
INCLUDES = -I Snap-6.0/snap-core -I Snap-6.0/glib-core -I Snap-6.0/snap-adv
LDFLAGS = -L Snap-6.0/snap-core
LDLIBS = -l Snap.o
SOURCES = main.cpp Snap-6.0/snap-core/Snap.cpp
TARGET = main

all:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET)

