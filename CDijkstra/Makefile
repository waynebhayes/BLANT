BUILD_DIR = build
TARGET = $(BUILD_DIR)/cdijkstra
HEADERS = $(wildcard include/*.h)
SOURCES = $(wildcard src/*.cpp)
LIBWAYNE_HOME = ../libwayne

CXX = g++
CXXFLAGS = -g -std=c++11 -Iinclude -I$(LIBWAYNE_HOME)/include
LDFLAGS = -L$(LIBWAYNE_HOME) -lwayne -lm

.PHONY: all clean

all: $(TARGET)

clean:
	-rm -rf $(BUILD_DIR)

$(TARGET): $(BUILD_DIR) $(HEADERS) $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

libwayne:
	cd $(LIBWAYNE_HOME) && $(MAKE)

$(BUILD_DIR):
	mkdir -p $@
