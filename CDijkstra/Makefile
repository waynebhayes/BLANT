BUILD_DIR=build
TARGET=$(BUILD_DIR)/cdijkstra
HEADERS=$(wildcard include/*.h)
SOURCES=$(wildcard src/*.cpp) include/matrix.h
LIBWAYNE_HOME=../libwayne

CXX=g++
CXXFLAGS=-g -std=c++11 -Iinclude

.PHONY: all clean

all: $(TARGET)

clean:
	-rm -rf $(BUILD_DIR)

$(TARGET): $(BUILD_DIR) $(HEADERS) $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES)

$(BUILD_DIR):
	mkdir -p $@
