BUILD_DIR=build
TARGET=$(BUILD_DIR)/cdijkstra
SRCS=$(wildcard *.cpp)

CC=g++
CFLAGS=-g

.PHONY: all

all: $(TARGET)

$(TARGET): $(BUILD_DIR) $(SRCS)
	$(CC) $(CFLAGS) -o $@ $(SRCS)

$(BUILD_DIR):
	mkdir -p $@
