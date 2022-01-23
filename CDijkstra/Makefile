BUILD_DIR=build
TARGET=$(BUILD_DIR)/cdijkstra
SRCS=$(wildcard *.cpp)

CC=g++
CFLAGS=-g

.PHONY: all clean

all: $(TARGET)

clean:
	-rm -rf $(BUILD_DIR)

$(TARGET): $(BUILD_DIR) $(SRCS)
	$(CC) $(CFLAGS) -o $@ $(SRCS)

$(BUILD_DIR):
	mkdir -p $@
