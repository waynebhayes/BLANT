BUILD_DIR=build
TARGET=$(BUILD_DIR)/cdijkstra
HEADERS=$(wildcard include/*.h)
SOURCES=$(wildcard src/*.c)
LIBWAYNE_HOME=../libwayne

CC=gcc
CFLAGS=-g -Iinclude -I$(LIBWAYNE_HOME)/include
CLIBS=-L$(LIBWAYNE_HOME) -lwayne -lm

.PHONY: all clean

all: $(TARGET)

clean:
	-rm -rf $(BUILD_DIR)

$(TARGET): $(BUILD_DIR) libwayne $(HEADERS) $(SOURCES) 
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(CLIBS)

libwayne:
	cd .. && $(MAKE) $@

$(BUILD_DIR):
	mkdir -p $@
