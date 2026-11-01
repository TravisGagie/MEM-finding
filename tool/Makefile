# Builds bin/katka and the tests.  No dependencies besides zlib (libsais is in third_party/).
#   make            (CXXFLAGS=... to change the flags; the default uses -march=native)
CXX ?= g++
CC ?= gcc
CXXFLAGS ?= -std=c++17 -O3 -march=native -DNDEBUG -Wall -Wno-unused-function
CFLAGS ?= -O3 -march=native -DNDEBUG
override CXXFLAGS += -Ithird_party/libsais/include -pthread
LDLIBS = -lz -pthread
HDRS = $(wildcard src/*.hpp)
LIBSAIS = build/libsais.o build/libsais64.o

all: bin/katka bin/katka-test

build/%.o: third_party/libsais/src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -Ithird_party/libsais/include -c $< -o $@

bin/katka: src/main.cpp $(HDRS) $(LIBSAIS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) src/main.cpp $(LIBSAIS) -o $@ $(LDLIBS)

bin/katka-test: tests/test.cpp $(HDRS) $(LIBSAIS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Isrc tests/test.cpp $(LIBSAIS) -o $@ $(LDLIBS)

clean:
	rm -rf bin build

.PHONY: all clean
