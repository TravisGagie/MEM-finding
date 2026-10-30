# Builds the six programs of the paper's pipeline into bin/.
#   make SDSL=/path/to/sdsl     (the prefix sdsl-lite's install.sh was given; default ../deps/sdsl)
SDSL ?= ../deps/sdsl
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -DNDEBUG
override CXXFLAGS += -isystem $(SDSL)/include -pthread   # -isystem: no warnings from sdsl's own headers
LDLIBS = -L$(SDSL)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lz

PROGS = katka-prep katka-build katka-tagbuild katka-csabuild katka-kbbuild katka-classify
HDRS = $(wildcard src/*.hpp)

all: $(addprefix bin/,$(PROGS))

bin/katka-prep: src/prep.cpp $(HDRS)
bin/katka-build: src/build.cpp $(HDRS)
bin/katka-tagbuild: src/tagbuild.cpp $(HDRS)
bin/katka-csabuild: src/csabuild.cpp $(HDRS)
bin/katka-kbbuild: src/kbbuild.cpp $(HDRS)
bin/katka-classify: src/classify.cpp $(HDRS)

bin/%:
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -rf bin

.PHONY: all clean
