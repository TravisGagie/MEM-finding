// Small helpers shared by all the programs.
#pragma once
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace katka {

typedef uint64_t u64;

// wall-clock time in seconds
inline double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// A, C, G, T -> 0, 1, 2, 3; anything else -> -1
inline int base_code(unsigned char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}

// complement of a base; anything else (e.g. the separator X) is its own complement
inline unsigned char complement(unsigned char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return c; }
}

[[noreturn]] inline void die(const std::string &msg) {
    fprintf(stderr, "%s\n", msg.c_str());
    exit(1);
}

}  // namespace katka
