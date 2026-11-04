// Small helpers shared by the whole tool: types, bases, errors, timing, binary I/O.
#pragma once
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace katka {

typedef uint64_t u64;

inline double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// A, C, G, T -> 0, 1, 2, 3; anything else -> -1
inline int base_code(unsigned char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}

// complement of a base; anything else is its own complement
inline unsigned char complement(unsigned char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return c; }
}

inline std::string reverse_complement(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (auto &c : r) c = (char)complement((unsigned char)c);
    return r;
}

[[noreturn]] inline void die(const std::string &msg) {
    fprintf(stderr, "katka: %s\n", msg.c_str());
    exit(1);
}

// binary I/O on FILE*, throwing on short reads
template <class T> inline void wr(FILE *f, const T &x) { fwrite(&x, sizeof(T), 1, f); }
template <class T> inline void rd(FILE *f, T &x) {
    if (fread(&x, sizeof(T), 1, f) != 1) throw std::runtime_error("short read");
}
template <class V> inline void wvec(FILE *f, const V &v) {
    uint64_t k = v.size(); fwrite(&k, 8, 1, f);
    if (k) fwrite(v.data(), sizeof(v[0]), k, f);
}
template <class V> inline void rvec(FILE *f, V &v) {
    uint64_t k; rd(f, k);
    v.resize(k);
    if (k && fread(v.data(), sizeof(v[0]), k, f) != k) throw std::runtime_error("short read");
}
inline void wmagic(FILE *f, const char *m) { fwrite(m, 1, 8, f); }
inline void rmagic(FILE *f, const char *m, const std::string &file) {
    char b[8];
    if (fread(b, 1, 8, f) != 8 || memcmp(b, m, 8)) die(file + ": not a " + std::string(m, 8) + " file (rebuild it)");
}
inline FILE *open_or_die(const std::string &file, const char *mode) {
    FILE *f = fopen(file.c_str(), mode);
    if (!f) die("cannot open " + file);
    return f;
}
inline bool file_exists(const std::string &f) {
    FILE *x = fopen(f.c_str(), "rb");
    if (x) fclose(x);
    return x != nullptr;
}
inline u64 file_size(const std::string &file) {
    FILE *f = open_or_die(file, "rb");
    fseeko(f, 0, SEEK_END);
    u64 n = ftello(f);
    fclose(f);
    return n;
}

// number of equal leading bytes of a[0..len) and b[0..len)
inline u64 lce_bytes(const uint8_t *a, const uint8_t *b, u64 len) {
    u64 i = 0;
    for (; i + 8 <= len; i += 8) {
        uint64_t x, y;
        memcpy(&x, a + i, 8); memcpy(&y, b + i, 8);
        if (x != y) return i + (__builtin_ctzll(x ^ y) >> 3);
    }
    while (i < len && a[i] == b[i]) ++i;
    return i;
}

}  // namespace katka
