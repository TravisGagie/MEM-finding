// Backward search (Section 2, "BWTs"): a BWT interval [sp, ep) and a step that prepends one character,
// on the RLBWT or on the RLCSA (both give the same intervals).  The text holds both strands, so the
// forward extension of a match R[a..e) is the backward extension of rc(R[a..e)) by the complement of R[e].
//
// The lookup table (-F t; KATKA App. A, "Lookup tables") stores the interval of every string of t bases,
// so a search can start with one lookup instead of t steps.  When the string does not occur, the caller
// falls back to stepping base by base, so BML still learns exactly where the search fails, and the
// answers are the same with and without a table.
#pragma once
#include "rlcsa.hpp"

namespace katka {

struct interval {
    u64 sp = 0, ep = 0;
    bool operator==(const interval &o) const { return sp == o.sp && ep == o.ep; }
};

struct stepper {
    const rlbwt *bwt = nullptr;
    const rlcsa *csa = nullptr;     // if set, steps use it instead of bwt
    u64 n() const { return csa ? csa->n : bwt->n; }
    interval full() const { return {0, n()}; }
    // [sp, ep) for P -> [sp, ep) for cP; returns false, leaving I unchanged, if cP does not occur
    bool step(interval &I, unsigned char c) const {
        u64 a = I.sp, b = I.ep;
        if (csa) csa->extend(a, b, c); else bwt->extend(a, b, c);
        if (a >= b) return false;
        I = {a, b};
        return true;
    }
};

// The table is over an alphabet: A, C, G, T by default, or the symbols of a digest (set_alphabet).
struct lookup_table {
    int t = 0;                      // 0: no table
    std::vector<unsigned char> alpha{'A', 'C', 'G', 'T'};
    int16_t code[256];
    std::vector<interval> iv;       // indexed by the string's code in base |alpha|, first symbol most significant
    std::vector<uint8_t> occurs;

    lookup_table() { set_alphabet(alpha); }
    void set_alphabet(const std::vector<unsigned char> &a) {
        alpha = a;
        std::fill(code, code + 256, -1);
        for (size_t i = 0; i < alpha.size(); ++i) code[alpha[i]] = (int16_t)i;
    }
    void build(const stepper &S, int tt) {
        t = tt;
        const u64 sig = alpha.size();
        long double total = 1;
        for (int i = 0; i < t; ++i) total *= sig;
        if (total > 5e8) die("lookup table too large (" + std::to_string((double)total) + " entries)");
        std::vector<interval> a{S.full()};
        std::vector<uint8_t> ok{1};
        for (int i = 0; i < t; ++i) {        // the strings of length i + 1 are those of length i with a symbol in front
            u64 m = a.size();
            std::vector<interval> a2(sig * m);
            std::vector<uint8_t> ok2(sig * m, 0);
            for (u64 c = 0; c < sig; ++c)
                for (u64 w = 0; w < m; ++w)
                    if (ok[w]) { interval I = a[w]; if (S.step(I, alpha[c])) { a2[c * m + w] = I; ok2[c * m + w] = 1; } }
            a.swap(a2); ok.swap(ok2);
        }
        iv.swap(a); occurs.swap(ok);
    }
    u64 bytes() const { return iv.size() * (sizeof(interval) + 1); }
    // the interval of p[0..t); false if it does not occur (or has a symbol outside the alphabet)
    bool get(const unsigned char *p, interval &I) const {
        u64 key = 0;
        for (int i = 0; i < t; ++i) { int c = code[p[i]]; if (c < 0) return false; key = alpha.size() * key + c; }
        if (!occurs[key]) return false;
        I = iv[key];
        return true;
    }
};

}  // namespace katka
