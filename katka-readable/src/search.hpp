// Backward search (Section 5): a BWT interval [sp, ep) and one step that prepends a character, on the
// RLBWT or, with katka-classify -C, on the RLCSA (both give the same intervals).  Because the text holds both
// strands, the forward extension of a match R[a..e) is the backward extension of rc(R[a..e)) by the
// complement of R[e].
//
// The lookup table (katka-classify -F t; Appendix A, "Lookup tables"), as in Bowtie 2 and Cliffy, stores the
// interval of every string of t bases, so a search can start with one lookup instead of t steps.  When the
// string does not occur, the caller falls back to stepping base by base, so BML still learns exactly where
// the search fails, and the answers are the same with and without a table.
#pragma once
#include "rlcsa.hpp"

namespace katka {

struct interval { u64 sp, ep; };

struct searcher {
    const rlbwt *bwt = nullptr;
    const rlcsa *csa = nullptr;     // if set, steps use it instead of bwt

    interval full() const { return {0, bwt->n}; }
    // [sp, ep) for P -> [sp, ep) for cP; returns false, leaving I unchanged, if cP does not occur
    bool step(interval &I, unsigned char c) const {
        u64 a = I.sp, b = I.ep;
        if (csa) csa->extend(a, b, c); else bwt->extend(a, b, c);
        if (a >= b) return false;
        I = {a, b};
        return true;
    }
};

struct lookup_table {
    int t = 0;                      // 0: no table
    std::vector<interval> iv;       // indexed by the string's 2-bit code, first base most significant
    std::vector<uint8_t> occurs;

    // the strings of length i + 1 are those of length i with one more base in front
    void build(const searcher &S, int tt) {
        t = tt;
        std::vector<interval> a{S.full()};
        std::vector<uint8_t> ok{1};
        for (int i = 0; i < t; ++i) {
            u64 m = a.size();
            std::vector<interval> a2(4 * m);
            std::vector<uint8_t> ok2(4 * m, 0);
            for (int c = 0; c < 4; ++c)
                for (u64 w = 0; w < m; ++w)
                    if (ok[w]) { interval I = a[w]; if (S.step(I, "ACGT"[c])) { a2[c * m + w] = I; ok2[c * m + w] = 1; } }
            a.swap(a2); ok.swap(ok2);
        }
        iv.swap(a); occurs.swap(ok);
    }
    u64 bytes() const { return iv.size() * (sizeof(interval) + 1); }
    // the interval of p[0..t) (bases A/C/G/T only); false if it does not occur
    bool get(const unsigned char *p, interval &I) const {
        u64 key = 0;
        for (int i = 0; i < t; ++i) key = 4 * key + base_code(p[i]);
        if (!occurs[key]) return false;
        I = iv[key];
        return true;
    }
};

}  // namespace katka
