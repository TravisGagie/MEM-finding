// Finding the long MEMs of a read with Boyer-Moore-Li (BML; Section 5, and Appendix A "BML without
// digests").
//
// A MEM of R is a substring R[a..e) that occurs in the text but cannot be extended left or right and still
// occur.  BML finds every MEM of length at least L in a segment R[0..n):
//
//   x = 0
//   while x + L <= n:
//     search backward for R[x..x+L), from its right end
//     if the step for R[j-1] fails, no MEM of length >= L starts in [x, j), so set x = j and repeat
//     otherwise extend backward from x to the leftmost a such that R[a..x+L) occurs,
//       then forward from a to the largest e such that R[a..e) occurs: R[a..e) is a MEM
//     report it, and set x = e - L + 1 (a MEM of length >= L starting before that would end inside R[a..e))
//
// Since the text has both strands, the forward search is a backward search for rc(R[a..e)), and its
// interval, which we report, has the same tags as R[a..e)'s.
//
// Options (all from Section 5):
//   longest (-X)   adaptive L: after reporting a MEM longer than L, raise L to its length, so that the
//                  longest MEMs of the read are among those reported (with L = 1 at the start, the longest
//                  overall); shorter ones found before L was raised are reported too
//   filter  (-K)   KeBaB: search only the pseudo-MEMs of the read (kmer_filter.hpp) instead of its maximal
//                  A/C/G/T runs.  With longest, search them from the longest down and stop at the first
//                  one shorter than L: every MEM lies inside a pseudo-MEM, so no longer MEM can be missed
//   min3    (-D)   the low-complexity filter: ignore MEMs with fewer than min3 distinct 3-mers (they do not
//                  raise L either)
#pragma once
#include "kmer_filter.hpp"
#include "search.hpp"
#include <algorithm>

namespace katka {

struct bml_options {
    u64 L = 0;
    bool longest = false;
    int min3 = 0;
    const kmer_filter *filter = nullptr;
    const lookup_table *table = nullptr;
};

// number of distinct 3-mers of s[0..m)
inline int distinct_3mers(const char *s, u64 m) {
    uint64_t seen = 0;
    for (u64 i = 0; i + 3 <= m; ++i) {
        int c0 = base_code(s[i]), c1 = base_code(s[i + 1]), c2 = base_code(s[i + 2]);
        if (c0 >= 0 && c1 >= 0 && c2 >= 0) seen |= 1ULL << (c0 * 16 + c1 * 4 + c2);
    }
    return __builtin_popcountll(seen);
}

// The segments of the read that BML searches, as (start, length): its maximal A/C/G/T runs, or with the
// KeBaB filter the maximal substrings of those runs all of whose k-mers pass the filter.
inline void segments(const std::string &read, const kmer_filter *F, std::vector<std::pair<u64, u64>> &segs) {
    segs.clear();
    u64 i = 0, len = read.size();
    while (i < len) {
        while (i < len && base_code(read[i]) < 0) ++i;
        u64 r0 = i;
        while (i < len && base_code(read[i]) >= 0) ++i;
        u64 n = i - r0;
        if (n == 0) continue;
        if (!F) { segs.push_back({r0, n}); continue; }
        if (n < F->k) continue;
        static thread_local std::vector<uint64_t> kmers;
        kmers.assign(n - F->k + 1, 0);
        F->for_each_kmer((const unsigned char *)read.data() + r0, n, [&](uint64_t p, uint64_t c) { kmers[p] = c; });
        u64 start = 0; bool open = false;
        for (u64 p = 0; p < kmers.size(); ++p) {
            if (p + 8 < kmers.size()) F->prefetch(kmers[p + 8]);
            bool ok = F->has(kmers[p]);
            if (ok && !open) { start = p; open = true; }
            if (!ok && open) { segs.push_back({r0 + start, p + F->k - 1 - start}); open = false; }  // k-mers start..p-1
        }
        if (open) segs.push_back({r0 + start, n - start});
    }
}

// BML on one segment R[0..n); calls report(a, e, I) for every MEM R[a..e) found, with I the interval of
// rc(R[a..e)).  L may be raised (opt.longest).
template <class F>
void bml_segment(const searcher &S, const char *R, u64 n, u64 &L, const bml_options &opt, F report) {
    const lookup_table *T = opt.table && opt.table->t ? opt.table : nullptr;
    const u64 t = T ? T->t : 0;
    u64 x = 0;
    while (x + L <= n) {
        // backward search for R[x..x+L), starting with a lookup of its last t bases if possible
        interval I = S.full();
        u64 j = x + L;
        if (T && L >= t) { interval I2; if (T->get((const unsigned char *)R + j - t, I2)) { I = I2; j -= t; } }
        bool failed = false;
        while (j > x) {
            if (!S.step(I, R[j - 1])) { failed = true; break; }
            --j;
        }
        if (failed) { x = j; continue; }         // R[j-1..x+L) does not occur
        // extend left to the leftmost a, then right (on the reverse complement) to the largest e
        while (j > 0 && S.step(I, R[j - 1])) --j;
        u64 a = j, e = a;
        interval J = S.full();
        if (T && e + t <= n) {
            unsigned char rc[64];
            for (u64 q = 0; q < t; ++q) rc[q] = complement(R[e + t - 1 - q]);
            interval J2; if (T->get(rc, J2)) { J = J2; e += t; }
        }
        while (e < n && S.step(J, complement(R[e]))) ++e;
        // the MEM R[a..e)
        if (opt.min3 && distinct_3mers(R + a, e - a) < opt.min3) { x = e - L + 1; continue; }
        report(a, e, J);
        if (opt.longest && e - a > L) L = e - a;  // no longer MEM can start at or before a
        x = e - L + 1;
    }
}

// BML on a whole read; calls report(a, e, I) as above, with positions in the read
template <class F>
void find_mems(const searcher &S, const std::string &read, const bml_options &opt, F report) {
    static thread_local std::vector<std::pair<u64, u64>> segs;
    segments(read, opt.filter, segs);
    const bool longest_first = opt.longest && opt.filter;
    if (longest_first)
        std::stable_sort(segs.begin(), segs.end(), [](const std::pair<u64, u64> &p, const std::pair<u64, u64> &q) { return p.second > q.second; });
    u64 L = opt.L;
    for (auto [r0, n] : segs) {
        if (n < L) { if (longest_first) break; continue; }
        bml_segment(S, read.data() + r0, n, L, opt, [&](u64 a, u64 e, const interval &I) { report(r0 + a, r0 + e, I); });
    }
}

}  // namespace katka
