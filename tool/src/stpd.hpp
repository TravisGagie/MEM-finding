// A suffix-tree path decomposition (STPD) backend: longest-suffix searches by binary searches over a
// sample of text positions, with the text itself as the oracle (Becker, Cenzato, Gagie, Groot Koerkamp,
// Kim, Manzini and Prezza, "Compressing Suffix Trees by Path Decompositions", ICALP 2026).
//
// The paper matches a pattern left to right with an order-preserving permutation pi, keeping the
// occurrence of the matched prefix that minimizes pi.  BML and forward-backward search only need the
// longest suffix of X[lo..end) that occurs (the forward searches are on the reverse complement, since S
// holds both strands), so we use the mirror image: match right to left, keeping the occurrence of the
// matched suffix Z = X[j..end) whose suffix of S is lexicographically smallest, i.e. t = SA[sp(Z)].  This is
// the paper's st-colex on the reversed text: pi = ISA, which is order-preserving for leftward extension
// because LF keeps the order of rows with the same BWT character.
//
//   start: Z empty, t = SA[0] = n (the suffix $)
//   step with c = X[j-1]:  if S[t-1] = c, then t - 1 starts the smallest suffix prefixed by cZ (LF maps
//     the first row of Z's interval to the first row of cZ's); otherwise the smallest suffix prefixed by
//     cZ starts at a sampled position, found by a binary search over the sample, sorted lexicographically
//     by suffix, comparing S[b..] with cZ = X[j-1..end) character by character; if no sampled suffix is
//     prefixed by cZ, then cZ does not occur and the search stops.
//
// The sample.  If BWT[sp(Z)] != c but cZ occurs, the first row x of Z's interval with BWT[x] = c starts a
// BWT run (BWT[x-1] != c, and x - 1 is in the interval), and SA[x] - 1 starts the smallest suffix prefixed
// by cZ.  Not every run start is needed: with y the last row before x with BWT[y] = c, x is the first c in
// some interval of a prefix of its suffix iff lcp(SA[y], SA[x]) < lcp(SA[x-1], SA[x]) (or there is no y).
// The sample is { SA[x] - 1 : x > 0 a needed run start }, sorted by the row of SA[x] - 1, which is LF(x):
// at most r positions, which we compute from the LF pass of the r-index (SA at the ends of the runs) and
// longest-common-extension queries on the text.  The pruned sample corresponds to the irreducible
// positions of the paper; the full set of run starts is also correct (--stpd-all).
//
// An STPD search returns positions, not BWT intervals: the interval of a MEM, when it is needed (tags,
// the mems format), is computed by an ordinary backward search over the MEM.
#pragma once
#include "tag_build.hpp"
#include "twolevel.hpp"

namespace katka {

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

struct stpd_index {
    u64 n = 0;                        // |S$|
    std::vector<uint8_t> S;           // the text S$ (the oracle), with $ = byte 1
    PackedInts pos;                   // the sample, sorted lexicographically by suffix
    u64 runs = 0;                     // r, for the statistics

    // compares the suffix S[b..] with x[0..len): returns the length of their common prefix and, if it is
    // shorter than len, whether the suffix is smaller
    inline u64 cmp(u64 b, const uint8_t *x, u64 len, bool &less) const {
        u64 m = std::min<u64>(len, n - b);
        u64 l = lce_bytes(S.data() + b, x, m);
        less = l < len && (l == n - b || S[b + l] < x[l]);
        return l;
    }
    // the sampled position starting the smallest suffix prefixed by x[0..len), or UINT64_MAX
    u64 find(const uint8_t *x, u64 len, u64 &compares) const {
        u64 lo = 0, hi = pos.size(), llo = 0, lhi = 0;       // common prefixes with the bounds
        while (lo < hi) {
            u64 mid = (lo + hi) >> 1, skip = std::min(llo, lhi), b = pos[mid];
            bool less;
            ++compares;
            u64 l = skip + cmp(b + skip, x + skip, len - skip, less);
            if (l >= len) { hi = mid; lhi = l; }                // prefixed by x: go left
            else if (less) { lo = mid + 1; llo = l; }
            else { hi = mid; lhi = l; }
        }
        if (lo == pos.size()) return UINT64_MAX;
        u64 b = pos[lo];
        bool less;
        return cmp(b, x, len, less) >= len ? b : UINT64_MAX;
    }

    // build from the BWT file of S$ (one terminator) and the text
    void build(const std::string &Sfile, const std::string &bwtfile, const rlbwt &B, bool prune) {
        if (B.C[2] - B.C[1] != 1) die("the STPD index needs a BWT with one terminator (BWT(S$))");
        double t0 = now();
        n = B.n; runs = B.R;
        S.resize(n);
        { FILE *f = open_or_die(Sfile, "rb"); if (fread(S.data(), 1, n - 1, f) != n - 1) die("cannot read " + Sfile); fclose(f); }
        S[n - 1] = 1;
        if (file_size(Sfile) != n - 1) die(Sfile + " does not match the BWT");
        // SA at the first and last row of every run, by one pass of LF from the end of S$
        move_lf M; M.build(bwtfile, B);
        const u64 R = B.R;
        std::vector<u64> ssa(R, UINT64_MAX), esa(R, UINT64_MAX);
        {
            u64 row = 0, k = M.run_containing(0);              // row 0 is the suffix $, position n - 1
            for (u64 p = n - 1;; --p) {
                if (row == M.runs[k].start()) ssa[k] = p;
                if (row + 1 == M.runs[k + 1].start()) esa[k] = p;
                if (p == 0) break;
                M.LF(row, k);
            }
        }
        for (u64 k = 0; k < R; ++k) if (ssa[k] == UINT64_MAX || esa[k] == UINT64_MAX) die("internal error: a run without samples");
        // the needed run starts, with their rows LF(x)
        std::vector<std::pair<u64, u64>> v;                     // (row of SA[x] - 1, SA[x] - 1)
        std::vector<int64_t> prev_of(256, -1);                  // the last run so far of each character
        u64 kept = 0, cand = 0;
        auto lce_pos = [&](u64 a, u64 b) { return lce_bytes(S.data() + a, S.data() + b, n - std::max(a, b)); };
        for (u64 k = 0; k < R; ++k) {
            unsigned char c = M.head[k];
            int64_t y = prev_of[c];
            prev_of[c] = (int64_t)k;
            if (k == 0 || ssa[k] == 0) continue;               // row 0, or the row of position 0 (BWT = $)
            ++cand;
            if (prune && y >= 0) {
                u64 lx = lce_pos(esa[k - 1], ssa[k]), ly = lce_pos(esa[y], ssa[k]);
                if (ly >= lx) continue;
            }
            v.push_back({M.runs[k].lf_start(), ssa[k] - 1});
            ++kept;
        }
        std::vector<move_lf::run>().swap(M.runs);
        std::sort(v.begin(), v.end());
        std::vector<u64> p(v.size());
        for (u64 i = 0; i < v.size(); ++i) p[i] = v[i].second;
        pos.build(p);
        fprintf(stderr, "STPD: n=%lu r=%lu, %lu run starts, %lu sampled positions (%.3f of r), %.3f MB on top of the text, %.1f s\n",
                n, R, cand, kept, (double)kept / R, pos.bytes() / 1e6, now() - t0);
    }
    u64 bytes() const { return pos.bytes(); }
    void save(const std::string &f) const { FILE *o = open_or_die(f, "wb"); wmagic(o, "KTSTPD01"); wr(o, n); wr(o, runs); pos.save(o); fclose(o); }
    void load(const std::string &f, const std::string &Sfile) {
        FILE *i = open_or_die(f, "rb"); rmagic(i, "KTSTPD01", f); rd(i, n); rd(i, runs); pos.load(i); fclose(i);
        if (file_size(Sfile) != n - 1) die(Sfile + " does not match " + f);
        S.resize(n);
        FILE *s = open_or_die(Sfile, "rb");
        if (fread(S.data(), 1, n - 1, s) != n - 1) die("cannot read " + Sfile);
        fclose(s);
        S[n - 1] = 1;
    }
};

// The state of an STPD search of X[j..end): t starts the smallest suffix of S prefixed by it.
struct stpd_state { u64 end = 0, j = 0, t = 0; int dir = -1; };

// extends st (X[st.j..st.end) occurring at st.t) to the left, down to lo, as far as it occurs
inline void stpd_extend(const stpd_index &D, const char *X, u64 lo, stpd_state &st, search_counters &cnt) {
    const uint8_t *x = (const uint8_t *)X;
    while (st.j > lo) {
        unsigned char c = x[st.j - 1];
        ++cnt.char_steps;
        if (st.t > 0 && D.S[st.t - 1] == c) { --st.t; --st.j; continue; }
        ++cnt.stpd_jumps;
        u64 b = D.find(x + st.j - 1, st.end - st.j + 1, cnt.stpd_compares);
        if (b == UINT64_MAX) break;
        st.t = b; --st.j;
    }
}

}  // namespace katka
