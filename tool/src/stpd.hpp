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
#include "parse_text.hpp"

namespace katka {

struct stpd_index {
    u64 n = 0;                        // |S$|
    std::vector<uint8_t> S;           // the text S$ (the oracle), with $ = byte 1
    const parse_text *PT = nullptr;   // or, if set (and S empty), S$ read through the two-level parse
    // the text oracle: S$[x], and the longest common prefix of S$[b..b+len) and x[0..len)
    inline uint8_t at(u64 x) const { return PT ? PT->at(x) : S[x]; }
    inline u64 lce_at(u64 b, const uint8_t *x, u64 len) const { return PT ? PT->lce(b, x, len) : lce_bytes(S.data() + b, x, len); }
    PackedInts pos;                   // the sample, sorted lexicographically by suffix
    u64 runs = 0;                     // r, for the statistics
    // phi^-1 (SA[x] -> SA[x+1]), for locating every occurrence: the SA values at the ends of the BWT runs
    // (sorted), the SA value at the start of the next run for each, and SA[n-1] (no successor)
    EliasFano ikeys; PackedInts ivals; u64 last_sa = 0;

    // phi^-1(i) = SA[ISA[i] + 1] (i != SA[n-1]): from the run end with the largest SA value <= i, as the
    // r-index computes phi (Gagie, Navarro and Prezza)
    inline u64 phi_inv(u64 i) const { u64 t = ikeys.count_less(i + 1) - 1; return ivals[t] + (i - ikeys.access(t)); }

    // Every occurrence of x[0..len), given its first in SA order, t = SA[sp]: phi^-1 steps in batches of
    // doubling size, checking against the text only the last position of each batch (the occurrences are a
    // prefix of the sequence SA[sp], SA[sp+1], ...), then a binary search in the batch that overshot, as in
    // Becker et al.'s locating with st-colex.  Appends to out; returns the number of occurrences.
    u64 locate(const uint8_t *x, u64 len, u64 t, std::vector<u64> &out, u64 &checks) const {
        const u64 first = out.size();
        out.push_back(t);
        auto occurs_at = [&](u64 p) { ++checks; return p + len <= n - 1 && lce_at(p, x, len) >= len; };
        u64 batch = 1;
        for (;;) {
            const u64 b0 = out.size();
            bool end = false;
            for (u64 q = 0; q < batch; ++q) {
                if (out.back() == last_sa) { end = true; break; }
                out.push_back(phi_inv(out.back()));
            }
            if (out.size() == b0) break;                        // the last row was reached
            if (!occurs_at(out.back())) {                       // the boundary is in this batch: binary search
                u64 lo = b0, hi = out.size() - 1;               // out[lo..hi]: hi fails; find the first failing
                while (lo < hi) { u64 mid = (lo + hi) / 2; if (occurs_at(out[mid])) lo = mid + 1; else hi = mid; }
                out.resize(lo);
                break;
            }
            if (end) break;
            batch *= 2;
        }
        return out.size() - first;
    }

    // compares the suffix S[b..] with x[0..len): returns the length of their common prefix and, if it is
    // shorter than len, whether the suffix is smaller
    inline u64 cmp(u64 b, const uint8_t *x, u64 len, bool &less) const {
        u64 m = std::min<u64>(len, n - b);
        u64 l = lce_at(b, x, m);
        less = l < len && (l == n - b || at(b + l) < x[l]);
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
        {   // phi^-1: SA at the end of run k -> SA at the start of run k + 1
            std::vector<std::pair<u64, u64>> kv(R - 1);
            for (u64 k = 0; k + 1 < R; ++k) kv[k] = {esa[k], ssa[k + 1]};
            std::sort(kv.begin(), kv.end());
            std::vector<u64> kk(R - 1), vv(R - 1);
            for (u64 k = 0; k + 1 < R; ++k) { kk[k] = kv[k].first; vv[k] = kv[k].second; }
            ikeys.build(kk, n); ivals.build(vv); last_sa = esa[R - 1];
        }
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
        fprintf(stderr, "STPD: n=%lu r=%lu, %lu run starts, %lu sampled positions (%.3f of r), %.3f MB on top of the text (+%.3f MB for phi^-1), %.1f s\n",
                n, R, cand, kept, (double)kept / R, pos.bytes() / 1e6, locate_bytes() / 1e6, now() - t0);
    }
    u64 bytes() const { return pos.bytes(); }
    u64 locate_bytes() const { return ikeys.bytes() + ivals.bytes(); }
    void save(const std::string &f) const {
        FILE *o = open_or_die(f, "wb"); wmagic(o, "KTSTPD02"); wr(o, n); wr(o, runs); pos.save(o);
        ikeys.save(o); ivals.save(o); wr(o, last_sa); fclose(o);
    }
    // with Sfile empty, the text is not loaded (set PT before searching)
    void load(const std::string &f, const std::string &Sfile) {
        FILE *i = open_or_die(f, "rb"); rmagic(i, "KTSTPD02", f); rd(i, n); rd(i, runs); pos.load(i);
        ikeys.load(i); ivals.load(i); rd(i, last_sa); fclose(i);
        if (Sfile.empty()) return;
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
        if (st.t > 0 && D.at(st.t - 1) == c) { --st.t; --st.j; continue; }
        ++cnt.stpd_jumps;
        u64 b = D.find(x + st.j - 1, st.end - st.j + 1, cnt.stpd_compares);
        if (b == UINT64_MAX) break;
        st.t = b; --st.j;
    }
}

}  // namespace katka
