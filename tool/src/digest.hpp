// Symmetric minimizer digests (App. B of the paper, "Minimizer digests"; App. A, "Digests"), from KATKA's
// rz-index/digest.hpp.
//
// Number A, C, G, T as 0..3; a k-mer's code is its 2k-bit integer and its hash is the splitmix64 finalizer
// of the smaller of its code and its reverse complement's (a canonical hash).  In each maximal A/C/G/T run,
// every window of w bases (w - k + 1 consecutive k-mers) selects ALL its k-mers of minimum hash; the
// selected k-mers, in text order and in their own orientation, form the digest of the run, after
// collapsing equal consecutive codes.  Digests of consecutive runs are joined by X.  Ties are broken
// symmetrically, so digest(rc(G)) = rc(digest(G)), where rc reverses and maps each symbol to that of the
// reverse-complemented k-mer: the digested text stays closed under reverse complement.
//
// Symbols are bytes: 0..2 are terminators and X separates, so 252 bytes are free; with k = 4 all 256 codes
// occur, and the rarest codes, taken in reverse-complement pairs until the rest fit (6 codes), share one
// byte.  The others get bytes 3, 4, ... (skipping X) in increasing order of code.  The map is
// PREFIX.map (KATKA's format, "RZDGMAP1").
#pragma once
#include "common.hpp"
#include <algorithm>
#include <deque>

namespace katka {

inline uint64_t dg_mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline uint32_t dg_rc_code(uint32_t f, int k) {
    uint32_t r = 0;
    for (int i = 0; i < k; ++i) { r = (r << 2) | (3 - (f & 3)); f >>= 2; }
    return r;
}

// The selected k-mers of one A/C/G/T run s[0..n), n >= k: calls sel(p) for each selected k-mer start, in
// order, and, if win is given, win(y, first, last) for each full window y with the first and last selected
// k-mer starts in it.  A run shorter than w is one window over all its k-mers.
template <class S, class Wn>
inline void dg_select(const char *s, uint64_t n, int k, int w, std::vector<uint32_t> &f, S sel, Wn win) {
    uint64_t nk = n - k + 1, W = (uint64_t)(w - k + 1);
    if (W > nk) W = nk;
    static thread_local std::vector<uint64_t> h;
    static thread_local std::vector<char> mark;
    h.resize(nk); f.resize(nk); mark.assign(nk, 0);
    uint32_t mask = (uint32_t)((1ull << (2 * k)) - 1), fw = 0;
    for (uint64_t i = 0; i < n; ++i) {
        fw = ((fw << 2) | (uint32_t)base_code(s[i])) & mask;
        if (i + 1 >= (uint64_t)k) {
            uint64_t p = i + 1 - k;
            uint32_t rc = dg_rc_code(fw, k);
            f[p] = fw; h[p] = dg_mix(fw < rc ? fw : rc);
        }
    }
    std::deque<uint64_t> dq;
    for (uint64_t p = 0; p < nk; ++p) {
        while (!dq.empty() && h[dq.back()] > h[p]) dq.pop_back();
        dq.push_back(p);
        if (p + 1 >= W) {
            uint64_t y = p + 1 - W;
            while (dq.front() < y) dq.pop_front();
            uint64_t mn = h[dq.front()], last = dq.front();
            for (uint64_t q : dq) { if (h[q] != mn) break; mark[q] = 1; last = q; }
            win(y, dq.front(), last);
        }
    }
    for (uint64_t p = 0; p < nk; ++p) if (mark[p]) sel(p);
}

struct digest_map {
    int k = 0, w = 0;
    unsigned char code2byte[257];      // codes 0..4^k-1, and 4^k for the separator; 0 = never seen
    unsigned char comp[256];           // byte -> byte of the reverse-complemented k-mer (others: themselves)
    bool on() const { return k > 0; }

    // the digest of a genome (non-ACGT characters separate runs) as codes, runs joined by 4^k
    void codes(const std::string &g, std::vector<uint32_t> &out) const {
        out.clear();
        const uint32_t sep = 1u << (2 * k);
        static thread_local std::vector<uint32_t> f;
        uint64_t i = 0, n = g.size();
        while (i < n) {
            while (i < n && base_code(g[i]) < 0) ++i;
            uint64_t j = i;
            while (j < n && base_code(g[j]) >= 0) ++j;
            if (j - i >= (uint64_t)k) {
                bool first = true;
                dg_select(g.data() + i, j - i, k, w, f, [&](uint64_t p) {
                    if (first) { if (!out.empty()) out.push_back(sep); first = false; }
                    if (out.empty() || out.back() != f[p]) out.push_back(f[p]);
                }, [](uint64_t, uint64_t, uint64_t) {});
            }
            i = j;
        }
    }
    std::string digest(const std::string &g) const {
        std::vector<uint32_t> c;
        codes(g, c);
        std::string out(c.size(), 0);
        for (size_t i = 0; i < c.size(); ++i) out[i] = (char)code2byte[c[i]];
        return out;
    }
    std::string rc(const std::string &s) const {
        std::string r(s.rbegin(), s.rend());
        for (auto &c : r) c = (char)comp[(unsigned char)c];
        return r;
    }
    // pass 1 over the genomes: which codes occur (on either strand); give each a byte
    template <class ForEachGenome> void build(int kk, int ww, ForEachGenome each) {
        k = kk; w = ww;
        const uint32_t K = 1u << (2 * k);
        std::vector<uint64_t> seen(K + 1, 0);
        std::vector<uint32_t> c;
        each([&](const std::string &g) {
            codes(g, c);
            for (uint32_t x : c) if (x < K) { seen[x]++; seen[dg_rc_code(x, k)]++; }
        });
        std::fill(code2byte, code2byte + 257, 0);
        for (int i = 0; i < 256; ++i) comp[i] = (unsigned char)i;
        std::vector<uint32_t> occ;
        uint64_t total = 0;
        for (uint32_t x = 0; x < K; ++x) if (seen[x]) { occ.push_back(x); total += seen[x]; }
        std::vector<char> merged(K, 0);
        uint64_t nmerged = 0, merged_occ = 0;
        if (occ.size() > 252) {
            std::vector<uint32_t> byf(occ);
            std::sort(byf.begin(), byf.end(), [&](uint32_t a, uint32_t b) { return seen[a] < seen[b]; });
            for (uint32_t x : byf) {
                if (occ.size() - nmerged + 1 <= 252) break;
                for (uint32_t y : {x, dg_rc_code(x, k)}) if (!merged[y]) { merged[y] = 1; ++nmerged; merged_occ += seen[y]; }
            }
        }
        int next = 3, shared = 0;
        for (uint32_t x : occ) {
            if (merged[x]) continue;
            if (next == 'X') ++next;
            code2byte[x] = (unsigned char)next++;
        }
        if (nmerged) { if (next == 'X') ++next; shared = next; for (uint32_t x : occ) if (merged[x]) code2byte[x] = (unsigned char)shared; }
        if (next > 255) die("digest: too many distinct minimizers");
        code2byte[K] = 'X';
        for (uint32_t x = 0; x < K; ++x) if (code2byte[x]) comp[code2byte[x]] = code2byte[dg_rc_code(x, k)];
        fprintf(stderr, "digest k=%d w=%d: %zu distinct minimizers", k, w, occ.size());
        if (nmerged) fprintf(stderr, "; the %lu rarest share byte %d (%.3f%% of occurrences)", nmerged, shared, 100.0 * merged_occ / std::max<uint64_t>(1, total));
        fprintf(stderr, "\n");
    }
    void save(const std::string &f) const {
        FILE *o = open_or_die(f, "wb"); wmagic(o, "RZDGMAP1"); wr(o, k); wr(o, w);
        fwrite(code2byte, 1, 257, o); fwrite(comp, 1, 256, o); fclose(o);
    }
    void load(const std::string &f) {
        FILE *i = open_or_die(f, "rb"); rmagic(i, "RZDGMAP1", f); rd(i, k); rd(i, w);
        if (fread(code2byte, 1, 257, i) != 257 || fread(comp, 1, 256, i) != 256) die("cannot read " + f);
        fclose(i);
    }
};

// The digest of one A/C/G/T run R[0..n) of a read, n >= w, with, for each window y (bases y..y+w-1), the
// positions lo[y] and hi[y] in the digest of the first and last symbols it selects (both non-decreasing).
struct read_digest {
    std::string dg;
    std::vector<uint64_t> lo, hi;
    void build(const digest_map &M, const char *R, uint64_t n) {
        static thread_local std::vector<uint32_t> f;
        static thread_local std::vector<uint64_t> mlo, mhi, idx;
        const uint64_t nk = n - M.k + 1, nwin = nk - (M.w - M.k + 1) + 1;
        mlo.resize(nwin); mhi.resize(nwin); idx.assign(nk, 0);
        dg.clear();
        uint32_t lastcode = UINT32_MAX;
        dg_select(R, n, M.k, M.w, f, [&](uint64_t p) {
            if (dg.empty() || lastcode != f[p]) { dg.push_back((char)M.code2byte[f[p]]); lastcode = f[p]; }
            idx[p] = dg.size() - 1;
        }, [&](uint64_t y, uint64_t a, uint64_t b) { mlo[y] = a; mhi[y] = b; });
        lo.resize(nwin); hi.resize(nwin);
        for (uint64_t y = 0; y < nwin; ++y) { lo[y] = idx[mlo[y]]; hi[y] = idx[mhi[y]]; }
    }
    // the span in bases of the digested match dg[a..b]: from the first window with lo >= a to just before
    // the end of the first window with hi > b (to the end of the run if none)
    void span(uint64_t a, uint64_t b, uint64_t n, int w, uint64_t &left, uint64_t &right, uint64_t &yR) const {
        left = std::lower_bound(lo.begin(), lo.end(), a) - lo.begin();
        yR = std::upper_bound(hi.begin(), hi.end(), b) - hi.begin();
        right = yR >= lo.size() ? n - 1 : yR + w - 2;
    }
};

}  // namespace katka
