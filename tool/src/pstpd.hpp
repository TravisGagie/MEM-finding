// STPD over the parse (NOTES.md, "STPD over the parse"): the structures on top of the parse oracle
// (parse_text.hpp: phrase ids Z[0..m), boundaries, dictionary strings) and the two-level index (B, the
// dictionary's hash table).  The parse is Z$ = Z[0..m) followed by 0; parse row 0 is the suffix 0.
//   pos     the parse sample: parse positions SA_Z[x] - 1 for the needed run starts x of BWT_Z, sorted by
//           parse suffix (the right-to-left STPD sample of Z$, as stpd.hpp's for S$), pruned with the
//           character-level LCP of the suffixes of S$ at the boundaries (see below)
//   entry   entry[y] = SA_Z[C[y]]: the parse position of the smallest parse suffix starting with id y
//   C       C[y] = the first parse row of id y (C[0] = 0, the row of the terminator; C[sigma+1] = m+1)
//   phi^-1  over the parse (SA_Z at the ends of the parse BWT's runs -> SA_Z at the next run's start)
//
// Pruning.  The character-level rule (stpd.hpp) keeps a run start x of BWT_Z iff lcp(SA[y], SA[x]) <
// lcp(SA[x-1], SA[x]), y the last row before x with the same BWT symbol.  Over the parse, the ranges a
// search meets are the parse rows whose suffixes, expanded to characters, are prefixed by a string
// (whole phrase ids, then a prefix of a phrase), not ranges of id strings; such a range contains x-1, x
// and not y only if the character-level LCP of SA[x-1] and SA[x] exceeds that of SA[y] and SA[x].  So
// the LCPs are measured in characters of the expansions: the shared phrases contribute len - k each (the
// next starts k characters before the end of the previous one), then the LCP of the first two different
// phrases' strings (the dictionary is prefix-free, so it ends inside both).
#pragma once
#include "parse_text.hpp"
#include "libsais.h"

namespace katka {

struct parse_stpd {
    u64 m = 0, runs = 0;
    uint32_t sigma = 0, K = 0, SM = 0;
    PackedInts pos, entry;
    EliasFano C;
    EliasFano ikeys; PackedInts ivals; u64 last_sa = 0;

    inline u64 c_of(u64 y) const { return C.access(y); }
    inline u64 id_of_row(u64 rho) const { return C.count_less(rho + 1) - 1; }
    inline u64 phi_inv(u64 v) const { u64 t = ikeys.count_less(v + 1) - 1; return ivals.get(t) + (v - ikeys.access(t)); }
    static inline u64 zat(const parse_text &PT, u64 v) { return v < PT.m ? PT.ids.get(v) : 0; }

    // Compares the parse suffix Z[b..] with the pattern element string: ids e[0..h), then the characters
    // xr[0..lr) as a prefix of the next phrase's expansion.  Returns -1, 0 (Z[b..] is in the range: its
    // expansion is prefixed by the pattern's) or +1, and the number of leading ids that matched (lcp).
    static int cmp(const parse_text &PT, u64 b, const uint32_t *e, u64 h, const uint8_t *xr, u64 lr, u64 skip, u64 &lcp) {
        for (u64 g = skip; g < h; ++g) {
            u64 z = zat(PT, b + g);
            if (z != e[g]) { lcp = g; return z < e[g] ? -1 : 1; }
        }
        lcp = h;
        u64 z = zat(PT, b + h);
        if (z == 0) return -1;                                   // the terminator is smallest
        const uint8_t *p = PT.dict.data() + PT.doff.get(z - 1);
        u64 lp = PT.doff.get(z) - PT.doff.get(z - 1), q = std::min(lp, lr);
        u64 l = lce_bytes(p, xr, q);
        if (l < q) return p[l] < xr[l] ? -1 : 1;
        if (l == lr) return 0;
        // the phrase is a proper prefix of xr: compare the expansion through the oracle (not expected:
        // a phrase ends with a syncmer, and xr holds no syncmer after its first position)
        u64 s = PT.starts.access(b + h), l2 = PT.lce(s, xr, std::min(lr, PT.n - s));
        if (l2 == lr) return 0;
        return (s + l2 == PT.n || PT.at(s + l2) < xr[l2]) ? -1 : 1;
    }
    // The sampled parse position of the smallest parse suffix in the range of (e, xr), or UINT64_MAX.
    u64 find(const parse_text &PT, const uint32_t *e, u64 h, const uint8_t *xr, u64 lr, u64 &compares) const {
        u64 lo = 0, hi = pos.size(), llo = 0, lhi = 0;
        while (lo < hi) {
            u64 mid = (lo + hi) >> 1, l;
            ++compares;
            int r = cmp(PT, pos.get(mid), e, h, xr, lr, std::min(llo, lhi), l);
            if (r < 0) { lo = mid + 1; llo = l; } else { hi = mid; lhi = l; }
        }
        if (lo == pos.size()) return UINT64_MAX;
        u64 b = pos.get(lo), l;
        ++compares;
        return cmp(PT, b, e, h, xr, lr, 0, l) == 0 ? b : UINT64_MAX;
    }

    // From the parse oracle (the same parse as the two-level index).  Memory: the parse as 32-bit ids and
    // its suffix array (8 bytes per phrase), plus O(r_Z).
    void build(const parse_text &PT, bool prune = true) {
        double t0 = now();
        m = PT.m; sigma = PT.sigma; K = PT.K; SM = PT.SM;
        if (m + 1 >= (1ull << 31)) die("too many phrases for a 32-bit parse suffix array");
        std::vector<uint32_t> Z(m + 1), SA(m + 1);
        for (u64 i = 0; i < m; ++i) Z[i] = (uint32_t)PT.ids.get(i);
        Z[m] = 0;
        if (libsais_int(reinterpret_cast<int32_t *>(Z.data()), reinterpret_cast<int32_t *>(SA.data()), (int32_t)(m + 1), (int32_t)sigma + 1, 0) != 0)
            die("libsais_int failed");
        for (u64 i = 0; i < m; ++i) Z[i] = (uint32_t)PT.ids.get(i);      // (libsais_int overwrote it)
        Z[m] = 0;
        const u64 N = m + 1;
        auto bwt = [&](u64 x) -> uint32_t { return SA[x] == 0 ? 0 : Z[SA[x] - 1]; };
        // C
        {
            std::vector<uint64_t> cnt(sigma + 2, 0);
            for (u64 i = 0; i < N; ++i) cnt[Z[i] + 1]++;
            for (u64 y = 1; y <= sigma + 1; ++y) cnt[y] += cnt[y - 1];
            for (u64 y = 0; y <= sigma; ++y) if (cnt[y + 1] == cnt[y]) die("parse STPD: an id does not occur in the parse");
            C.build(cnt, N + 1);
            std::vector<uint64_t> ent(sigma + 1, 0);
            for (u64 y = 1; y <= sigma; ++y) ent[y] = SA[cnt[y]];
            entry.build(ent);
        }
        auto len_of = [&](u64 id) { return PT.doff.get(id) - PT.doff.get(id - 1); };
        auto clcp = [&](u64 a, u64 b) -> u64 {               // character LCP of the boundary suffixes at parse positions a, b
            u64 acc = 0;
            while (Z[a] == Z[b] && Z[a] != 0) { acc += len_of(Z[a]) - K; ++a; ++b; }
            if (Z[a] == 0 || Z[b] == 0) return acc;
            const uint8_t *pa = PT.dict.data() + PT.doff.get(Z[a] - 1), *pb = PT.dict.data() + PT.doff.get(Z[b] - 1);
            return acc + lce_bytes(pa, pb, std::min(len_of(Z[a]), len_of(Z[b])));
        };
        // runs, SA at their ends, phi^-1, and the sample (rows LF(x) from running counts)
        std::vector<uint64_t> rank(sigma + 1, 0);
        std::vector<int64_t> prev_end(sigma + 1, -1);            // SA at the last row so far of each symbol's last run
        std::vector<std::pair<u64, u64>> kv, smp;                // (esa, next ssa); (row, position)
        u64 cand = 0;
        runs = 0;
        u64 prev_sa_end = 0;
        for (u64 x = 0; x < N; ++x) {
            uint32_t c = bwt(x);
            bool start = x == 0 || c != bwt(x - 1);
            if (start) {
                ++runs;
                if (x > 0) kv.push_back({SA[x - 1], SA[x]});
                if (x > 0 && SA[x] != 0) {
                    ++cand;
                    bool keep = true;
                    if (prune && prev_end[c] >= 0) keep = clcp((u64)prev_end[c], SA[x]) < clcp(SA[x - 1], SA[x]);
                    if (keep) smp.push_back({C.access(c) + rank[c], SA[x] - 1});
                }
            }
            ++rank[c];
            bool end = x + 1 == N || bwt(x + 1) != c;
            if (end) prev_end[c] = SA[x];
            prev_sa_end = SA[x];
        }
        last_sa = prev_sa_end;
        std::vector<uint32_t>().swap(Z);
        std::vector<uint32_t>().swap(SA);
        std::sort(kv.begin(), kv.end());
        {
            std::vector<uint64_t> kk(kv.size()), vv(kv.size());
            for (size_t q = 0; q < kv.size(); ++q) { kk[q] = kv[q].first; vv[q] = kv[q].second; }
            ikeys.build(kk, N); ivals.build(vv);
        }
        std::sort(smp.begin(), smp.end());
        std::vector<uint64_t> p(smp.size());
        for (size_t q = 0; q < smp.size(); ++q) p[q] = smp[q].second;
        pos.build(p);
        fprintf(stderr, "parse STPD: m=%lu, sigma=%u, r_Z=%lu parse BWT runs, %lu run starts, %lu sampled (%.3f of r_Z): sample %.3f MB, entry %.3f MB, C %.3f MB, phi^-1 %.3f MB (%.1f s)\n",
                m, sigma, runs, cand, (u64)p.size(), (double)p.size() / runs, pos.bytes() / 1e6, entry.bytes() / 1e6, C.bytes() / 1e6,
                (ikeys.bytes() + ivals.bytes()) / 1e6, now() - t0);
    }
    u64 bytes_locate() const { return ikeys.bytes() + ivals.bytes(); }
    void save(const std::string &f) const {
        FILE *o = open_or_die(f, "wb"); wmagic(o, "KTPSTP01");
        wr(o, m); wr(o, runs); wr(o, sigma); wr(o, K); wr(o, SM);
        pos.save(o); entry.save(o); C.save(o); ikeys.save(o); ivals.save(o); wr(o, last_sa);
        fclose(o);
    }
    void load(const std::string &f) {
        FILE *i = open_or_die(f, "rb"); rmagic(i, "KTPSTP01", f);
        rd(i, m); rd(i, runs); rd(i, sigma); rd(i, K); rd(i, SM);
        pos.load(i); entry.load(i); C.load(i); ikeys.load(i); ivals.load(i); rd(i, last_sa);
        fclose(i);
    }
};

}  // namespace katka
