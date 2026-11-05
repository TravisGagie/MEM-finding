// Building twolevel_locate (twolevel.hpp): the samples for the toehold of two-level search with blocks, and
// phi^-1, from the BWT file, the RLBWT, the two-level index (its parse BWT and block BWT, B and BS) and the
// parse oracle (the same parse: its ids, for a suffix array of the parse, and its phrase starts).
#pragma once
#include "tag_build.hpp"
#include "parse_text.hpp"
#include "libsais.h"

namespace katka {

inline void build_twolevel_locate(const std::string &bwtfile, const rlbwt &B, const two_level &L, const parse_text &PT, twolevel_locate &T) {
    if (L.maxblock) die("two-level locating is not supported with --max-block (a cut breaks the toehold's entry into the blocks)");
    double t0 = now();
    if (!L.has_blocks) die("two-level locating needs the blocks (an index without --no-blocks)");
    if (PT.m != L.m || PT.n != L.n || (int32_t)PT.K != L.K || (int32_t)PT.SM != L.SM) die("the parse oracle does not match the two-level index");
    T.n = B.n; T.m = L.m; T.nblocks = L.nblocks;
    // ---- characters: SA at the start and end of every run, by one LF pass (as stpd_index::build)
    {
        move_lf M; M.build(bwtfile, B);
        const u64 R = B.R, n = B.n;
        std::vector<u64> ssa(R, UINT64_MAX), esa(R, UINT64_MAX);
        {
            u64 row = 0, k = M.run_containing(0);
            for (u64 p = n - 1;; --p) {
                if (row == M.runs[k].start()) ssa[k] = p;
                if (row + 1 == M.runs[k + 1].start()) esa[k] = p;
                if (p == 0) break;
                M.LF(row, k);
            }
        }
        for (u64 k = 0; k < R; ++k) if (ssa[k] == UINT64_MAX || esa[k] == UINT64_MAX) die("internal error: a run without samples");
        // csamp: SA[x] - 1 at each run start x, in each character's run order
        u64 nr[256] = {0};
        for (u64 k = 0; k < R; ++k) nr[M.head[k]]++;
        T.Cr[0] = 0;
        for (int c = 0; c < 256; ++c) T.Cr[c + 1] = T.Cr[c] + nr[c];
        std::vector<u64> cs(R), next(T.Cr, T.Cr + 256);
        for (u64 k = 0; k < R; ++k) cs[next[M.head[k]]++] = ssa[k] ? ssa[k] - 1 : n - 1;
        std::vector<move_lf::run>().swap(M.runs);
        T.csamp.build(cs);
        std::vector<u64>().swap(cs);
        // phi^-1: SA at the end of run k -> SA at the start of run k + 1
        std::vector<std::pair<u64, u64>> kv(R - 1);
        for (u64 k = 0; k + 1 < R; ++k) kv[k] = {esa[k], ssa[k + 1]};
        T.last_sa = esa[R - 1];
        std::vector<u64>().swap(ssa); std::vector<u64>().swap(esa);
        std::sort(kv.begin(), kv.end());
        std::vector<u64> kk(R - 1), vv(R - 1);
        for (u64 k = 0; k + 1 < R; ++k) { kk[k] = kv[k].first; vv[k] = kv[k].second; }
        std::vector<std::pair<u64, u64>>().swap(kv);
        T.ikeys.build(kk, n); T.ivals.build(vv);
    }
    fprintf(stderr, "two-level locate: characters: %lu run samples %.3f MB, phi^-1 %.3f MB (%.1f s)\n", B.R, T.csamp.bytes() / 1e6,
            (T.ikeys.bytes() + T.ivals.bytes()) / 1e6, now() - t0);
    // ---- the parse: its suffix array
    const u64 m = L.m, N = m + 1;
    std::vector<uint32_t> Z(N), SA(N);
    for (u64 i = 0; i < m; ++i) Z[i] = (uint32_t)PT.zid(i);
    Z[m] = 0;
    if (libsais_int(reinterpret_cast<int32_t *>(Z.data()), reinterpret_cast<int32_t *>(SA.data()), (int32_t)N, (int32_t)L.sigma + 1, 0) != 0) die("libsais_int failed");
    for (u64 i = 0; i < m; ++i) Z[i] = (uint32_t)PT.zid(i);  // (libsais_int overwrote it)
    Z[m] = 0;
    // the parse BWT's runs, from the EFRank (key c * N + start), checked against the parse
    {
        const EFRank &E = L.pbwt;
        if (E.n != N) die("the parse BWT does not match the parse");
        std::vector<std::pair<u64, u64>> st(E.nruns);
        for (u64 j = 0; j < E.nruns; ++j) {
            u64 key = E.KS.access(j), c = key / N, s = key % N;
            u64 bw = SA[s] == 0 ? 0 : Z[SA[s] - 1];
            if (bw != c) die("internal error: the parse BWT's runs do not match the parse's suffix array");
            st[j] = {s, SA[s]};
        }
        std::sort(st.begin(), st.end());
        std::vector<u64> rows(st.size()), sv(st.size());
        for (size_t q = 0; q < st.size(); ++q) { rows[q] = st[q].first; sv[q] = st[q].second; }
        T.prow.build(rows, N); T.psamp.build(sv);
    }
    // the blocks: the parse position of each block's first phrase (from BS's rows), and the block BWT's runs
    {
        std::vector<u64> bs;
        bs.reserve(L.nblocks);
        for (u64 r = 1; r <= L.nblocks; ++r) bs.push_back(SA[L.BS.access(r)]);
        std::sort(bs.begin(), bs.end());
        if (bs.empty() || bs[0] != 0) die("internal error: the first block does not start at phrase 0");
        T.bstarts.build(bs, m);
        const EFRank &E = L.bbwt;
        const u64 NB = L.nblocks + 1;
        if (E.n != NB) die("the block BWT does not match the blocks");
        std::vector<std::pair<u64, u64>> st(E.nruns);
        for (u64 j = 0; j < E.nruns; ++j) {
            u64 key = E.KS.access(j), s = key % NB;
            st[j] = {s, T.block_of(SA[L.BS.access(s)])};
        }
        std::sort(st.begin(), st.end());
        std::vector<u64> rows(st.size()), sv(st.size());
        for (size_t q = 0; q < st.size(); ++q) { rows[q] = st[q].first; sv[q] = st[q].second; }
        T.brow.build(rows, NB); T.bsamp.build(sv);
    }
    std::vector<uint32_t>().swap(Z); std::vector<uint32_t>().swap(SA);
    {   // the phrase starts
        std::vector<u64> ps(m);
        for (u64 i = 0; i < m; ++i) ps[i] = PT.start(i);
        T.pstarts.build(ps, L.n);
    }
    fprintf(stderr, "two-level locate: parse %lu run samples %.3f MB, blocks %lu run samples %.3f MB, block -> phrase %.3f MB, phrase -> text %.3f MB; "
                    "in all %.3f MB (characters %.3f, parse %.3f, blocks %.3f) (%.1f s)\n",
            L.pbwt.nruns, (T.prow.bytes() + T.psamp.bytes()) / 1e6, L.bbwt.nruns, (T.brow.bytes() + T.bsamp.bytes()) / 1e6, T.bstarts.bytes() / 1e6,
            T.pstarts.bytes() / 1e6, (T.bytes_char() + T.bytes_parse() + T.bytes_blocks()) / 1e6, T.bytes_char() / 1e6, T.bytes_parse() / 1e6,
            T.bytes_blocks() / 1e6, now() - t0);
}

}  // namespace katka
