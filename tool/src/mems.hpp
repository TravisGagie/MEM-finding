// Finding MEMs (Section 2 of the paper, "MEMs", "Boyer-Moore-Li", "Longest MEMs", "KeBaB").
//
// Every algorithm here uses only one primitive: the longest suffix X[j..end) of a substring X[lo..end) of
// the pattern, or of its reverse complement, that occurs in the text, with its BWT interval.  The
// segment_search below provides it by plain backward search (on the RLBWT or the RLCSA, optionally
// starting with a lookup table) or by two-level search (twolevel.hpp); both return the same j and
// interval, so BML and forward-backward search give the same output with either.
//
// BML (bml_segment), as in KATKA:
//   x = 0
//   while x + L <= n:
//     [inline KeBaB] if a k-mer of R[x..x+L) fails the filter, jump past the rightmost one that fails
//     search backward for R[x..x+L); if it stops at j > x, set x = j and repeat
//     otherwise extend left to the smallest a such that R[a..x+L) occurs, then forward from a to the
//     largest e such that R[a..e) occurs: R[a..e) is a MEM, with the interval of rc(R[a..e))
//     set x = e - L + 1
// Options:
//   adaptive L   after a MEM longer than L, set L to its length (so that the longest MEMs, with ties, are
//                among those reported; start with L = 1 to get the longest overall)
//   top t        keep the t longest MEMs found so far; once there are t, L = 1 + the shortest of them
//   KeBaB        (a) pseudo-MEMs: search only the maximal substrings whose k-mers all pass the filter,
//                longest first with adaptive L or top t, stopping at the first one shorter than L;
//                (b) inline: the window check above, testing the window's k-mers lazily from its right end
//                and stopping at the first that fails (O(k + 1/p) hash steps per failed window, for a
//                fraction p of failing k-mers), before any backward step
//   longest t    the t longest MEMs with ties: every MEM at least as long as the t-th longest (keeping, once
//                t are found, L = the length of the t-th longest so far, so ties are found too)
//   KeBaB + longest t, three rules for which pseudo-MEMs (sorted longest first) to search with BML:
//                (in pseudo-MEMs, L is at least k: a shorter match there need not be a MEM)
//                streaming: until the next one is shorter than L (exact if L >= k: every MEM of length >= k
//                  lies in a pseudo-MEM, and is no longer than it);
//                top (KeBaB's heuristic, --kebab-top): only the t longest pseudo-MEMs (no guarantee);
//                batch (--kebab-batch): the 2t longest, then batches of doubling size (4t, 8t, ...); after
//                  each, with lambda the length of the next pseudo-MEM (k - 1 if none is left), stop if at
//                  least t MEMs found are longer than lambda (exact: the MEMs not yet found are at most
//                  lambda long; strictly longer, so ties are complete); if the pseudo-MEMs run out with
//                  fewer than t MEMs of length >= k, BML on the whole read (only if L < k: with L >= k,
//                  every MEM it could find lies in a pseudo-MEM already searched)
//                check (--kebab-check T t): the T longest pseudo-MEMs (T >= t), then accept if there are at
//                  least T + 1 pseudo-MEMs and at least t MEMs found are at least as long as the (T+1)st
//                  (--strict: longer), else BML on the whole read (exact; with --strict, ties complete too),
//                  from L, or (--check-seed) from the t-th longest MEM found if at least t were found
//   KeBaB top with a fixed L (--kebab-top t --fixed-L, KeBaB's pipeline): every MEM of length >= max(L, k)
//                in the t longest pseudo-MEMs, nothing else
//   min3         ignore MEMs with fewer than min3 distinct 3-mers (they do not raise L either)
// Forward-backward search (fb_segment), as in KATKA's -W: from x, search forward to the largest e such
// that R[x..e) occurs (a MEM, as x is left-maximal by construction); report it if e - x >= L; then the
// next x is the start of the longest suffix of R[0..e+1) that occurs.
#pragma once
#include "digest.hpp"
#include "kmer_filter.hpp"
#include "stpd.hpp"
#include "pstpd.hpp"
#include "twolevel.hpp"
#include <algorithm>

namespace katka {

struct mem { u64 a, e; interval I; u64 occ = 0, haps = 0; };   // haps: distinct sequences among the located occurrences
   // occ: occurrences located (STPD with locate), else 0    // R[a..e) and the interval of rc(R[a..e)) (or of R[a..e) with FB on chr19)

// longest-suffix searches over one segment of a pattern and over its reverse complement
struct segment_search {
    const stepper *S = nullptr;
    const lookup_table *T = nullptr;      // plain search only
    const two_level *TL = nullptr;        // if set, two-level search
    const twolevel_locate *TLC = nullptr; // with TL and locate (no STPD): the toehold samples and phi^-1
    bool jumps = true;
    const unsigned char *comp = nullptr;  // complement of each symbol (digest symbols), if not DNA
    const stpd_index *SP = nullptr;       // if set, STPD search (stpd.hpp)
    bool locate = false;                  // with STPD: locate every occurrence of each reported MEM by phi^-1
    u64 last_occ = 0, locate_checks = 0;  //   (the number of occurrences of the last one, and text checks so far)
    std::vector<u64> occbuf;              //   (and their positions: of rc(MEM) in S, for BML and FB)
    bool intervals = true;                // with STPD: compute the intervals of the reported MEMs (by two-level
                                          // search if TL is set too, else by plain backward search)
    search_counters cnt;
    stpd_state sst[2];
    // STPD over the parse (PS, with SP's character-level sample and its parse oracle, and TL's B and
    // dictionary): the state of the last search in each direction
    const parse_stpd *PS = nullptr;
    const std::vector<u64> *gstarts = nullptr;   // with locate: the starts of the sequences in S (haplotypes)
    u64 last_haps = 0;
    struct pstate {
        int dir = -1; u64 end = 0, j = 0; bool failed = false;
        int mode = 0;                     // 0: characters (sst[dir]); 1: parse (u, i)
        u64 u = 0; int64_t i = 0;         // parse mode: Z[u..] is the smallest parse suffix in X[S[i]..end)'s range
        int64_t ridx = -1;                // the last syncmer start S[ridx] with S[ridx] + k <= end
        u64 core_u = 0; int64_t core_i = -1;   // at the last exit from parse mode
    } pst[2];

    const char *R = nullptr;
    u64 n = 0;
    std::string rc;
    pattern_parse pi, qi;

    void prepare(const char *r, u64 len) {
        R = r; n = len;
        rc.resize(n);
        for (u64 t = 0; t < n; ++t) {
            unsigned char c = (unsigned char)R[n - 1 - t];
            rc[t] = (char)(comp ? comp[c] : complement(c));
        }
        if (TL) parse_pattern(*TL, R, n, pi, qi);
    }
    // the longest suffix of X[lo..end), X = R (dir 0) or rc(R) (dir 1); returns its start j and interval I
    u64 back(int dir, u64 lo, u64 end, interval &I) {
        const char *X = dir ? rc.data() : R;
        if (PS) return ps_back(dir, lo, end);
        if (SP) {
            ++cnt.searches;
            stpd_state &st = sst[dir];
            st.dir = dir; st.end = end; st.j = end; st.t = SP->n - 1;     // t = SA[0], the suffix $
            stpd_extend(*SP, X, lo, st, cnt);
            I = {0, 0};
            return st.j;
        }
        if (TL) return search2(*S, *TL, dir ? qi : pi, X, (int64_t)lo, (int64_t)end, jumps, I, cnt);
        return plain_back(dir, lo, end, I);
    }
    u64 plain_back(int dir, u64 lo, u64 end, interval &I) {
        const char *X = dir ? rc.data() : R;
        ++cnt.searches;
        I = S->full();
        u64 j = end;
        if (T && T->t && end - lo >= (u64)T->t) {
            interval J;
            if (T->get((const unsigned char *)X + end - T->t, J)) { I = J; j = end - T->t; }
        }
        while (j > lo) {
            ++cnt.char_steps;
            if (!S->step(I, (unsigned char)X[j - 1])) break;
            --j;
        }
        return j;
    }
    // X[j..end) occurs with interval I: extend it to the left down to lo
    u64 extend_left(int dir, u64 lo, u64 j, u64 end, interval &I) {
        if (PS) {
            pstate &P = pst[dir];
            if (P.dir != dir || P.end != end || P.j != j || P.failed) return ps_back(dir, lo, end);
            ps_extend(dir, lo);
            return P.j;
        }
        if (SP) {
            stpd_state &st = sst[dir];
            if (st.dir != dir || st.end != end || st.j != j) return back(dir, lo, end, I);
            stpd_extend(*SP, dir ? rc.data() : R, lo, st, cnt);
            return st.j;
        }
        if (TL) return back(dir, lo, end, I);
        const char *X = dir ? rc.data() : R;
        while (j > lo) {
            ++cnt.char_steps;
            if (!S->step(I, (unsigned char)X[j - 1])) break;
            --j;
        }
        return j;
    }
    // the largest e <= n such that R[a..e) occurs, with the interval of rc(R[a..e))
    // (with STPD, J is computed only by rc_interval, for the MEMs reported)
    u64 forward(u64 a, interval &J) { return n - back(1, 0, n - a, J); }
    // ---- STPD over the parse (NOTES.md, "STPD over the parse")
    // the longest suffix of X[lo..end): right end by backward steps to the last syncmer start, entry into
    // the parse, then ps_extend
    u64 ps_back(int dir, u64 lo, u64 end) {
        const char *X = dir ? rc.data() : R;
        const pattern_parse &pp = dir ? qi : pi;
        const auto &S = pp.S;
        const int64_t K = TL->K;
        ++cnt.searches;
        pstate &P = pst[dir];
        P = pstate(); P.dir = dir; P.end = end; P.j = end;
        P.ridx = (int64_t)(std::upper_bound(S.begin(), S.end(), (int64_t)end - K) - S.begin()) - 1;
        if (P.ridx >= 0 && S[P.ridx] >= (int64_t)lo) {
            interval J = this->S->full();
            u64 j = end;
            while (j > (u64)S[P.ridx]) {
                ++cnt.ps_right;
                if (!this->S->step(J, (unsigned char)X[j - 1])) break;
                --j;
            }
            if (j > (u64)S[P.ridx]) {                      // failed at the right end
                P.j = j; P.failed = true; P.mode = 0;
                if (locate && dir == 1) {                    // (the smallest occurrence, for locating: forward searches only)
                    sst[dir] = {end, end, SP->n - 1, dir};
                    stpd_extend(*SP, X, j, sst[dir], cnt);
                }
                return j;
            }
            const u64 rho = TL->B.rank1(J.sp) + 1;          // the first parse row in X[S[ridx]..end)'s range
            const u64 y = PS->id_of_row(rho);
            if (PS->c_of(y) != rho) die("internal error: parse STPD entry (the range does not start an id's block)");
            P.u = PS->entry.get(y); P.i = P.ridx; P.mode = 1; P.j = (u64)S[P.ridx];
        } else {
            P.mode = 0;
            sst[dir] = {end, end, SP->n - 1, dir};          // t = SA[0], the suffix $
        }
        ps_extend(dir, lo);
        return P.j;
    }
    // continues the search of pst[dir] down to lo
    void ps_extend(int dir, u64 lo) {
        const char *X = dir ? rc.data() : R;
        const uint8_t *x = (const uint8_t *)X;
        const pattern_parse &pp = dir ? qi : pi;
        const auto &S = pp.S;
        const parse_text &PT = *SP->PT;
        pstate &P = pst[dir];
        for (;;) {
            if (P.mode == 1) {
                bool fail = false;
                const u64 xr = (u64)S[P.ridx], lr = P.end - xr;
                while (P.i > 0 && S[P.i - 1] >= (int64_t)lo && pp.pid[P.i - 1]) {
                    const uint32_t a = pp.pid[P.i - 1];
                    ++cnt.ps_mid;
                    if (P.u > 0 && PT.zid(P.u - 1) == a) { --P.u; --P.i; continue; }
                    ++cnt.ps_jumps;
                    u64 v = PS->find(PT, &pp.pid[P.i - 1], (u64)(P.ridx - P.i + 1), x + xr, lr, cnt.ps_compares);
                    if (v == UINT64_MAX) { fail = true; break; }
                    P.u = v; --P.i;
                }
                P.core_u = P.u; P.core_i = P.i;
                sst[dir] = {P.end, (u64)S[P.i], PT.start(P.u), dir};
                P.j = (u64)S[P.i]; P.mode = 0;
                if (fail) {                                  // the match ends inside phrase i - 1
                    u64 c0 = cnt.char_steps;
                    stpd_extend(*SP, X, lo, sst[dir], cnt);
                    cnt.ps_left += cnt.char_steps - c0;
                    P.j = sst[dir].j; P.failed = P.j > lo;
                    return;
                }
            }
            // characters: switch to the parse at a syncmer start whose phrase to the left can be matched
            const int64_t ic = std::lower_bound(S.begin(), S.end(), (int64_t)P.j) - S.begin();
            if (ic < (int64_t)S.size() && S[ic] == (int64_t)P.j && ic <= P.ridx && ic > 0 && S[ic - 1] >= (int64_t)lo && pp.pid[ic - 1]) {
                const u64 t = sst[dir].t, u = PT.pred(t);
                if (PT.start(u) != t) die("internal error: parse STPD (an occurrence of a syncmer not at a boundary)");
                P.u = u; P.i = ic; P.mode = 1;
                continue;
            }
            if (P.j <= lo) return;
            const int64_t ii = std::min<int64_t>(ic - 1, P.ridx);
            const u64 stop = ii >= 0 && S[ii] >= (int64_t)lo ? (u64)S[ii] : lo;
            u64 c0 = cnt.char_steps;
            stpd_extend(*SP, X, stop, sst[dir], cnt);
            cnt.ps_left += cnt.char_steps - c0;
            P.j = sst[dir].j;
            if (P.j > stop) { P.failed = true; return; }
        }
    }
    // every occurrence of rc(R[a..e)) (= X[j..end), X = rc, after the forward search), into occbuf: by
    // phi^-1 over the parse from the smallest occurrence of its core X[S[i0]..end) (i0: its first syncmer),
    // checking phrase ids and the right end's characters, then the left end's characters; without a core,
    // the character-level locate.  Returns the number of occurrences.
    u64 ps_locate(u64 a, u64 e) {
        const pattern_parse &pp = qi;
        const auto &S = pp.S;
        const parse_text &PT = *SP->PT;
        const uint8_t *x = (const uint8_t *)rc.data();
        pstate &P = pst[1];
        const u64 j = n - e, end = n - a;
        occbuf.clear();
        const int64_t i0 = std::lower_bound(S.begin(), S.end(), (int64_t)j) - S.begin();
        if (P.end != end || P.j != j) die("internal error: parse STPD locate without its search");
        if (P.core_i < 0 || P.core_i != i0 || i0 > P.ridx) {
            if (sst[1].end != end || sst[1].j != j) die("internal error: parse STPD locate (character state)");
            return SP->locate(x + j, end - j, sst[1].t, occbuf, locate_checks);
        }
        const u64 h = (u64)(P.ridx - i0), xr = (u64)S[P.ridx], lr = end - xr, ll = (u64)S[i0] - j;
        const uint32_t *ids = &pp.pid[i0];
        auto occurs = [&](u64 v) { u64 l; ++locate_checks; return parse_stpd::cmp(PT, v, ids, h, x + xr, lr, 0, l) == 0; };
        static thread_local std::vector<u64> core;
        core.clear();
        core.push_back(P.core_u);
        for (u64 batch = 1;; batch *= 2) {                  // as stpd_index::locate, over the parse
            const u64 b0 = core.size();
            bool last = false;
            for (u64 q = 0; q < batch; ++q) {
                if (core.back() == PS->last_sa) { last = true; break; }
                core.push_back(PS->phi_inv(core.back()));
            }
            if (core.size() == b0) break;
            if (!occurs(core.back())) {
                u64 lo2 = b0, hi2 = core.size() - 1;
                while (lo2 < hi2) { u64 mid = (lo2 + hi2) / 2; if (occurs(core[mid])) lo2 = mid + 1; else hi2 = mid; }
                core.resize(lo2);
                break;
            }
            if (last) break;
        }
        for (u64 v : core) {                                 // the left end
            u64 p = PT.start(v);
            ++locate_checks;
            if (p >= ll && PT.lce(p - ll, x + j, ll) == ll) occbuf.push_back(p - ll);
        }
        return occbuf.size();
    }
    // with STPD, the interval of rc(R[a..e)) for a MEM about to be reported, by plain backward search
    void rc_interval(u64 a, u64 e, interval &J) {
        if (!SP && !(locate && TLC)) return;
        interval I2;
        if (locate) {                                   // the forward search left the first occurrence of rc(R[a..e)) in sst[1]
            occbuf.clear();
            if (!SP) {                                  // two-level search with a toehold: SA[sp], then phi^-1
                u64 t = 0;
                u64 j2 = search2(*S, *TL, qi, rc.data(), (int64_t)(n - e), (int64_t)(n - a), jumps, I2, cnt, TLC, S->bwt, &t);
                if (j2 != n - e || I2.sp >= I2.ep) die("internal error: two-level locate (the MEM does not occur)");
                for (u64 x = I2.sp; x < I2.ep; ++x) { occbuf.push_back(t); if (x + 1 < I2.ep) t = TLC->phi_inv(t); }
                last_occ = occbuf.size();
            } else
            last_occ = PS ? ps_locate(a, e) : SP->locate((const uint8_t *)rc.data() + (n - e), e - a, sst[1].t, occbuf, locate_checks);
            if (gstarts) {                          // distinct sequences among the occurrences
                static thread_local std::vector<u64> hs;
                hs.clear();
                for (u64 p : occbuf) hs.push_back(std::upper_bound(gstarts->begin(), gstarts->end(), p) - gstarts->begin());
                std::sort(hs.begin(), hs.end());
                last_haps = std::unique(hs.begin(), hs.end()) - hs.begin();
            }
        }
        if (!SP) { J = I2; return; }                    // (two-level locating gave the interval too)
        if (!intervals) { J = {0, 0}; return; }
        // by two-level search (with blocks) if there is a two-level index, else by plain backward search
        u64 j = TL && TL->pbwt.n ? search2(*S, *TL, qi, rc.data(), (int64_t)(n - e), (int64_t)(n - a), jumps, J, cnt) : plain_back(1, n - e, n - a, J);
        if (j != n - e) die("internal error: STPD and backward search disagree");
    }
};

struct mem_options {
    u64 L = 1;
    bool adaptive = false;          // -X
    u64 top = 0;                    // keep the t longest (0: off)
    u64 longest = 0;                // the t longest with ties (0: off)
    enum { STREAM, TOP, BATCH, CHECK } kebab_rule = STREAM;   // with pseudo and longest
    u64 check_T = 0;                // CHECK: search the T longest pseudo-MEMs, then test (T >= t)
    bool check_strict = false;      //   accept only with t MEMs longer than the (T+1)st pseudo-MEM
    bool check_seed = false;        //   the fallback starts with L = the t-th longest MEM found, if t were found
    u64 fixed_top = 0;              // with pseudo, no longest: every MEM >= L in the t longest pseudo-MEMs only
    int min3 = 0;                   // -D
    const kmer_filter *filter = nullptr;
    bool pseudo = false;            // KeBaB pseudo-MEMs (-K)
    bool inline_check = false;      // KeBaB inline window check
    bool fb = false;                // forward-backward instead of BML
    const digest_map *digest = nullptr;   // BML over minimizer digests (the index is of the digested text)
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

// The segments to search, as (start, length): the maximal A/C/G/T runs of the read, or with pseudo-MEMs
// the maximal substrings of those runs all of whose k-mers pass the filter.
inline void segments(const std::string &read, const kmer_filter *F, std::vector<std::pair<u64, u64>> &segs, u64 *tests = nullptr) {
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
        if (tests) *tests += kmers.size();
        u64 start = 0; bool open = false;
        for (u64 p = 0; p < kmers.size(); ++p) {
            if (p + 8 < kmers.size()) F->prefetch(kmers[p + 8]);
            bool ok = F->has(kmers[p]);
            if (ok && !open) { start = p; open = true; }
            if (!ok && open) { segs.push_back({r0 + start, p + F->k - 1 - start}); open = false; }
        }
        if (open) segs.push_back({r0 + start, n - start});
    }
}

// The state of one read: L, and the t longest MEMs kept so far.
struct read_state {
    u64 L;
    std::vector<mem> kept;          // top t: a min-heap by length (then by position, for determinism)
    std::vector<mem> best;          // longest t: every MEM found at least as long as the t-th longest so far
    static bool longer(const mem &x, const mem &y) { return x.e - x.a != y.e - y.a ? x.e - x.a > y.e - y.a : x.a < y.a; }
};

// Records a MEM; returns false if it is dropped by the low-complexity filter.
template <class F>
inline void found(const char *R, u64 r0, u64 a, u64 e, const interval &J, const mem_options &opt, read_state &st, F report, u64 occ = 0, u64 haps = 0) {
    if (opt.min3 && distinct_3mers(R + a, e - a) < opt.min3) return;
    mem M{r0 + a, r0 + e, J, occ, haps};
    if (opt.longest) {
        st.best.push_back(M);
        if (st.best.size() >= opt.longest) {
            static thread_local std::vector<u64> lens;
            lens.clear();
            for (auto &x : st.best) lens.push_back(x.e - x.a);
            std::nth_element(lens.begin(), lens.begin() + (opt.longest - 1), lens.end(), std::greater<u64>());
            const u64 lt = lens[opt.longest - 1];       // the t-th longest so far
            st.best.erase(std::remove_if(st.best.begin(), st.best.end(), [&](const mem &x) { return x.e - x.a < lt; }), st.best.end());
            st.L = std::max(st.L, lt);
        }
        return;
    }
    if (opt.top) {
        if (st.kept.size() < opt.top) { st.kept.push_back(M); std::push_heap(st.kept.begin(), st.kept.end(), read_state::longer); }
        else if (e - a > st.kept.front().e - st.kept.front().a) {
            std::pop_heap(st.kept.begin(), st.kept.end(), read_state::longer);
            st.kept.back() = M;
            std::push_heap(st.kept.begin(), st.kept.end(), read_state::longer);
        }
        if (st.kept.size() == opt.top) st.L = std::max(st.L, st.kept.front().e - st.kept.front().a + 1);
        return;
    }
    report(M);
    if (opt.adaptive && e - a > st.L) st.L = e - a;
}

// BML on the prepared segment Q.R[0..n), at offset r0 in the read
template <class F>
void bml_segment(segment_search &Q, u64 r0, const mem_options &opt, read_state &st, F report) {
    const u64 n = Q.n;
    const char *R = Q.R;
    const kmer_filter *KF = opt.inline_check ? opt.filter : nullptr;
    u64 x = 0;
    while (x + st.L <= n) {
        const u64 L = st.L;
        if (KF && L >= KF->k) {
            // inline KeBaB, lazily: the k-mers of R[x..x+L) from the right end leftward, with the rolling
            // codes restarted at the window's last k-mer; jump past the first one that fails
            const u64 k = KF->k, mask = k == 32 ? ~0ULL : (1ULL << (2 * k)) - 1;
            u64 p = x + L - k, fw = 0, rcv = 0;
            for (u64 i = 0; i < k; ++i) {
                u64 c = (u64)base_code(R[p + i]);
                fw = ((fw << 2) | c) & mask;
                rcv = (rcv >> 2) | ((3 - c) << (2 * (k - 1)));
            }
            Q.cnt.kmer_restarts += 1;
            int64_t fail = -1;
            for (;;) {
                ++Q.cnt.kmer_tests;
                if (!KF->has(fw < rcv ? fw : rcv)) { fail = (int64_t)p; break; }
                if (p == x) break;
                --p;
                u64 c = (u64)base_code(R[p]);
                fw = (fw >> 2) | (c << (2 * (k - 1)));
                rcv = ((rcv << 2) | (3 - c)) & mask;
            }
            if (fail >= 0) { x = (u64)fail + 1; continue; }
        }
        interval I;
        u64 j = Q.back(0, x, x + L, I);
        if (j > x) { x = j; continue; }                     // R[j-1..x+L) does not occur
        u64 a = Q.extend_left(0, 0, x, x + L, I);
        interval J;
        u64 e = Q.forward(a, J);
        Q.rc_interval(a, e, J);
        found(R, r0, a, e, J, opt, st, report, Q.locate ? Q.last_occ : 0, Q.last_haps);
        // the first window that ends after the MEM (L may have been raised; with top t, past e + 1)
        x = e + 1 >= st.L ? e + 1 - st.L : 0;
    }
}

// forward-backward search on the prepared segment
template <class F>
void fb_segment(segment_search &Q, u64 r0, const mem_options &opt, read_state &st, F report) {
    const u64 m = Q.n;
    u64 s0 = 0;
    while (s0 < m) {
        interval J, I;
        u64 e = Q.forward(s0, J);                           // R[s0..e) is the MEM starting at s0
        if (e <= s0) { ++s0; continue; }
        if (e - s0 >= st.L) { Q.rc_interval(s0, e, J); found(Q.R, r0, s0, e, J, opt, st, report, Q.locate ? Q.last_occ : 0, Q.last_haps); }
        if (e >= m) break;
        u64 s2 = Q.back(0, 0, e + 1, I);                   // leftmost start of a match ending at R[e]
        s0 = s2 > e ? e + 1 : s2;
    }
}

// BML over minimizer digests (App. B "Minimizer digests", App. A "BML with digests"), on one A/C/G/T run
// R[0..n) of a read, n >= max(L, w), with L in bases.  Window y covers bases y..y+w-1 and selects digest
// symbols lo[y]..hi[y]; a DNA match R[x..x+L) contains exactly the windows x..x+L-w, so its digest is
// dg[lo[x]..hi[x+L-w]], which must occur if R[x..x+L) does.
//   search backward for dg[lo[x]..hi[x+L-w]]; if it fails at symbol j-1, jump to the first window y with
//     lo[y] > j-1 (every x' with lo[x'] <= j-1 fails too)
//   otherwise extend to a maximal digested MEM dg[a..b] and take its span in bases (digest.hpp); continue
//     from the first window start whose last window has hi > b
// Finally, MEMs whose span lies inside another's are dropped (of equal spans, the first is kept).  This
// finds every digested MEM whose span is at least L and is not contained in another's span.
template <class F>
void bml_digest_run(segment_search &Q, const digest_map &M, const char *R, u64 n, u64 r0, const mem_options &opt, F report) {
    static thread_local read_digest D;
    D.build(M, R, n);
    const u64 L = opt.L, w = M.w, nwin = D.lo.size(), nd = D.dg.size();
    Q.prepare(D.dg.data(), nd);
    struct found_mem { u64 l, r; interval I; };
    static thread_local std::vector<found_mem> fr;
    fr.clear();
    u64 x = 0;
    while (x + L <= n) {
        u64 a = D.lo[x], b = D.hi[x + L - w];
        interval I;
        u64 j = Q.back(0, a, b + 1, I);
        if (j > a) {                                        // dg[j-1..b] does not occur
            u64 nx = std::upper_bound(D.lo.begin(), D.lo.end(), j - 1) - D.lo.begin();
            if (nx >= nwin) break;
            x = nx; continue;
        }
        u64 aa = Q.extend_left(0, 0, a, b + 1, I);
        interval J;
        u64 e = Q.forward(aa, J);
        Q.rc_interval(aa, e, J);
        u64 left, right, yR;
        D.span(aa, e - 1, n, (int)w, left, right, yR);
        fr.push_back({left, right, J});
        if (yR >= nwin) break;
        u64 nx = yR >= L - w ? yR - (L - w) : 0;
        x = nx > x ? nx : x + 1;
    }
    for (u64 p = 0; p < fr.size(); ++p) {
        bool in = false;
        for (u64 q = 0; q < fr.size() && !in; ++q)
            if (q != p && fr[q].l <= fr[p].l && fr[p].r <= fr[q].r && (fr[q].l < fr[p].l || fr[p].r < fr[q].r || q < p)) in = true;
        if (in) continue;
        if (opt.min3 && distinct_3mers(R + fr[p].l, fr[p].r + 1 - fr[p].l) < opt.min3) continue;
        report(mem{r0 + fr[p].l, r0 + fr[p].r + 1, fr[p].I});
    }
}

// All the MEMs of a read, by BML or forward-backward search; calls report(mem) in the order found (with
// top t, the t longest at the end, sorted by position).
template <class F>
void find_mems(segment_search &Q, const std::string &read, const mem_options &opt, F report) {
    static thread_local std::vector<std::pair<u64, u64>> segs;
    if (opt.digest) {
        segments(read, nullptr, segs);
        for (auto [r0, n] : segs)
            if (n >= opt.L && n >= (u64)opt.digest->w) bml_digest_run(Q, *opt.digest, read.data() + r0, n, r0, opt, report);
        return;
    }
    segments(read, opt.pseudo ? opt.filter : nullptr, segs, &Q.cnt.kmer_tests);
    const bool longest_first = (opt.adaptive || opt.top || opt.longest || opt.fixed_top) && opt.pseudo;
    if (longest_first)
        std::stable_sort(segs.begin(), segs.end(), [](const std::pair<u64, u64> &p, const std::pair<u64, u64> &q) { return p.second > q.second; });
    read_state st;
    st.L = opt.L;
    auto search = [&](u64 r0, u64 n) {
        Q.prepare(read.data() + r0, n);
        if (opt.fb) fb_segment(Q, r0, opt, st, report);
        else bml_segment(Q, r0, opt, st, report);
    };
    if (opt.longest && opt.pseudo) {
        const u64 t = opt.longest, ns = segs.size(), k = opt.filter->k;
        Q.cnt.kb_pseudo += ns;
        // inside pseudo-MEMs, only matches of length >= k are sure to be MEMs (one that ends at a pseudo-MEM's
        // boundary could extend past it if shorter than k): shorter MEMs come only from the whole-read fallback
        st.L = std::max<u64>(st.L, k);
        auto search_pseudo = [&](u64 i) {                // (one shorter than L holds no MEM of length >= L)
            if (segs[i].second < st.L) return;
            ++Q.cnt.kb_searched;
            search(segs[i].first, segs[i].second);
        };
        if (opt.kebab_rule == mem_options::CHECK) {
            const u64 T = opt.check_T;
            for (u64 i = 0; i < ns && i < T; ++i) search_pseudo(i);
            Q.cnt.kb_batches += 1;
            int why = 0;                                 // 0: accepted; 1: at most T pseudo-MEMs; 2: too few long MEMs
            if (ns < T + 1) why = 1;
            else {
                const u64 lambda = segs[T].second;
                u64 c = 0;
                for (auto &M : st.best) c += opt.check_strict ? M.e - M.a > lambda : M.e - M.a >= lambda;
                if (c < t) why = 2;
            }
            if (why) {
                ++Q.cnt.kb_fallbacks;
                (why == 1 ? Q.cnt.kb_fb_few : Q.cnt.kb_fb_short) += 1;
                // --check-seed: every MEM found is a MEM of the read (of length >= k), so the t-th longest found
                // is at most the t-th longest overall, and the whole-read search can start there
                u64 seed = opt.L;
                if (opt.check_seed && st.best.size() >= t) {
                    std::vector<u64> ls;
                    for (auto &M : st.best) ls.push_back(M.e - M.a);
                    std::nth_element(ls.begin(), ls.begin() + (t - 1), ls.end(), std::greater<u64>());
                    seed = std::max(seed, ls[t - 1]);
                    ++Q.cnt.kb_seeded;
                }
                st = read_state(); st.L = seed;
                segments(read, nullptr, segs);
                for (auto [r0, n] : segs) if (n >= st.L) search(r0, n);
            }
        } else if (opt.kebab_rule == mem_options::STREAM) {
            for (u64 i = 0; i < ns && segs[i].second >= st.L; ++i) search_pseudo(i);
        } else if (opt.kebab_rule == mem_options::TOP) {
            for (u64 i = 0; i < ns && i < t; ++i) search_pseudo(i);
        } else {
            u64 i = 0, size = 2 * t, batches = 0;
            for (;;) {
                const u64 end = std::min(ns, i + size);
                for (; i < end; ++i) search_pseudo(i);
                ++batches;
                const u64 lambda = i < ns ? segs[i].second : k - 1;
                u64 longer = 0;
                for (auto &M : st.best) longer += M.e - M.a > lambda;
                if (longer >= t) break;
                if (i >= ns) {                           // the pseudo-MEMs ran out (and fewer than t are >= k)
                    if (opt.L < k) {
                        ++Q.cnt.kb_fallbacks;
                        st = read_state(); st.L = opt.L;
                        segments(read, nullptr, segs);
                        for (auto [r0, n] : segs) if (n >= st.L) search(r0, n);
                    }
                    break;
                }
                size *= 2;
            }
            Q.cnt.kb_batches += batches;
            Q.cnt.kb_multibatch += batches > 1;
        }
    } else if (opt.fixed_top && opt.pseudo) {             // KeBaB's pipeline: the t longest pseudo-MEMs, fixed L
        Q.cnt.kb_pseudo += segs.size();
        st.L = std::max<u64>(st.L, opt.filter->k);       // (in pseudo-MEMs, only matches >= k are sure to be MEMs)
        for (u64 i = 0; i < segs.size() && i < opt.fixed_top; ++i)
            if (segs[i].second >= st.L) { ++Q.cnt.kb_searched; search(segs[i].first, segs[i].second); }
    } else {
        for (auto [r0, n] : segs) {
            if (n < st.L) { if (longest_first) break; continue; }
            search(r0, n);
        }
    }
    if (opt.top || opt.longest) {
        std::vector<mem> &v = opt.top ? st.kept : st.best;
        std::sort(v.begin(), v.end(), [](const mem &x, const mem &y) { return x.a < y.a; });
        for (auto &M : v) report(M);
    }
}

}  // namespace katka
