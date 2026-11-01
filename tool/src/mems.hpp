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
//   min3         ignore MEMs with fewer than min3 distinct 3-mers (they do not raise L either)
// Forward-backward search (fb_segment), as in KATKA's -W: from x, search forward to the largest e such
// that R[x..e) occurs (a MEM, as x is left-maximal by construction); report it if e - x >= L; then the
// next x is the start of the longest suffix of R[0..e+1) that occurs.
#pragma once
#include "digest.hpp"
#include "kmer_filter.hpp"
#include "stpd.hpp"
#include "twolevel.hpp"
#include <algorithm>

namespace katka {

struct mem { u64 a, e; interval I; };    // R[a..e) and the interval of rc(R[a..e)) (or of R[a..e) with FB on chr19)

// longest-suffix searches over one segment of a pattern and over its reverse complement
struct segment_search {
    const stepper *S = nullptr;
    const lookup_table *T = nullptr;      // plain search only
    const two_level *TL = nullptr;        // if set, two-level search
    bool jumps = true;
    const unsigned char *comp = nullptr;  // complement of each symbol (digest symbols), if not DNA
    const stpd_index *SP = nullptr;       // if set, STPD search (stpd.hpp)
    bool intervals = true;                // with STPD: compute the intervals of the reported MEMs
    search_counters cnt;
    stpd_state sst[2];

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
        if (TL) return search2(*S, *TL, dir ? qi : pi, X, (int64_t)lo, (int64_t)end, jumps, I, cnt);
        if (SP) {
            ++cnt.searches;
            stpd_state &st = sst[dir];
            st.dir = dir; st.end = end; st.j = end; st.t = SP->n - 1;     // t = SA[0], the suffix $
            stpd_extend(*SP, X, lo, st, cnt);
            I = {0, 0};
            return st.j;
        }
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
        if (TL) return back(dir, lo, end, I);
        if (SP) {
            stpd_state &st = sst[dir];
            if (st.dir != dir || st.end != end || st.j != j) return back(dir, lo, end, I);
            stpd_extend(*SP, dir ? rc.data() : R, lo, st, cnt);
            return st.j;
        }
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
    // with STPD, the interval of rc(R[a..e)) for a MEM about to be reported, by plain backward search
    void rc_interval(u64 a, u64 e, interval &J) {
        if (!SP) return;
        if (!intervals) { J = {0, 0}; return; }
        if (plain_back(1, n - e, n - a, J) != n - e) die("internal error: STPD and backward search disagree");
    }
};

struct mem_options {
    u64 L = 1;
    bool adaptive = false;          // -X
    u64 top = 0;                    // keep the t longest (0: off)
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
            if (!ok && open) { segs.push_back({r0 + start, p + F->k - 1 - start}); open = false; }
        }
        if (open) segs.push_back({r0 + start, n - start});
    }
}

// The state of one read: L, and the t longest MEMs kept so far.
struct read_state {
    u64 L;
    std::vector<mem> kept;          // top t: a min-heap by length (then by position, for determinism)
    static bool longer(const mem &x, const mem &y) { return x.e - x.a != y.e - y.a ? x.e - x.a > y.e - y.a : x.a < y.a; }
};

// Records a MEM; returns false if it is dropped by the low-complexity filter.
template <class F>
inline void found(const char *R, u64 r0, u64 a, u64 e, const interval &J, const mem_options &opt, read_state &st, F report) {
    if (opt.min3 && distinct_3mers(R + a, e - a) < opt.min3) return;
    mem M{r0 + a, r0 + e, J};
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
        found(R, r0, a, e, J, opt, st, report);
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
        if (e - s0 >= st.L) { Q.rc_interval(s0, e, J); found(Q.R, r0, s0, e, J, opt, st, report); }
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
    segments(read, opt.pseudo ? opt.filter : nullptr, segs);
    const bool longest_first = (opt.adaptive || opt.top) && opt.pseudo;
    if (longest_first)
        std::stable_sort(segs.begin(), segs.end(), [](const std::pair<u64, u64> &p, const std::pair<u64, u64> &q) { return p.second > q.second; });
    read_state st;
    st.L = opt.L;
    for (auto [r0, n] : segs) {
        if (n < st.L) { if (longest_first) break; continue; }
        Q.prepare(read.data() + r0, n);
        if (opt.fb) fb_segment(Q, r0, opt, st, report);
        else bml_segment(Q, r0, opt, st, report);
    }
    if (opt.top) {
        std::sort(st.kept.begin(), st.kept.end(), [](const mem &x, const mem &y) { return x.a < y.a; });
        for (auto &M : st.kept) report(M);
    }
}

}  // namespace katka
