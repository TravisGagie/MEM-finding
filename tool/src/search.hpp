// Backward search (Section 2, "BWTs"): a BWT interval [sp, ep) and a step that prepends one character,
// on the RLBWT or on the RLCSA (both give the same intervals).  The text holds both strands, so the
// forward extension of a match R[a..e) is the backward extension of rc(R[a..e)) by the complement of R[e].
//
// The lookup table (-F t; KATKA App. A, "Lookup tables") stores the interval of every string of t bases,
// so a search can start with one lookup instead of t steps.  When the string does not occur, the caller
// falls back to stepping base by base, so BML still learns exactly where the search fails (or, with --absent,
// reads that from the table), and the answers are the same with and without a table.
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
    // step on the RLCSA, within c's runs rl..rh (see lookup_table::narrow)
    bool step_in(interval &I, unsigned char c, u64 rl, u64 rh) const {
        u64 a = I.sp, b = I.ep;
        csa->extend_in(a, b, c, rl, rh);
        if (a >= b) return false;
        I = {a, b};
        return true;
    }
};

// The code of the last t symbols prepended, kept as a search moves left (lookup_table::narrow)
struct table_key {
    u64 key = 0;
    int have = 0;                   // how many symbols of X[j..) it holds, up to t
};

// The table is over an alphabet: A, C, G, T by default, or the symbols of a digest (set_alphabet).
struct lookup_table {
    int t = 0;                      // 0: no table
    std::vector<unsigned char> alpha{'A', 'C', 'G', 'T'};
    int16_t code[256];
    std::vector<interval> iv;       // indexed by the string's code in base |alpha|, first symbol most significant
    std::vector<uint8_t> occurs;
    // --narrow, on the RLCSA: when a step extends a match X (|X| >= t - 1) by c, the first of c's runs at or after either end of
    // X's interval lies between the two for X[0..t-1)'s interval, which contains it.  Both are stored for each t-mer c X[0..t-1),
    // so the binary searches start there.  A t-mer that doesn't occur means the step fails, with no search.
    bool narrow = false;
    std::vector<std::pair<uint32_t, uint32_t>> runs;
    // --absent: for a t-mer that doesn't occur, iv[key].sp is the length of its longest suffix that does (it costs no space).  So
    // when BML with t < L finds its window's last t-mer absent, the failure point comes from that one lookup, and the window's
    // search needs no backward steps.
    bool absent = false;
    // --compact: the intervals of the occurring t-mers are disjoint and appear down the BWT in the t-mers' order, so instead of the
    // arrays above: OC, a bitvector over the t-mers, with rank; ES and EE, Elias-Fano lists of the r-th occurring t-mer's
    // interval start and end; AB, the failure length of each absent t-mer, in 4 bits (t <= 13); and for --narrow, ER1 and ER2,
    // Elias-Fano lists of its two run bounds, offset by the first run of its first symbol so that they increase too.  An
    // interval takes one rank and two selects, instead of one array read.
    bool compact = false;
    BitRank OC;
    EliasFano ES, EE, ER1, ER2;
    PackedInts AB;
    u64 nocc = 0;
    const rlcsa *csa_ = nullptr;

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
                    if (ok[w]) {
                        interval I = a[w];
                        if (S.step(I, alpha[c])) { a2[c * m + w] = I; ok2[c * m + w] = 1; }
                        else a2[c * m + w].sp = (u64)i;       // absent: its longest occurring suffix is w, of length i
                    } else a2[c * m + w].sp = a[w].sp;         // or that of w
            if (narrow && i == t - 1) {                  // a holds the intervals of the (t-1)-mers
                if (!S.csa) die("--narrow needs an RLCSA");
                runs.assign(sig * m, {0, 0});
                for (u64 c = 0; c < sig; ++c)
                    for (u64 w = 0; w < m; ++w)
                        if (ok2[c * m + w]) {
                            u64 lo = S.csa->first_run(alpha[c], a[w].sp), hi = S.csa->first_run(alpha[c], a[w].ep);
                            if (hi >> 32) die("--narrow: too many runs");
                            runs[c * m + w] = {(uint32_t)lo, (uint32_t)hi};
                        }
            }
            a.swap(a2); ok.swap(ok2);
        }
        iv.swap(a); occurs.swap(ok);
        top = 1; for (int i = 1; i < t; ++i) top *= sig;
        nocc = 0; for (uint8_t x : occurs) nocc += x;
        if (compact) make_compact(S);
    }
    void make_compact(const stepper &S) {
        const u64 m = iv.size();
        OC.init(m);
        for (u64 k = 0; k < m; ++k) if (occurs[k]) OC.set(k);
        OC.finish();
        AB.init(m - nocc, 4);
        std::vector<uint64_t> st, en, r1, r2;
        st.reserve(nocc); en.reserve(nocc);
        for (u64 k = 0, ab = 0; k < m; ++k) {
            if (occurs[k]) { st.push_back(iv[k].sp); en.push_back(iv[k].ep); }
            else AB.set(ab++, iv[k].sp);
        }
        const u64 n = S.n();
        ES.build(st, n + 1); EE.build(en, n + 1);
        if (narrow) {
            csa_ = S.csa;
            r1.reserve(nocc); r2.reserve(nocc);
            for (u64 k = 0; k < m; ++k) if (occurs[k]) {
                const u64 b = csa_->run_base(alpha[k / top]);
                r1.push_back(b + runs[k].first); r2.push_back(b + runs[k].second);
            }
            ER1.build(r1, csa_->run_universe()); ER2.build(r2, csa_->run_universe());
        }
        std::vector<interval>().swap(iv); std::vector<uint8_t>().swap(occurs);
        std::vector<std::pair<uint32_t, uint32_t>>().swap(runs);
        ntable = m;
    }
    u64 ntable = 0;                 // with compact: the number of t-mers
    u64 top = 1;                    // |alpha|^(t-1)
    u64 bytes() const {
        if (compact) return OC.bytes() + ES.bytes() + EE.bytes() + AB.bytes() + ER1.bytes() + ER2.bytes();
        return iv.size() * (sizeof(interval) + 1) + runs.size() * sizeof(runs[0]);
    }
    std::string describe() const {
        char b[200];
        if (compact) snprintf(b, sizeof b, "compact: occurring bitvector %.3f, starts %.3f, ends %.3f, absent lengths %.3f, run bounds %.3f MB", OC.bytes() / 1e6, ES.bytes() / 1e6, EE.bytes() / 1e6, AB.bytes() / 1e6, (ER1.bytes() + ER2.bytes()) / 1e6);
        else snprintf(b, sizeof b, "plain: intervals %.3f, occurs %.3f, run bounds %.3f MB", iv.size() * sizeof(interval) / 1e6, occurs.size() / 1e6, runs.size() * sizeof(runs[0]) / 1e6);
        return b;
    }
    // the interval of p[0..t); false if it does not occur (or has a symbol outside the alphabet)
    bool get(const unsigned char *p, interval &I) const {
        u64 key = 0;
        for (int i = 0; i < t; ++i) { int c = code[p[i]]; if (c < 0) return false; key = alpha.size() * key + c; }
        if (compact) {
            if (!OC.get(key)) return false;
            const u64 r = OC.rank1(key);
            I = {ES.access(r), EE.access(r)};
            return true;
        }
        if (!occurs[key]) return false;
        I = iv[key];
        return true;
    }
    // as get, but on a t-mer that doesn't occur sets *miss to the length of its longest occurring suffix; -1: a symbol outside the alphabet
    int find(const unsigned char *p, interval &I, u64 &miss) const {
        u64 key = 0;
        for (int i = 0; i < t; ++i) { int c = code[p[i]]; if (c < 0) return -1; key = alpha.size() * key + c; }
        if (compact) {
            const u64 r = OC.rank1(key);
            if (!OC.get(key)) { miss = AB.get(key - r); return 0; }
            I = {ES.access(r), EE.access(r)};
            return 1;
        }
        if (!occurs[key]) { miss = iv[key].sp; return 0; }
        I = iv[key];
        return 1;
    }
    // the code of X[j..end), or its first t symbols
    void key_of(const char *X, u64 j, u64 end, table_key &K) const {
        K = table_key();
        const u64 e = end - j > (u64)t ? j + t : end;
        for (u64 p = e; p > j; --p) push(K, (unsigned char)X[p - 1]);
    }
    void push(table_key &K, unsigned char c) const {
        const int d = code[c];
        if (d < 0) { K.have = 0; K.key = 0; return; }
        K.key = K.key / alpha.size() + (u64)d * top;
        if (K.have < t) ++K.have;
    }
    // a step by c from a match with code K (narrow)
    bool step(const stepper &S, interval &I, unsigned char c, table_key &K) const {
        push(K, c);
        if (K.have < t) return S.step(I, c);
        if (compact) {
            if (!OC.get(K.key)) return false;
            const u64 r = OC.rank1(K.key), b = csa_->run_base(c);
            return S.step_in(I, c, ER1.access(r) - b, ER2.access(r) - b);
        }
        if (!occurs[K.key]) return false;
        return S.step_in(I, c, runs[K.key].first, runs[K.key].second);
    }
};

}  // namespace katka
