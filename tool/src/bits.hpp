// Succinct building blocks, with no outside dependencies:
//   EliasFano   a non-decreasing integer sequence (equivalently a sparse bitvector): access, count_less
//               (= rank), with a rank directory over its high bits for fast selects;
//   PackedInts  fixed-width integers;
//   BitRank     a plain bitvector with rank and select;
//   RLBits      a bitvector stored as runs of ones;
//   EFRank      rank over a sequence with a large alphabet, by Elias-Fano lists of its runs grouped by
//               symbol (Brown et al., "Faster run-length compressed suffix arrays"): the RLCSA of
//               Section 2 of the paper, used for the BWTs of the parse and of the block sequence;
//   WaveletMatrix  rank and access over a small alphabet (the run heads of the RLBWT).
#pragma once
#include "common.hpp"
#include <algorithm>
#ifdef __BMI2__
#include <immintrin.h>
#endif

namespace katka {

inline uint64_t select_in_word(uint64_t w, uint64_t k) {      // position of the k-th (0-based) one in w
#ifdef __BMI2__
    return __builtin_ctzll(_pdep_u64(1ull << k, w));
#else
    for (uint64_t i = 0;; ++i) if ((w >> i) & 1) { if (!k) return i; --k; }
#endif
}

// A non-decreasing sequence of n integers below u: l low bits each, the high parts in unary.
struct EliasFano {
    uint64_t n = 0, u = 0;
    uint32_t l = 0;
    std::vector<uint64_t> low, high;
    uint64_t hbits = 0;
    std::vector<uint64_t> s0, s1;       // positions of every SAMPLE-th zero and one in high
    std::vector<uint64_t> ones;         // ones before each 512-bit block of high
    static constexpr uint64_t SAMPLE = 256;
#ifndef KATKA_EF_SAMPLE1
#define KATKA_EF_SAMPLE1 64
#endif
    // finer samples for select-one, built at load (not saved): the position of every SAMPLE1-th one, 32 bits each
    std::vector<uint32_t> f1;
    static constexpr uint64_t SAMPLE1 = KATKA_EF_SAMPLE1;
    void build_fine() {
        f1.clear();
        if (!SAMPLE1 || hbits >= (1ull << 32)) return;
        f1.reserve(n / SAMPLE1 + 1);
        uint64_t c1 = 0;
        for (uint64_t w = 0; w < high.size(); ++w) {
            uint64_t word = high[w], pc = __builtin_popcountll(word);
            while (f1.size() * SAMPLE1 < c1 + pc) {
                f1.push_back((uint32_t)((w << 6) + select_in_word(word, f1.size() * SAMPLE1 - c1)));
            }
            c1 += pc;
        }
    }

    void build(const std::vector<uint64_t> &v, uint64_t universe) {
        begin(v.size(), universe);
        for (uint64_t x : v) push(x);
        finish();
    }
    // streaming: begin(count, universe), push the values in order, finish()
    uint64_t pushed = 0, lastv = 0;
    void begin(uint64_t count, uint64_t universe) {
        n = count; u = universe; l = 0; pushed = 0; lastv = 0;
        if (n && u / n > 1) while ((u / n) >> (l + 1)) ++l;
        low.assign((n * l + 63) / 64 + 1, 0);
        hbits = n + (u >> l) + 2;
        high.assign(hbits / 64 + 2, 0);
    }
    inline void push(uint64_t x) {
        const uint64_t j = pushed++;
        if (j >= n) throw std::runtime_error("Elias-Fano: too many values");
        if (j && x < lastv) throw std::runtime_error("Elias-Fano: not sorted");
        if (x >= u && u) throw std::runtime_error("Elias-Fano: value out of range");
        lastv = x;
        if (l) set_low(j, x & ((1ull << l) - 1));
        uint64_t p = (x >> l) + j;
        high[p >> 6] |= 1ull << (p & 63);
    }
    void finish() {
        if (pushed != n) throw std::runtime_error("Elias-Fano: too few values");
        ones.assign(hbits / 512 + 2, 0);
        for (uint64_t b = 1; b < ones.size(); ++b) {
            uint64_t c = 0;
            for (uint64_t w = (b - 1) * 8; w < b * 8 && w < high.size(); ++w) c += __builtin_popcountll(high[w]);
            ones[b] = ones[b - 1] + c;
        }
        s0.clear(); s1.clear();
        uint64_t c0 = 0, c1 = 0;
        for (uint64_t p = 0; p < hbits; ++p) {
            if (hbit(p)) { if (c1 % SAMPLE == 0) s1.push_back(p); ++c1; }
            else { if (c0 % SAMPLE == 0) s0.push_back(p); ++c0; }
        }
        build_fine();
    }
    inline void set_low(uint64_t j, uint64_t x) {
        uint64_t b = j * l;
        low[b >> 6] |= x << (b & 63);
        if ((b & 63) + l > 64) low[(b >> 6) + 1] |= x >> (64 - (b & 63));
    }
    inline uint64_t get_low(uint64_t j) const {
        if (!l) return 0;
        uint64_t b = j * l, w = b >> 6, o = b & 63, x = low[w] >> o;
        if (o + l > 64) x |= low[w + 1] << (64 - o);
        return x & ((1ull << l) - 1);
    }
    inline bool hbit(uint64_t p) const { return (high[p >> 6] >> (p & 63)) & 1; }

    template <bool ONE> inline uint64_t select_(uint64_t k) const {
        if (ONE && !f1.empty()) {
            uint64_t i = k / SAMPLE1, p = f1[i], rem = k - i * SAMPLE1, w = p >> 6;
            uint64_t word = high[w] & (~0ull << (p & 63));
            for (;;) {
                uint64_t c = __builtin_popcountll(word);
                if (rem < c) return (w << 6) + select_in_word(word, rem);
                rem -= c; word = high[++w];
            }
        }
        const std::vector<uint64_t> &smp = ONE ? s1 : s0;
        uint64_t i = k / SAMPLE;
        uint64_t lo = smp[i] >> 9;
        uint64_t hi = (i + 1 < smp.size()) ? (smp[i + 1] >> 9) + 1 : ones.size() - 1;
        auto before = [&](uint64_t blk) { return ONE ? ones[blk] : (blk << 9) - ones[blk]; };
        while (hi - lo > 1) {
            uint64_t mid = (lo + hi) >> 1;
            if (before(mid) <= k) lo = mid; else hi = mid;
        }
        uint64_t rem = k - before(lo);
        for (uint64_t w = lo * 8;; ++w) {
            uint64_t word = ONE ? high[w] : ~high[w];
            uint64_t c = __builtin_popcountll(word);
            if (rem < c) return (w << 6) + select_in_word(word, rem);
            rem -= c;
        }
    }
    inline uint64_t access(uint64_t j) const { return ((select_<true>(j) - j) << l) | get_low(j); }
    // access(j) and access(j + 1), j + 1 < n
    inline void access2(uint64_t j, uint64_t &a, uint64_t &b) const {
        uint64_t p = select_<true>(j), q = p + 1;
        a = ((p - j) << l) | get_low(j);
        while (!hbit(q)) ++q;
        b = ((q - j - 1) << l) | get_low(j + 1);
    }
    inline uint64_t next_one(uint64_t q) const {             // first one at or after position q
        uint64_t w = q >> 6, word = high[w] & (~0ull << (q & 63));
        while (!word) word = high[++w];
        return (w << 6) + __builtin_ctzll(word);
    }
    inline uint64_t next_zero(uint64_t q) const {
        uint64_t w = q >> 6, word = ~high[w] & (~0ull << (q & 63));
        while (!word) word = ~high[++w];
        return (w << 6) + __builtin_ctzll(word);
    }
    // the indices [lo, hi) whose high part is h
    inline void bucket(uint64_t h, uint64_t &lo, uint64_t &hi) const {
        if (h == 0) { lo = 0; hi = select_<false>(0); }
        else {                                   // the bucket ends at the next zero after the (h-1)-th
            uint64_t p = select_<false>(h - 1);
            lo = p - (h - 1);
            hi = next_zero(p + 1) - h;
        }
    }
    inline uint64_t in_bucket(uint64_t lo, uint64_t hi, uint64_t xl) const {
        while (hi - lo > 8) {
            uint64_t mid = (lo + hi) >> 1;
            if (get_low(mid) < xl) lo = mid + 1; else hi = mid;
        }
        while (lo < hi && get_low(lo) < xl) ++lo;
        return lo;
    }
    // number of elements < x
    inline uint64_t count_less(uint64_t x) const {
        if (x >= u) return n;
        uint64_t lo, hi; bucket(x >> l, lo, hi);
        return in_bucket(lo, hi, x & ((1ull << l) - 1));
    }
    // count_less at x1 <= x2: one select-zero when their high parts are equal or adjacent
    inline void count_less2(uint64_t x1, uint64_t x2, uint64_t &r1, uint64_t &r2) const {
        const uint64_t h1 = x1 >> l, h2 = x2 >> l, m = (1ull << l) - 1;
        if (x2 >= u || h2 > h1 + 1) { r1 = count_less(x1); r2 = count_less(x2); return; }
        uint64_t lo, hi; bucket(h1, lo, hi);
        r1 = in_bucket(lo, hi, x1 & m);
        if (h2 == h1) { r2 = in_bucket(r1, hi, x2 & m); return; }
        r2 = in_bucket(hi, next_zero(hi + h1 + 1) - h2, x2 & m);   // bucket h1 ends with the zero at hi + h1
    }
    // access(j), access(j + 1), access(j + 2) (b, c = u if past the end): one select
    inline void access3(uint64_t j, uint64_t &a, uint64_t &b, uint64_t &c) const {
        uint64_t p = select_<true>(j);
        a = ((p - j) << l) | get_low(j);
        if (j + 1 >= n) { b = c = u; return; }
        p = next_one(p + 1); b = ((p - j - 1) << l) | get_low(j + 1);
        if (j + 2 >= n) { c = u; return; }
        p = next_one(p + 1); c = ((p - j - 2) << l) | get_low(j + 2);
    }
    // as a bitvector over [0, u) with ones at the values (distinct values)
    inline uint64_t rank1(uint64_t i) const { return count_less(i); }      // ones in [0, i)
    inline uint64_t select1(uint64_t k) const { return access(k); }         // position of the k-th one, 0-based
    inline uint64_t pred_index(uint64_t x) const { return count_less(x + 1) - 1; }   // last element <= x (exists)

    uint64_t bytes() const { return (low.size() + high.size() + s0.size() + s1.size() + ones.size()) * 8 + f1.size() * 4 + 40; }
    void save(FILE *f) const {
        wr(f, n); wr(f, u); wr(f, l); wr(f, hbits);
        wvec(f, low); wvec(f, high); wvec(f, s0); wvec(f, s1); wvec(f, ones);
    }
    void load(FILE *f) {
        rd(f, n); rd(f, u); rd(f, l); rd(f, hbits);
        rvec(f, low); rvec(f, high); rvec(f, s0); rvec(f, s1); rvec(f, ones);
        build_fine();
    }
};

// Fixed-width packed integers.
struct PackedInts {
    uint64_t n = 0;
    uint32_t w = 1;
    std::vector<uint64_t> d;
    void init(uint64_t count, uint32_t width) { n = count; w = std::max<uint32_t>(1, width); d.assign((n * w + 63) / 64 + 1, 0); }
    static uint32_t width_for(uint64_t mx) { uint32_t w = 1; while (w < 64 && (mx >> w)) ++w; return w; }
    void build(const std::vector<uint64_t> &v) {
        uint64_t mx = 0;
        for (auto x : v) mx = std::max(mx, x);
        init(v.size(), width_for(mx));
        for (uint64_t i = 0; i < n; ++i) set(i, v[i]);
    }
    template <class V> void build_from(const V &v) {
        uint64_t mx = 0;
        for (auto x : v) mx = std::max<uint64_t>(mx, x);
        init(v.size(), width_for(mx));
        for (uint64_t i = 0; i < n; ++i) set(i, v[i]);
    }
    inline void set(uint64_t i, uint64_t x) {
        uint64_t b = i * w, k = b >> 6, o = b & 63;
        uint64_t mask = w == 64 ? ~0ull : (1ull << w) - 1;
        x &= mask;
        d[k] = (d[k] & ~(mask << o)) | (x << o);
        if (o + w > 64) { uint64_t hm = mask >> (64 - o); d[k + 1] = (d[k + 1] & ~hm) | (x >> (64 - o)); }
    }
    inline uint64_t get(uint64_t i) const {
        uint64_t b = i * w, k = b >> 6, o = b & 63, x = d[k] >> o;
        if (o + w > 64) x |= d[k + 1] << (64 - o);
        return w == 64 ? x : x & ((1ull << w) - 1);
    }
    uint64_t operator[](uint64_t i) const { return get(i); }
    uint64_t size() const { return n; }
    uint64_t bytes() const { return d.size() * 8 + 12; }
    void save(FILE *f) const { wr(f, n); wr(f, w); wvec(f, d); }
    void load(FILE *f) { rd(f, n); rd(f, w); rvec(f, d); }
};

// A plain bitvector with rank (counts every 512 bits, and within those every 64) and select (sampled).
struct BitRank {
    uint64_t n = 0, nones = 0;
    std::vector<uint64_t> bits, cnt, s1, s0;
    std::vector<uint16_t> sub;             // ones before each word, within its 512-bit block
    static constexpr uint64_t SAMPLE = 512;
    void init(uint64_t size) { n = size; bits.assign(n / 64 + 2, 0); }
    inline void set(uint64_t i) { bits[i >> 6] |= 1ull << (i & 63); }
    inline bool get(uint64_t i) const { return (bits[i >> 6] >> (i & 63)) & 1; }
    void finish() {                       // after all set()s
        cnt.assign(bits.size() / 8 + 2, 0);
        uint64_t c = 0;
        sub.assign(bits.size(), 0);
        for (uint64_t w = 0; w < bits.size(); ++w) {
            if (w % 8 == 0) cnt[w / 8] = c;
            sub[w] = (uint16_t)(c - cnt[w / 8]);
            c += __builtin_popcountll(bits[w]);
        }
        cnt[bits.size() / 8 + 1] = c;
        if (bits.size() % 8 == 0) cnt[bits.size() / 8] = c;
        nones = rank1(n);
        s1.clear(); s0.clear();
        uint64_t c1 = 0, c0 = 0;
        for (uint64_t i = 0; i < n; ++i) {
            if (get(i)) { if (c1++ % SAMPLE == 0) s1.push_back(i >> 9); }
            else { if (c0++ % SAMPLE == 0) s0.push_back(i >> 9); }
        }
    }
    inline uint64_t rank1(uint64_t i) const {          // ones in [0, i)
        uint64_t w = i >> 6, r = cnt[w >> 3] + sub[w];
        if (i & 63) r += __builtin_popcountll(bits[w] & ((1ull << (i & 63)) - 1));
        return r;
    }
    inline uint64_t rank0(uint64_t i) const { return i - rank1(i); }
    template <bool ONE> uint64_t select_(uint64_t k) const {   // position of the k-th (0-based) one / zero
        const std::vector<uint64_t> &smp = ONE ? s1 : s0;
        uint64_t lo = smp[k / SAMPLE], hi = k / SAMPLE + 1 < smp.size() ? smp[k / SAMPLE + 1] + 1 : cnt.size() - 1;
        auto before = [&](uint64_t b) { return ONE ? cnt[b] : (b << 9) - cnt[b]; };
        while (hi - lo > 1) { uint64_t m = (lo + hi) >> 1; if (before(m) <= k) lo = m; else hi = m; }
        uint64_t rem = k - before(lo);
        for (uint64_t w = lo * 8;; ++w) {
            uint64_t word = ONE ? bits[w] : ~bits[w];
            uint64_t c = __builtin_popcountll(word);
            if (rem < c) return (w << 6) + select_in_word(word, rem);
            rem -= c;
        }
    }
    inline uint64_t select1(uint64_t k) const { return select_<true>(k); }
    inline uint64_t select0(uint64_t k) const { return select_<false>(k); }
    uint64_t bytes() const { return (bits.size() + cnt.size() + s0.size() + s1.size()) * 8 + sub.size() * 2 + 16; }
    void save(FILE *f) const { wr(f, n); wr(f, nones); wvec(f, bits); wvec(f, cnt); wvec(f, s1); wvec(f, s0); wvec(f, sub); }
    void load(FILE *f) { rd(f, n); rd(f, nones); rvec(f, bits); rvec(f, cnt); rvec(f, s1); rvec(f, s0); rvec(f, sub); }
};

// A bitvector stored as runs of ones: rank and select by Elias-Fano lists of the run starts and of the
// number of ones before each run.
struct RLBits {
    uint64_t n = 0, ones = 0, nruns = 0;
    EliasFano RS, RC;
    void build(const std::vector<std::pair<uint64_t, uint64_t>> &runs, uint64_t size) {   // (start, length), sorted
        n = size; nruns = runs.size();
        std::vector<uint64_t> s(nruns), c(nruns + 1);
        uint64_t acc = 0;
        for (uint64_t k = 0; k < nruns; ++k) { s[k] = runs[k].first; c[k] = acc; acc += runs[k].second; }
        c[nruns] = acc; ones = acc;
        RS.build(s, n + 1); RC.build(c, ones + 1);
    }
    inline uint64_t rank1(uint64_t pos) const {        // ones in [0, pos)
        uint64_t k = RS.count_less(pos);
        if (k == 0) return 0;
        --k;
        uint64_t s = RS.access(k), c0, c1;
        RC.access2(k, c0, c1);
        uint64_t d = pos - s;
        return c0 + (d < c1 - c0 ? d : c1 - c0);
    }
    inline uint64_t select1(uint64_t r) const {        // position of the r-th one (0-based)
        uint64_t k = RC.count_less(r + 1) - 1;
        return RS.access(k) + (r - RC.access(k));
    }
    uint64_t bytes() const { return RS.bytes() + RC.bytes() + 24; }
    void save(FILE *f) const { wr(f, n); wr(f, ones); wr(f, nruns); RS.save(f); RC.save(f); }
    void load(FILE *f) { rd(f, n); rd(f, ones); rd(f, nruns); RS.load(f); RC.load(f); }
};

// rank over a sequence s[0, n) of symbols below sigma, by its runs (see the top of the file)
struct EFRank {
    uint64_t n = 0;
    uint32_t sigma = 0;
    uint64_t nruns = 0;
    EliasFano KS, FP;     // KS[j] = c * n + start of run j (runs grouped by c); FP[j] = its position in F
    // For an RLCSA (rlcsa::build_ef): slot[c], the rank of c among the characters that occur, so that KS's universe is (characters + 1) n,
    // not 256 n, and its low parts are about 6 bits shorter (the buckets of the predecessor search hold 1 or 2 runs, not about 80); a
    // character that doesn't occur maps past the last.  Empty: c itself.  first[c] and ffirst[c]: c's first run in KS and its position
    // in F, for every c (they depend only on c).
    std::vector<uint32_t> slot;
    std::vector<uint64_t> first, ffirst;
    inline uint64_t base_of(uint32_t c) const { return (uint64_t)(slot.empty() ? c : slot[c]) * n; }
    void prepare() {                      // after build or load: the per-character values
        if (sigma > 65536) { first.clear(); ffirst.clear(); return; }
        first.assign(sigma, 0); ffirst.assign(sigma, 0);
        for (uint32_t c = 0; c < sigma; ++c) { first[c] = KS.count_less(base_of(c)); ffirst[c] = first[c] < nruns ? FP.access(first[c]) : n; }
    }

    template <class V> void build(const V &s, uint32_t sig) {
        n = s.size(); sigma = sig;
        if ((unsigned __int128)sigma * (n + 1) >= ((unsigned __int128)1 << 63)) throw std::runtime_error("sigma*n too large");
        std::vector<uint64_t> cnt(sigma + 1, 0), runsof(sigma + 1, 0);
        nruns = 0;
        for (uint64_t i = 0; i < n; ++i) {
            cnt[s[i]]++;
            if (i == 0 || s[i] != s[i - 1]) { runsof[s[i]]++; nruns++; }
        }
        std::vector<uint64_t> Cf(sigma + 1, 0), slot(sigma + 1, 0);
        for (uint32_t c = 0; c < sigma; ++c) { Cf[c + 1] = Cf[c] + cnt[c]; slot[c + 1] = slot[c] + runsof[c]; }
        std::vector<uint64_t> ks(nruns), fp(nruns), seen(sigma, 0);
        for (uint64_t i = 0; i < n;) {
            uint64_t j = i + 1;
            while (j < n && s[j] == s[i]) ++j;
            uint32_t c = s[i];
            uint64_t k = slot[c]++;
            ks[k] = (uint64_t)c * n + i; fp[k] = Cf[c] + seen[c]; seen[c] += j - i;
            i = j;
        }
        KS.build(ks, (uint64_t)sigma * n);
        FP.build(fp, n);
    }
    inline uint64_t rank_from(uint64_t fb, uint64_t a, uint64_t base, uint64_t i) const {   // fb: c's first run's position in F
        uint64_t j = a - 1, start = KS.access(j) - base, f, fn;
        if (j + 1 < nruns) FP.access2(j, f, fn); else { f = FP.access(j); fn = n; }
        uint64_t len = fn - f, d = i - start;
        return f - fb + (d < len ? d : len);
    }
#ifndef KATKA_EF_PAIR
#define KATKA_EF_PAIR 1
#endif
    // r1, r2 for positions i1 <= i2 whose runs are a1 - 1 and a2 - 1 (a1, a2 > b): when the runs are the same or
    // neighbours, one select in KS and one in FP serve both
    inline void rank_pair(uint64_t fb, uint64_t a1, uint64_t a2, uint64_t base, uint64_t i1, uint64_t i2, uint64_t &r1, uint64_t &r2) const {
        if (KATKA_EF_PAIR && a2 == a1) {
            uint64_t j = a1 - 1, start = KS.access(j) - base, f, fn;
            if (j + 1 < nruns) FP.access2(j, f, fn); else { f = FP.access(j); fn = n; }
            uint64_t len = fn - f, d1 = i1 - start, d2 = i2 - start;
            r1 = f - fb + (d1 < len ? d1 : len); r2 = f - fb + (d2 < len ? d2 : len);
        } else if (KATKA_EF_PAIR && a2 == a1 + 1) {
            uint64_t j = a1 - 1, s1, s2, f1, f2, f3, x;
            KS.access3(j, s1, s2, x); FP.access3(j, f1, f2, f3);
            if (j + 2 >= nruns) f3 = n;
            uint64_t d1 = i1 - (s1 - base), d2 = i2 - (s2 - base), l1 = f2 - f1, l2 = f3 - f2;
            r1 = f1 - fb + (d1 < l1 ? d1 : l1); r2 = f2 - fb + (d2 < l2 ? d2 : l2);
        } else { r1 = rank_from(fb, a1, base, i1); r2 = rank_from(fb, a2, base, i2); }
    }
    // rank_c at two positions, finding c's first run once
    inline void rank2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t &r1, uint64_t &r2) const {
        if (c >= sigma) { r1 = r2 = 0; return; }
        const uint64_t base = base_of(c), b = first.empty() ? KS.count_less(base) : first[c];
        uint64_t a1, a2;
        if (KATKA_EF_PAIR && i1 && i1 <= i2) KS.count_less2(base + i1, base + i2, a1, a2);
        else { a1 = i1 ? KS.count_less(base + i1) : b; a2 = i2 ? KS.count_less(base + i2) : b; }
        if (a1 == b && a2 == b) { r1 = r2 = 0; return; }
        const uint64_t fb = first.empty() ? FP.access(b) : ffirst[c];
        if (a1 == b) { r1 = 0; r2 = rank_from(fb, a2, base, i2); }
        else if (i1 <= i2) rank_pair(fb, a1, a2, base, i1, i2, r1, r2);
        else { r1 = rank_from(fb, a1, base, i1); r2 = rank_from(fb, a2, base, i2); }
    }
    // rank2, knowing that both positions' runs lie among runs rl..rh (positions in KS): binary search over KS within them
    inline uint64_t count_less_in(uint64_t x, uint64_t lo, uint64_t hi) const {
        while (lo < hi) {
            uint64_t mid = (lo + hi) >> 1;
            if (KS.access(mid) < x) lo = mid + 1; else hi = mid;
        }
        return lo;
    }
    inline void rank2_in(uint32_t c, uint64_t i1, uint64_t i2, uint64_t rl, uint64_t rh, uint64_t &r1, uint64_t &r2) const {
        if (c >= sigma) { r1 = r2 = 0; return; }
        const uint64_t base = base_of(c), b = first.empty() ? KS.count_less(base) : first[c];
        uint64_t a1 = i1 ? count_less_in(base + i1, rl, rh) : b, a2 = i2 ? count_less_in(base + i2, rl, rh) : b;
        if (a1 == b && a2 == b) { r1 = r2 = 0; return; }
        const uint64_t fb = first.empty() ? FP.access(b) : ffirst[c];
        if (a1 == b) { r1 = 0; r2 = rank_from(fb, a2, base, i2); }
        else if (i1 <= i2) rank_pair(fb, a1, a2, base, i1, i2, r1, r2);
        else { r1 = rank_from(fb, a1, base, i1); r2 = rank_from(fb, a2, base, i2); }
    }
    // LF without a C array: C[c] + rank_c(i)
    inline uint64_t lf_one(uint64_t base, uint64_t b, uint64_t i) const {
        uint64_t a = i ? KS.count_less(base + (i > n ? n : i)) : b;
        if (a == b) return b < nruns ? FP.access(b) : n;
        uint64_t j = a - 1, start = KS.access(j) - base, f, fn;
        if (j + 1 < nruns) FP.access2(j, f, fn); else { f = FP.access(j); fn = n; }
        uint64_t len = fn - f, d = i - start;
        return f + (d < len ? d : len);
    }
    // a backward step: [i1, i2) -> [o1, o2) (empty if o1 >= o2)
    inline void lf2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t &o1, uint64_t &o2) const {
        if (c >= sigma) { o1 = o2 = n; return; }
        uint64_t base = (uint64_t)c * n, b = KS.count_less(base);
        if (b == nruns || KS.access(b) >= base + n) { o1 = o2 = n; return; }    // c does not occur
        o1 = lf_one(base, b, i1); o2 = lf_one(base, b, i2);
    }
    uint64_t bytes() const { return KS.bytes() + FP.bytes() + slot.size() * 4 + (first.size() + ffirst.size()) * 8; }
    void save(FILE *f) const { wr(f, n); wr(f, sigma); wr(f, nruns); KS.save(f); FP.save(f); }
    void load(FILE *f) { rd(f, n); rd(f, sigma); rd(f, nruns); KS.load(f); FP.load(f); }
};

// Wavelet matrix over symbols 0 .. sigma-1 (sigma small): access and rank.
struct WaveletMatrix {
    uint64_t n = 0;
    uint32_t levels = 0;
    std::vector<BitRank> bv;
    std::vector<uint64_t> zeros;
    std::vector<uint64_t> bstart;     // where the symbol's positions start after the last level

    template <class V> void build(const V &s, uint32_t sigma) {
        n = s.size(); levels = 1;
        while ((1u << levels) < sigma) ++levels;
        bv.assign(levels, BitRank()); zeros.assign(levels, 0);
        std::vector<uint32_t> cur(s.begin(), s.end()), nxt(n);
        for (uint32_t lv = 0; lv < levels; ++lv) {
            uint32_t sh = levels - 1 - lv;
            bv[lv].init(n);
            uint64_t z = 0;
            for (uint64_t i = 0; i < n; ++i) if (!((cur[i] >> sh) & 1)) ++z;
            zeros[lv] = z;
            uint64_t a = 0, b = z;
            for (uint64_t i = 0; i < n; ++i) {
                if ((cur[i] >> sh) & 1) { bv[lv].set(i); nxt[b++] = cur[i]; } else nxt[a++] = cur[i];
            }
            bv[lv].finish();
            cur.swap(nxt);
        }
        set_bstart();
    }
    void set_bstart() {
        bstart.assign(1u << levels, 0);
        for (uint32_t c = 0; c < (1u << levels); ++c) {
            uint64_t s = 0;
            for (uint32_t lv = 0; lv < levels; ++lv)
                s = ((c >> (levels - 1 - lv)) & 1) ? zeros[lv] + bv[lv].rank1(s) : bv[lv].rank0(s);
            bstart[c] = s;
        }
    }
    inline uint32_t access(uint64_t i) const {
        uint32_t c = 0;
        for (uint32_t lv = 0; lv < levels; ++lv) {
            const BitRank &b = bv[lv];
            if (b.get(i)) { c = (c << 1) | 1; i = zeros[lv] + b.rank1(i); }
            else { c <<= 1; i = b.rank0(i); }
        }
        return c;
    }
    // occurrences of c in [0, i)
    inline uint64_t rank(uint32_t c, uint64_t i) const {
        for (uint32_t lv = 0; lv < levels; ++lv) {
            const BitRank &b = bv[lv];
            i = ((c >> (levels - 1 - lv)) & 1) ? zeros[lv] + b.rank1(i) : b.rank0(i);
        }
        return i - bstart[c];
    }
    // (rank of s[i] in [0, i), s[i]) in one traversal
    inline std::pair<uint64_t, uint32_t> inverse_select(uint64_t i) const {
        uint32_t c = 0;
        for (uint32_t lv = 0; lv < levels; ++lv) {
            const BitRank &b = bv[lv];
            if (b.get(i)) { c = (c << 1) | 1; i = zeros[lv] + b.rank1(i); }
            else { c <<= 1; i = b.rank0(i); }
        }
        return {i - bstart[c], c};
    }
    uint64_t bytes() const { uint64_t b = 16 + zeros.size() * 8; for (auto &x : bv) b += x.bytes(); return b; }
    void save(FILE *f) const { wr(f, n); wr(f, levels); wvec(f, zeros); for (auto &x : bv) x.save(f); }
    void load(FILE *f) { rd(f, n); rd(f, levels); rvec(f, zeros); bv.assign(levels, BitRank()); for (auto &x : bv) x.load(f); set_bstart(); }
};

}  // namespace katka
