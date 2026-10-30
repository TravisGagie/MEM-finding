// Elias-Fano coded run lists for the compressed RLCSA (rlcsa.hpp), after Brown, Gagie, Manzini, Navarro and
// Sciortino, "Faster run-length compressed suffix arrays".  Only a variant of the explicit run lists; the
// paper's default configuration uses the explicit ones.
//
// The runs of the BWT are listed grouped by character and, within a character, by start.  Two Elias-Fano
// sequences over them:
//   KS[j] = c * n + start_j  (increasing): a predecessor query finds the last run of c starting before i;
//   FP[j] = position in F (the sorted BWT) of the first character of run j.
// rank_c(i) = FP[j] - FP[first run of c] + min(length of run j, i - start_j), with j that last run.
#pragma once
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <vector>
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

template <class V> inline void write_vec(FILE *f, const V &v) {
    uint64_t k = v.size(); fwrite(&k, 8, 1, f);
    if (k) fwrite(v.data(), sizeof(v[0]), k, f);
}
template <class V> inline void read_vec(FILE *f, V &v) {
    uint64_t k;
    if (fread(&k, 8, 1, f) != 1) throw std::runtime_error("read");
    v.resize(k);
    if (k && fread(v.data(), sizeof(v[0]), k, f) != k) throw std::runtime_error("read");
}

// a non-decreasing sequence of n integers below u: l low bits each, and the high parts in unary
struct EliasFano {
    uint64_t n = 0, u = 0;
    uint32_t l = 0;
    std::vector<uint64_t> low, high;
    uint64_t hbits = 0;
    std::vector<uint64_t> s0, s1;           // positions of every SAMPLE-th zero and one in high
    static constexpr uint64_t SAMPLE = 256;

    void build(const std::vector<uint64_t> &v, uint64_t universe) {
        n = v.size(); u = universe; l = 0;
        if (n && u / n > 1) while ((u / n) >> (l + 1)) ++l;
        low.assign((n * l + 63) / 64 + 1, 0);
        hbits = n + (u >> l) + 2;
        high.assign(hbits / 64 + 2, 0);
        for (uint64_t j = 0; j < n; ++j) {
            if (j && v[j] < v[j - 1]) throw std::runtime_error("Elias-Fano: not sorted");
            if (l) set_low(j, v[j] & ((1ull << l) - 1));
            uint64_t p = (v[j] >> l) + j;
            high[p >> 6] |= 1ull << (p & 63);
        }
        s0.clear(); s1.clear();
        uint64_t c0 = 0, c1 = 0;
        for (uint64_t p = 0; p < hbits; ++p) {
            if (hbit(p)) { if (c1 % SAMPLE == 0) s1.push_back(p); ++c1; }
            else { if (c0 % SAMPLE == 0) s0.push_back(p); ++c0; }
        }
    }
    void set_low(uint64_t j, uint64_t x) {
        uint64_t b = j * l;
        low[b >> 6] |= x << (b & 63);
        if ((b & 63) + l > 64) low[(b >> 6) + 1] |= x >> (64 - (b & 63));
    }
    uint64_t get_low(uint64_t j) const {
        if (!l) return 0;
        uint64_t b = j * l, w = b >> 6, o = b & 63, x = low[w] >> o;
        if (o + l > 64) x |= low[w + 1] << (64 - o);
        return x & ((1ull << l) - 1);
    }
    bool hbit(uint64_t p) const { return (high[p >> 6] >> (p & 63)) & 1; }

    uint64_t select(uint64_t k, bool one) const {       // position of the k-th (0-based) one or zero in high
        uint64_t p = (one ? s1 : s0)[k / SAMPLE], rem = k % SAMPLE, w = p >> 6;
        uint64_t word = (one ? high[w] : ~high[w]) & (~0ull << (p & 63));
        for (;;) {
            uint64_t c = __builtin_popcountll(word);
            if (rem < c) return (w << 6) + select_in_word(word, rem);
            rem -= c; ++w; word = one ? high[w] : ~high[w];
        }
    }
    uint64_t access(uint64_t j) const { return ((select(j, true) - j) << l) | get_low(j); }
    void access2(uint64_t j, uint64_t &a, uint64_t &b) const {    // access(j) and access(j + 1), j + 1 < n
        uint64_t p = select(j, true), q = p + 1;
        a = ((p - j) << l) | get_low(j);
        while (!hbit(q)) ++q;
        b = ((q - j - 1) << l) | get_low(j + 1);
    }
    uint64_t count_less(uint64_t x) const {             // number of elements < x
        if (x >= u) return n;
        uint64_t h = x >> l, xl = x & ((1ull << l) - 1), idx = 0, pos = 0;
        if (h) { uint64_t p = select(h - 1, false); idx = p - (h - 1); pos = p + 1; }
        while (pos < hbits && hbit(pos) && get_low(idx) < xl) { ++idx; ++pos; }
        return idx;
    }
    uint64_t bytes() const { return (low.size() + high.size() + s0.size() + s1.size()) * 8 + 40; }
    void save(FILE *f) const {
        fwrite(&n, 8, 1, f); fwrite(&u, 8, 1, f); fwrite(&l, 4, 1, f); fwrite(&hbits, 8, 1, f);
        write_vec(f, low); write_vec(f, high); write_vec(f, s0); write_vec(f, s1);
    }
    void load(FILE *f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&u, 8, 1, f) != 1 || fread(&l, 4, 1, f) != 1 || fread(&hbits, 8, 1, f) != 1)
            throw std::runtime_error("read");
        read_vec(f, low); read_vec(f, high); read_vec(f, s0); read_vec(f, s1);
    }
};

struct EFRank {
    uint64_t n = 0;
    uint32_t sigma = 0;
    uint64_t nruns = 0;
    EliasFano KS, FP;

    // rank_c(i) when b (c's first run) and a (number of runs of c starting before i, plus b) are known
    uint64_t rank_from(uint64_t b, uint64_t a, uint32_t c, uint64_t i) const {
        uint64_t j = a - 1, start = KS.access(j) - (uint64_t)c * n, f, fn;
        if (j + 1 < nruns) FP.access2(j, f, fn); else { f = FP.access(j); fn = n; }
        uint64_t len = fn - f, d = i - start;
        return f - FP.access(b) + (d < len ? d : len);
    }
    // rank_c at two positions, finding c's first run once
    void rank2(uint32_t c, uint64_t i1, uint64_t i2, uint64_t &r1, uint64_t &r2) const {
        if (c >= sigma) { r1 = r2 = 0; return; }
        uint64_t base = (uint64_t)c * n, b = KS.count_less(base);
        uint64_t a1 = i1 ? KS.count_less(base + i1) : b, a2 = i2 ? KS.count_less(base + i2) : b;
        r1 = a1 == b ? 0 : rank_from(b, a1, c, i1);
        r2 = a2 == b ? 0 : rank_from(b, a2, c, i2);
    }
    uint64_t bytes() const { return KS.bytes() + FP.bytes(); }
    void save(FILE *f) const { fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nruns, 8, 1, f); KS.save(f); FP.save(f); }
    void load(FILE *f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&nruns, 8, 1, f) != 1) throw std::runtime_error("read");
        KS.load(f); FP.load(f);
    }
};

}  // namespace katka
