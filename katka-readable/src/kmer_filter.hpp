// The KeBaB k-mer filter (Section 5, "Longest MEMs"): a blocked Bloom filter of the canonical k-mers
// (k <= 32) of the text.  A k-mer of a read that the filter rejects occurs in no reference, so no MEM
// contains it; the maximal substrings of the read all of whose k-mers pass the filter (pseudo-MEMs) contain
// every MEM of length at least k.
//
// All h bits of a k-mer lie in one 512-bit block (one cache line), so a query costs one memory access
// whatever h is, which allows the optimal number of hash functions.
//
// File (katka-kbbuild): "RZKBLM01", k, h (32-bit), number of blocks (64-bit), the blocks.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace katka {

struct kmer_filter {
    uint32_t k = 0, h = 0;
    uint64_t nblocks = 0;
    std::vector<uint64_t> bits;                     // nblocks blocks of 8 words

    static uint64_t mix(uint64_t x) {               // the splitmix64 finalizer
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL; x ^= x >> 27; x *= 0x94d049bb133111ebULL; x ^= x >> 31; return x;
    }
    static int code(unsigned char c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; } }
    uint64_t block_of(uint64_t hv) const { return (uint64_t)(((unsigned __int128)hv * nblocks) >> 64); }

    // the h bit positions of a k-mer, all in its block, 9 bits each from a stream of hashes
    template <class F> void for_each_bit(uint64_t kmer, F f) const {
        uint64_t hv = mix(kmer), g = mix(hv ^ 0x9e3779b97f4a7c15ULL), blk = block_of(hv) * 8;
        for (uint32_t i = 0; i < h; ++i) {
            uint32_t b = (g >> (9 * (i % 7))) & 511;
            if (i % 7 == 6) g = mix(g);
            if (!f(blk + (b >> 6), b & 63)) return;
        }
    }
    void add(uint64_t kmer) { for_each_bit(kmer, [&](uint64_t w, uint32_t b) { bits[w] |= 1ULL << b; return true; }); }
    bool has(uint64_t kmer) const {
        bool ok = true;
        for_each_bit(kmer, [&](uint64_t w, uint32_t b) { return ok = (bits[w] >> b) & 1; });
        return ok;
    }
    void prefetch(uint64_t kmer) const { __builtin_prefetch(&bits[block_of(mix(kmer)) * 8]); }

    // calls f(i, canonical code) for every k-mer s[i..i+k) of A/C/G/T only
    template <class F> void for_each_kmer(const unsigned char *s, uint64_t n, F f) const {
        uint64_t fw = 0, rc = 0, len = 0, mask = k == 32 ? ~0ULL : (1ULL << (2 * k)) - 1;
        for (uint64_t i = 0; i < n; ++i) {
            int c = code(s[i]);
            if (c < 0) { len = 0; fw = rc = 0; continue; }
            fw = ((fw << 2) | (uint64_t)c) & mask;
            rc = (rc >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
            if (++len >= k) f(i + 1 - k, fw < rc ? fw : rc);
        }
    }

    void init(uint64_t nkeys, double bits_per_key, uint32_t hh) {
        h = hh;
        nblocks = std::max<uint64_t>(1, (uint64_t)std::ceil(nkeys * bits_per_key / 512.0));
        bits.assign(nblocks * 8, 0);
    }
    uint64_t bytes() const { return bits.size() * 8; }
    bool save(const std::string &f) const {
        FILE *o = fopen(f.c_str(), "wb");
        if (!o) return false;
        fwrite("RZKBLM01", 1, 8, o); fwrite(&k, 4, 1, o); fwrite(&h, 4, 1, o); fwrite(&nblocks, 8, 1, o);
        fwrite(bits.data(), 8, bits.size(), o);
        return fclose(o) == 0;
    }
    bool load(const std::string &f) {
        FILE *i = fopen(f.c_str(), "rb"); char mg[8];
        if (!i) return false;
        bool ok = fread(mg, 1, 8, i) == 8 && !memcmp(mg, "RZKBLM01", 8) && fread(&k, 4, 1, i) == 1 && fread(&h, 4, 1, i) == 1
                  && fread(&nblocks, 8, 1, i) == 1;
        if (ok) { bits.resize(nblocks * 8); ok = fread(bits.data(), 8, bits.size(), i) == bits.size(); }
        fclose(i);
        return ok;
    }
};

}  // namespace katka
