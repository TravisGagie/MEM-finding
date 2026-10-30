// katka-kbbuild: build the KeBaB k-mer filter (kmer_filter.hpp) of the text, for katka-classify -K.
// usage: katka-kbbuild [-k k] [-e eps] [-H h] prefix.S out.kbf
//   k = 20 (as in KeBaB), eps = 0.1 the target false-positive rate, h = the number of bits set per k-mer
//   (default: the optimal round(ln 2 * bits per k-mer)).
// Two passes over S: the first estimates the number of distinct canonical k-mers with a HyperLogLog sketch
// (2^16 registers), which sizes the filter; the second inserts them.
#include "common.hpp"
#include "kmer_filter.hpp"
#include <unistd.h>

using namespace katka;

int main(int argc, char **argv) {
    uint32_t k = 20, H = 0; double eps = 0.1; int o;
    while ((o = getopt(argc, argv, "k:e:H:")) != -1) {
        if (o == 'k') k = atoi(optarg); else if (o == 'e') eps = atof(optarg); else if (o == 'H') H = atoi(optarg); else return 1;
    }
    if (argc - optind != 2 || k < 1 || k > 32 || eps <= 0 || eps >= 1) die("usage: katka-kbbuild [-k k<=32] [-e eps] [-H h] prefix.S out.kbf");
    kmer_filter F; F.k = k;
    double t0 = now();
    // calls f(position, canonical k-mer) for every k-mer of S, reading it in blocks (carrying k-1 bytes over)
    auto stream = [&](auto f) {
        FILE *in = fopen(argv[optind], "rb");
        if (!in) die(std::string("cannot open ") + argv[optind]);
        const size_t B = 1 << 24;
        std::vector<unsigned char> buf(B + 32);
        size_t keep = 0, got;
        while ((got = fread(buf.data() + keep, 1, B, in)) > 0) {
            size_t len = keep + got;
            F.for_each_kmer(buf.data(), len, f);
            keep = std::min<size_t>(k - 1, len);
            memmove(buf.data(), buf.data() + len - keep, keep);
        }
        fclose(in);
    };
    // HyperLogLog: register j keeps the largest number of leading zeros (+1) among hashes in bucket j
    const int P = 16;
    std::vector<uint8_t> reg(1u << P, 0);
    uint64_t total = 0;
    stream([&](uint64_t, uint64_t c) {
        uint64_t x = kmer_filter::mix(c ^ 0x5bd1e995ULL);
        uint32_t j = x >> (64 - P);
        uint8_t r = (uint8_t)(__builtin_clzll((x << P) | (1ULL << (P - 1))) + 1);
        if (r > reg[j]) reg[j] = r;
        ++total;
    });
    double m = 1u << P, z = 0; uint64_t zeros = 0;
    for (auto r : reg) { z += std::ldexp(1.0, -r); if (!r) ++zeros; }
    double est = 0.7213 / (1 + 1.079 / m) * m * m / z;
    if (est < 2.5 * m && zeros) est = m * std::log(m / zeros);           // small-range correction
    double bpk = -std::log(eps) / (std::log(2.0) * std::log(2.0)) * 1.07;  // optimal bits per k-mer, +7% for blocking
    uint32_t h = H ? H : std::max(1u, (uint32_t)std::lround(std::log(2.0) * bpk));
    F.init((uint64_t)est, bpk, h);
    stream([&](uint64_t, uint64_t c) { F.add(c); });
    if (!F.save(argv[optind + 1])) die(std::string("cannot write ") + argv[optind + 1]);
    fprintf(stderr, "k=%u: %lu k-mers, ~%.0f distinct canonical (HLL); eps=%.3f -> %.2f bits/k-mer, h=%u; %lu blocks, %.3f GB | %.0f s\n",
            k, total, est, eps, bpk, h, F.nblocks, F.bytes() / 1e9, now() - t0);
    return 0;
}
