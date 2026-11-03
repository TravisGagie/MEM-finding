// RLZ of the chr19 two-level parse against the CHM13 parse (an analysis program, not part of the tool).
//   rlz_parse PREFIX.ptxt REF_LEN [DUMP]   (DUMP: also write "start id" per phrase, for brute_force.py; REF_LEN: length of the first sequence, CHM13, column 4 of PREFIX.tbl)
// The parse is that of `katka index --parse-oracle` (syncmer parse of S = G1 X rc(G1) X G2 X rc(G2) X ...,
// lexicographic phrase ids 1..sigma; consecutive phrases overlap by K characters). The reference parse is the
// phrases starting in S[0 .. 2(REF_LEN+1)), i.e. CHM13 and its reverse complement; the target is the rest.
// Reports:
//   - H0 and H1 of the whole parse (bits per phrase);
//   - greedy RLZ (Kuruppu et al.): longest match in the reference, or a one-phrase literal if the phrase does
//     not occur there; factor count and mean length in phrases;
//   - GDC-style RLZ (Deorowicz & Grabowski): each factor is a longest match followed by one mismatch phrase;
//     with Ferrada et al.'s relative pointers (source - target position), choosing among equally long matches
//     the one continuing the previous pointer, and the number of runs of equal relative pointers;
//   - RLZAP-style (Cox et al.): after each match, look ahead j = 0..J phrases for a match of at least M phrases
//     at the relative offset rel + delta, |delta| <= D (smallest j, then smallest |delta|); if found, the j
//     skipped phrases are literals and delta a pointer change; otherwise a greedy longest match with an explicit
//     pointer. Run for (J, D, M) = (4, 4, 2) and (8, 8, 4);
//   - size estimates of an RLZ parse oracle (see the README for the formulas).
// Build (from tool/): g++ -std=c++17 -O3 -march=native -Isrc -Ithird_party/libsais/include
//   results/chr19_rlz/rlz_parse.cpp build/libsais.o build/libsais64.o -o bin/rlz_parse -lz -pthread
#include "tag_build.hpp"     // lce_bytes, which parse_text.hpp uses
#include "parse_text.hpp"
#include "libsais.h"
#include <cmath>
using namespace katka;
static double lg(double x) { return std::log2(x); }
static u64 cl2(u64 x) { u64 b = 0; while ((1ULL << b) < x) ++b; return b; }
static double ef_bits(double cnt, double uni) { return cnt * (2 + std::ceil(std::max(0.0, lg(uni / cnt)))); }
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: rlz_parse PREFIX.ptxt REF_LEN\n"); return 1; }
    double t0 = now();
    parse_text P; P.load(argv[1]);
    const u64 m = P.m, sigma = P.sigma, n = P.n, refend = 2 * (std::stoull(argv[2]) + 1);
    std::vector<int32_t> id(m);
    for (u64 i = 0; i < m; ++i) id[i] = (int32_t)P.ids.get(i);
    u64 R = 0; while (R < m && P.starts.access(R) < refend) ++R;
    if (argc > 3) { FILE *f = open_or_die(argv[3], "w"); fprintf(f, "%lu %lu %lu\n", n, R, sigma); for (u64 i = 0; i < m; ++i) fprintf(f, "%lu %d\n", P.starts.access(i), id[i]); fclose(f); }
    printf("parse: n = %lu characters, m = %lu phrases (%.1f characters per phrase step), sigma = %lu distinct, K = %u\n", n, m, (double)n / m, sigma, P.K);
    printf("reference parse (CHM13, both strands, S[0..%lu)): %lu phrases; target: %lu phrases\n", refend, R, m - R);
    // H0, H1
    { std::vector<u64> c(sigma + 1, 0); for (u64 i = 0; i < m; ++i) ++c[id[i]];
      double h = 0; for (u64 a = 1; a <= sigma; ++a) if (c[a]) h -= c[a] * lg((double)c[a] / m);
      std::vector<u64> pr(m - 1); for (u64 i = 1; i < m; ++i) pr[i - 1] = ((u64)id[i - 1] << 32) | (u64)id[i];
      std::sort(pr.begin(), pr.end());
      double h1 = 0; u64 i = 0, ctx = 0, distinct = 0;
      while (i < pr.size()) {           // per context: sum over next symbols of -n_cs log(n_cs / n_c)
          u64 j = i, nc = 0; while (j < pr.size() && (pr[j] >> 32) == (pr[i] >> 32)) ++j;
          nc = j - i; ++ctx;
          for (u64 a = i; a < j;) { u64 b = a; while (b < j && pr[b] == pr[a]) ++b; h1 -= (b - a) * lg((double)(b - a) / nc); ++distinct; a = b; }
          i = j;
      }
      printf("H0 = %.4f bits per phrase; H1 = %.4f bits per phrase (%lu contexts, %lu distinct pairs); log2(sigma) = %.2f\n",
             h / m, h1 / (m - 1), ctx, distinct, lg(sigma)); }
    // suffix array and inverse of the reference parse
    std::vector<int32_t> T(id.begin(), id.begin() + R), SA(R), ISA(R);
    if (libsais_int(T.data(), SA.data(), (int32_t)R, (int32_t)sigma + 1, 0) != 0) die("libsais_int failed");
    for (u64 x = 0; x < R; ++x) ISA[SA[x]] = (int32_t)x;
    auto at = [&](int64_t q) -> int64_t { return q < (int64_t)R ? T[q] : -1; };
    // the longest match of id[p..) in T: returns its length and leaves the SA interval in [lo, hi)
    auto match = [&](u64 p, u64 &lo, u64 &hi) {
        lo = 0; hi = R; u64 d = 0;
        while (p + d < m) {
            int64_t c = id[p + d];
            u64 a = lo, b = hi;                                    // first x in [lo,hi) with T[SA[x]+d] >= c
            while (a < b) { u64 mid = (a + b) / 2; if (at((int64_t)SA[mid] + d) < c) a = mid + 1; else b = mid; }
            u64 l2 = a; b = hi;                                    // first x with T[SA[x]+d] > c
            while (a < b) { u64 mid = (a + b) / 2; if (at((int64_t)SA[mid] + d) <= c) a = mid + 1; else b = mid; }
            if (l2 == a) break;
            lo = l2; hi = a; ++d;
        }
        return d;
    };
    // plain RLZ
    u64 z = 0, lit = 0, covered = 0; std::vector<u64> fstart_chr; fstart_chr.reserve(1 << 20);
    for (u64 p = R; p < m;) {
        u64 lo, hi, d = match(p, lo, hi);
        fstart_chr.push_back(P.starts.access(p));
        ++z; if (d == 0) { ++lit; p += 1; } else { covered += d; p += d; }
    }
    double tz = (double)(m - R);
    printf("RLZ (greedy longest match): %lu factors (%lu literals), %.2f phrases per factor (%.2f per non-literal factor)\n",
           z, lit, tz / z, (double)covered / std::max<u64>(1, z - lit));
    // GDC with relative pointers
    u64 zg = 0, runs = 0; int64_t prev = INT64_MIN;
    for (u64 p = R; p < m;) {
        u64 lo, hi, d = match(p, lo, hi); int64_t q = -1;
        if (d > 0) {
            q = SA[lo];
            if (prev != INT64_MIN) { int64_t c = (int64_t)p + prev; if (c >= 0 && c < (int64_t)R && (u64)ISA[c] >= lo && (u64)ISA[c] < hi) q = c; }
            int64_t rel = q - (int64_t)p; if (rel != prev) { ++runs; prev = rel; }
        }
        ++zg; p += d + 1;                                          // the match, then one mismatch phrase
    }
    printf("GDC-style RLZ (match + one mismatch phrase): %lu factors, %.2f phrases per factor; %lu runs of equal relative pointers\n",
           zg, tz / zg, runs);
    // sizes (bits)
    const double w = cl2(sigma + 1), lgR = cl2(R);
    double starts = ef_bits(z, n), offs = (double)z * lgR, lits = (double)lit * w, refp = (double)R * w, refst = ef_bits(R, refend);
    double dict = 8.0 * (P.dict.size() + P.doff.bytes());
    double cur_ids = 8.0 * P.ids.bytes(), cur_starts = 8.0 * P.starts.bytes();
    auto MB = [](double bits) { return bits / 8e6; };
    printf("sizes (MB):\n");
    printf("  current parse oracle: ids %.1f + phrase starts %.1f + dictionary %.1f = %.1f\n", MB(cur_ids), MB(cur_starts), MB(dict), MB(cur_ids + cur_starts + dict));
    printf("  RLZ parse oracle: factor starts (EF over S positions) %.1f + offsets (%u bits each) %.1f + literals %.1f\n"
           "                   + reference parse (ids %.1f + phrase starts EF %.1f) + dictionary %.1f = %.1f\n",
           MB(starts), (unsigned)lgR, MB(offs), MB(lits), MB(refp), MB(refst), MB(dict), MB(starts + offs + lits + refp + refst + dict));
    double gstarts = ef_bits(zg, n), gmis = (double)zg * w, grel = (double)runs * (lgR + 1 + 2 + std::ceil(std::max(0.0, lg((double)zg / runs))));
    printf("  GDC + relative pointers (estimate): factor starts %.1f + mismatch phrases %.1f + pointer runs (value + EF of run starts) %.1f\n"
           "                   + reference parse %.1f + dictionary %.1f = %.1f\n",
           MB(gstarts), MB(gmis), MB(grel), MB(refp + refst), MB(dict), MB(gstarts + gmis + grel + refp + refst + dict));
    // RLZAP-style parse and its size
    auto gamma = [](u64 x) -> double { u64 b = 0; while ((x >> b) > 1) ++b; return 2.0 * b + 1; };    // x >= 1
    auto rlzap = [&](u64 J, int64_t D, u64 M) {
        u64 zf = 0, lits = 0, changes = 0, explicit_ = 0, litfactors = 0; double litcnt_bits = 0, ptr_bits = 0; u64 gap = 0;
        bool have = false; int64_t rel = 0;
        auto lce = [&](u64 p, int64_t q) { u64 l = 0; while (p + l < m && q + (int64_t)l < (int64_t)R && id[p + l] == T[q + l]) ++l; return l; };
        u64 p = R;
        while (p < m) {
            bool done = false;
            if (have) {
                for (u64 j = 0; j <= J && !done && p + j < m; ++j)
                    for (int64_t a = 0; a <= D && !done; ++a)
                        for (int sg = 0; sg < (a ? 2 : 1) && !done; ++sg) {
                            int64_t dl = sg ? -a : a, src = (int64_t)(p + j) + rel + dl;
                            if (src < 0 || src >= (int64_t)R) continue;
                            u64 l = lce(p + j, src);
                            if (l >= M) {
                                ++zf; lits += j; litcnt_bits += gamma(j + 1);
                                if (dl != 0) { ++changes; ptr_bits += gamma(gap + 1) + 1 + 1 + gamma((u64)a); gap = 0; rel += dl; } else ++gap;
                                p += j + l; done = true;
                            }
                        }
            }
            if (done) continue;
            u64 lo, hi, d = match(p, lo, hi);
            ++zf;
            if (d == 0) { ++lits; ++litfactors; litcnt_bits += gamma(2); ++gap; p += 1; continue; }   // one literal, no match
            int64_t q = SA[lo];
            if (have) { int64_t c = (int64_t)p + rel; if (c >= 0 && c < (int64_t)R && (u64)ISA[c] >= lo && (u64)ISA[c] < hi) q = c; }
            int64_t nr = q - (int64_t)p;
            litcnt_bits += gamma(1);
            if (!have || nr != rel) { ++explicit_; ptr_bits += gamma(gap + 1) + 1 + lgR; gap = 0; rel = nr; have = true; } else ++gap;
            p += d;
        }
        double st = ef_bits(zf, n), lb = (double)lits * w;
        printf("RLZAP-style (J = %lu, D = %ld, M = %lu): %lu factors, %.2f phrases per factor; %lu literal phrases (%lu factors that are only a literal),"
               " %lu pointer changes (|delta| <= %ld), %lu explicit pointers\n", J, D, M, zf, tz / zf, lits, litfactors, changes, D, explicit_);
        printf("  size (MB): factor starts %.1f + literal counts %.1f + literals (%u bits) %.1f + pointer differences (run-length) %.1f"
               " + reference parse %.1f + dictionary %.1f = %.1f\n", MB(st), MB(litcnt_bits), (unsigned)w, MB(lb), MB(ptr_bits),
               MB(refp + refst), MB(dict), MB(st + litcnt_bits + lb + ptr_bits + refp + refst + dict));
    };
    rlzap(4, 4, 2);
    rlzap(8, 8, 4);
    printf("(%.1f s)\n", now() - t0);
    return 0;
}
