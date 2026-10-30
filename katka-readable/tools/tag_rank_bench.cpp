// tag_rank_bench: counting with rank queries (Appendix C) against scanning the tag runs (Section 6).
//
// usage: tag_rank_bench index.tag intervals.bin [reps]
//   index.tag      an unsampled tag array (s = 1) with its .rmq (katka-tagbuild without -n)
//   intervals.bin  the BWT intervals [sp, ep) of MEMs, as pairs of 64-bit integers (rz-classify with
//                  RZ_IVLOG=intervals.bin, or any other source)
//
// For every interval, it counts the rows of each tag in two ways and checks that they agree:
//   scan  the tag array's count(): scan the tag runs overlapping [sp, ep), adding up clipped lengths;
//   rank  list the distinct tags with Sadakane's version of Muthukrishnan's algorithm (tag_array::list),
//         then count each tag t as rank_t(ep) - rank_t(sp), where rank_t is a predecessor search in the
//         sorted list of t's tag runs with their cumulative lengths (run_lists from rlcsa.hpp, with jump
//         tables), as for the RLCSA.
// It reports the time per interval of each, of the listing alone, the space of the rank lists, and the
// times split by the number of tag runs in the interval.
#include "../src/search.hpp"
#include "../src/tag_array.hpp"
#include <map>

using namespace katka;

template <class W> struct tag_rank {
    run_lists<W> rl;
    // per-tag lists of the tag runs (start, rows of the tag before it), plus jump tables, as in rlcsa::build_explicit
    void build(tag_array &T, uint32_t ntags) {
        typedef typename run_lists<W>::Run Run;
        rl.n = T.n; rl.sigma = ntags; rl.nruns = T.rho;
        std::vector<uint64_t> nr(ntags, 0);
        for (u64 i = 0; i < T.rho; ++i) nr[T.tag(i)]++;
        rl.off.assign(ntags + 1, 0);
        for (uint32_t c = 0; c < ntags; ++c) rl.off[c + 1] = rl.off[c] + nr[c] + 1;
        rl.runs.assign(rl.off[ntags], Run{0, 0});
        std::vector<uint64_t> next(rl.off.begin(), rl.off.end() - 1), seen(ntags, 0);
        for (u64 i = 0; i < T.rho; ++i) {
            u64 t = T.tag(i), s = T.run_start(i);
            rl.runs[next[t]++] = Run{(W)s, (W)seen[t]}; seen[t] += T.run_end(i) - s;
        }
        for (uint32_t c = 0; c < ntags; ++c) rl.runs[next[c]] = Run{(W)T.n, (W)seen[c]};
        rl.joff.assign(ntags + 1, 0); rl.shift.assign(ntags, 0); rl.jumps.clear();
        for (uint32_t c = 0; c < ntags; ++c) {
            rl.joff[c] = rl.jumps.size();
            uint64_t Rc = rl.off[c + 1] - rl.off[c] - 1;
            if (Rc < run_lists<W>::JUMP_MIN_RUNS) continue;
            uint8_t sh = 0;
            while ((T.n >> sh) > Rc / 8) ++sh;
            rl.shift[c] = sh;
            const Run *r = &rl.runs[rl.off[c]];
            for (uint64_t b = 0, k = 0; b < (T.n >> sh) + 2; ++b) {
                while (k < Rc && r[k].start < (b << sh)) ++k;
                rl.jumps.push_back((uint32_t)k);
            }
        }
        rl.joff[ntags] = rl.jumps.size();
    }
    // (tag, count) pairs sorted by tag, as tag_array::count gives them
    void count(tag_array &T, u64 sp, u64 ep, std::vector<std::pair<u64, u64>> &out, std::vector<u64> &tags) const {
        tags.clear(); T.list(sp, ep, tags); std::sort(tags.begin(), tags.end());
        out.clear();
        for (u64 t : tags) out.push_back({t, rl.rank((uint32_t)t, ep) - rl.rank((uint32_t)t, sp)});
    }
};

template <class W> static int run(tag_array &T, const std::vector<interval> &iv, int reps) {
    uint32_t ntags = 0;
    for (u64 i = 0; i < T.rho; ++i) ntags = std::max<uint32_t>(ntags, (uint32_t)T.tag(i) + 1);
    double t0 = now();
    tag_rank<W> R; R.build(T, ntags);
    fprintf(stderr, "rank lists built in %.1f s: %u tags, %lu runs, %.3f GB (%.1f bits/run); tag array %.3f GB\n",
            now() - t0, ntags, T.rho, R.rl.bytes() / 1e9, 8.0 * R.rl.bytes() / T.rho, T.bytes() / 1e9);

    // check: both methods give the same counts on every interval
    std::vector<std::pair<u64, u64>> a, b; std::vector<u64> tags;
    u64 nruns_tot = 0, ntags_tot = 0;
    std::vector<u64> nruns(iv.size());
    for (size_t k = 0; k < iv.size(); ++k) {
        T.count(iv[k].sp, iv[k].ep, a); R.count(T, iv[k].sp, iv[k].ep, b, tags);
        if (a != b) { fprintf(stderr, "MISMATCH at interval %zu [%lu, %lu)\n", k, iv[k].sp, iv[k].ep); return 1; }
        nruns[k] = T.run_of(iv[k].ep - 1) - T.run_of(iv[k].sp) + 1; nruns_tot += nruns[k]; ntags_tot += a.size();
    }
    printf("intervals %zu: identical counts; mean %.2f tag runs and %.2f tags per interval\n",
           iv.size(), (double)nruns_tot / iv.size(), (double)ntags_tot / iv.size());

    // timing, overall and by the number of tag runs in the interval
    const u64 edges[] = {1, 2, 4, 16, 64, 256, 1024, ~0ULL};
    const int NB = 7;
    double ts[NB + 1] = {0}, tr[NB + 1] = {0}, tl[NB + 1] = {0}; u64 cnt[NB + 1] = {0};
    u64 sink = 0;
    auto bucket = [&](u64 r) { int q = 0; while (q + 1 < NB && r >= edges[q + 1]) ++q; return q; };
    for (int rep = 0; rep < reps; ++rep)
        for (size_t k = 0; k < iv.size(); ++k) {
            int q = bucket(nruns[k]); if (rep == 0) ++cnt[q];
            double x = now(); T.count(iv[k].sp, iv[k].ep, a); sink += a.size();
            double y = now(); R.count(T, iv[k].sp, iv[k].ep, b, tags); sink += b.size();
            double z = now(); tags.clear(); T.list(iv[k].sp, iv[k].ep, tags); sink += tags.size();
            double w = now();
            ts[q] += y - x; tr[q] += z - y; tl[q] += w - z;
        }
    double S = 0, Rr = 0, Lr = 0;
    for (int q = 0; q < NB; ++q) { S += ts[q]; Rr += tr[q]; Lr += tl[q]; }
    double den = 1e-6 * reps;
    printf("per interval (us):  scan %.3f   rank %.3f   (listing alone %.3f)   [%d reps, sink %lu]\n",
           S / den / iv.size(), Rr / den / iv.size(), Lr / den / iv.size(), reps, sink % 10);
    printf("  tag runs     intervals    scan us    rank us   list us\n");
    for (int q = 0; q < NB; ++q) {
        if (!cnt[q]) continue;
        char lab[32];
        if (edges[q + 1] == ~0ULL) snprintf(lab, sizeof lab, ">= %lu", edges[q]);
        else if (edges[q + 1] == edges[q] + 1) snprintf(lab, sizeof lab, "%lu", edges[q]);
        else snprintf(lab, sizeof lab, "%lu-%lu", edges[q], edges[q + 1] - 1);
        printf("  %-10s %11lu %10.3f %10.3f %9.3f\n", lab, cnt[q], ts[q] / den / cnt[q], tr[q] / den / cnt[q], tl[q] / den / cnt[q]);
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) die("usage: tag_rank_bench index.tag intervals.bin [reps]");
    int reps = argc > 3 ? std::max(1, atoi(argv[3])) : 3;
    tag_array T; T.load(argv[1], true);
    if (T.s != 1) die("tag_rank_bench needs an unsampled tag array (s = 1)");
    std::vector<interval> iv;
    {
        FILE *f = fopen(argv[2], "rb"); if (!f) die(std::string("cannot open ") + argv[2]);
        u64 v[2];
        while (fread(v, 8, 2, f) == 2) if (v[0] < v[1] && v[1] <= T.n) iv.push_back({v[0], v[1]});
        fclose(f);
    }
    if (iv.empty()) die("no intervals");
    fprintf(stderr, "tag array: n=%lu, %lu runs; %zu intervals\n", T.n, T.rho, iv.size());
    return T.n >> 32 ? run<uint64_t>(T, iv, reps) : run<uint32_t>(T, iv, reps);
}
