// Where the RLCSAs' backward steps spend their time: replays recorded steps (find built with -DKATKA_RECORD_STEPS, KATKA_RECORD=FILE)
// on the explicit and the Elias-Fano RLCSA, timing each part of a step separately.
//
//   rlcsa_profile PREFIX STEPS
//
// Elias-Fano (EFRank::rank2): b = KS.count_less(c n), the first run of c (the same for every step by c); a1, a2 = KS.count_less(c n + sp),
// KS.count_less(c n + ep), the runs holding sp and ep (a predecessor search over c's run starts); then rank_from for each, with one
// KS.access and an FP.access2, and FP.access(b).  Explicit (run_lists::rank): the jump table and a binary search over c's runs, then
// the rank from that run.
#include "../src/search.hpp"
#include <chrono>
using namespace katka;
static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: rlcsa_profile PREFIX STEPS [EFCSA]\n"); return 1; }
    const std::string P = argv[1];
    rlcsa X, Y; X.load(P + ".csa"); Y.load(argc > 3 ? argv[3] : (P + ".efcsa").c_str());
    struct st { unsigned char c; uint64_t sp, ep; };
    std::vector<st> S;
    { FILE *f = fopen(argv[2], "rb"); st x; while (fread(&x.c, 1, 1, f) == 1 && fread(&x.sp, 8, 1, f) == 1 && fread(&x.ep, 8, 1, f) == 1) if (Y.C[x.c + 1] > Y.C[x.c]) S.push_back(x); fclose(f); }
    const uint64_t N = S.size();
    printf("%lu steps (characters present)\n", N);
    volatile uint64_t sink = 0;
    const EFRank &E = Y.er;
    auto rep = [&](const char *name, auto f, double base_ns = 0) {
        double best = 1e30;
        for (int r = 0; r < 3; ++r) { double t = now_s(); uint64_t acc = 0; for (uint64_t i = 0; i < N; ++i) acc += f(S[i]); sink += acc; best = std::min(best, now_s() - t); }
        printf("  %-62s %7.1f ns per step%s\n", name, best / N * 1e9, base_ns > 0 ? "" : "");
        return best / N * 1e9;
    };
    {   // every step's answers against the explicit RLCSA's, and how often the two positions' runs, or their keys' buckets, coincide
        uint64_t bad = 0, same = 0, next = 0, hsame = 0, hnext = 0, both = 0;
        for (uint64_t i = 0; i < N; ++i) {
            const st &x = S[i]; uint64_t a, b; Y.er.rank2(x.c, x.sp, x.ep, a, b);
            if (a != X.rr.rank(x.c, x.sp) || b != X.rr.rank(x.c, x.ep)) ++bad;
            const EFRank &F = Y.er; uint64_t base = F.base_of(x.c), bb = F.first[x.c];
            uint64_t a1 = x.sp ? F.KS.count_less(base + x.sp) : bb, a2 = x.ep ? F.KS.count_less(base + x.ep) : bb;
            if (a1 == bb) continue;
            ++both; same += a1 == a2; next += a2 == a1 + 1;
            uint64_t h1 = (base + x.sp) >> F.KS.l, h2 = (base + x.ep) >> F.KS.l; hsame += h1 == h2; hnext += h2 == h1 + 1;
        }
        printf("check: %lu steps whose ranks differ from the explicit RLCSA's\n", (unsigned long)bad);
        printf("pairs (both positions past c's first run, %lu steps): same run %.1f%%, next run %.1f%%; keys in the same bucket %.1f%%, the next %.1f%%\n",
               (unsigned long)both, 100.0 * same / both, 100.0 * next / both, 100.0 * hsame / both, 100.0 * hnext / both);
        printf("fine select samples: every %lu ones, KS %.2f MB, FP %.2f MB\n", (unsigned long)EliasFano::SAMPLE1, Y.er.KS.f1.size() * 4 / 1e6, Y.er.FP.f1.size() * 4 / 1e6);
    }
    printf("Elias-Fano RLCSA (%.1f MB)\n", Y.bytes() / 1e6);
    double t_all = rep("whole step (rank2)", [&](const st &x) { uint64_t a, b; E.rank2(x.c, x.sp, x.ep, a, b); return a + b; });
    const bool cached = !E.first.empty();
    double t_b = rep(cached ? "first run of c and its F position: stored per character" : "first run of c: KS.count_less(c n), FP.access(b)",
                     [&](const st &x) { return cached ? E.first[x.c] + E.ffirst[x.c] : E.KS.count_less(E.base_of(x.c)) + E.FP.access(E.KS.count_less(E.base_of(x.c))); });
    double t_s = rep("runs of sp and ep: 2 KS.count_less (predecessor search)", [&](const st &x) { uint64_t base = E.base_of(x.c); return E.KS.count_less(base + x.sp) + E.KS.count_less(base + x.ep); });
    std::vector<uint64_t> A1(N), A2(N), BB(N), FB(N);
    for (uint64_t i = 0; i < N; ++i) { uint64_t base = E.base_of(S[i].c); BB[i] = E.KS.count_less(base); FB[i] = BB[i] < E.nruns ? E.FP.access(BB[i]) : E.n; A1[i] = S[i].sp ? E.KS.count_less(base + S[i].sp) : BB[i]; A2[i] = S[i].ep ? E.KS.count_less(base + S[i].ep) : BB[i]; }
    uint64_t k = 0;
    double t_r = rep("from the runs: 2 rank_from (KS.access, FP.access2)", [&](const st &x) { uint64_t i = k++ % N, base = E.base_of(x.c); uint64_t r1 = A1[i] == BB[i] ? 0 : E.rank_from(FB[i], A1[i], base, x.sp), r2 = A2[i] == BB[i] ? 0 : E.rank_from(FB[i], A2[i], base, x.ep); return r1 + r2; });
    printf("  sum of the parts %.1f ns against %.1f for the whole step\n", t_b + t_s + t_r, t_all);
    // the predecessor search's own work: elements compared within the high bucket (get_low), per count_less
    uint64_t lows = 0, buck = 0;
    for (uint64_t i = 0; i < N; ++i) for (uint64_t x : {S[i].sp, S[i].ep}) {
        uint64_t v = E.base_of(S[i].c) + x, h = v >> E.KS.l;
        uint64_t lo, hi;
        if (h == 0) { lo = 0; hi = E.KS.template select_<false>(0); }
        else { uint64_t p = E.KS.template select_<false>(h - 1); lo = p - (h - 1); uint64_t q = p + 1, w = q >> 6, word = ~E.KS.high[w] & (~0ull << (q & 63)); while (!word) word = ~E.KS.high[++w]; hi = ((w << 6) + __builtin_ctzll(word)) - h; }
        buck += hi - lo; uint64_t c = 0; while (hi - lo > 8) { ++c; uint64_t m = (lo + hi) >> 1; if (E.KS.get_low(m) < (v & ((1ull << E.KS.l) - 1))) lo = m + 1; else hi = m; } lows += c + std::min<uint64_t>(hi - lo, 8);
    }
    printf("  per search: %.2f entries in the high bucket, %.2f low-bit reads; KS: l = %u, %lu runs over universe %lu\n", (double)buck / (2 * N), (double)lows / (2 * N), E.KS.l, (unsigned long)E.KS.n, (unsigned long)E.KS.u);
    printf("Explicit RLCSA (%.1f MB)\n", X.bytes() / 1e6);
    double x_all = rep("whole step (2 run_lists::rank)", [&](const st &x) { return X.rr.rank(x.c, x.sp) + X.rr.rank(x.c, x.ep); });
    double x_f = rep("runs of sp and ep: 2 first_run (jump table, binary search)", [&](const st &x) { return X.rr.first_run(x.c, x.sp) + X.rr.first_run(x.c, x.ep); });
    printf("  rest (rank from the run): %.1f ns\n", x_all - x_f);
    printf("KS %.1f MB, FP %.1f MB, per-character values %.3f MB\n", E.KS.bytes() / 1e6, E.FP.bytes() / 1e6, (E.slot.size() * 4 + (E.first.size() + E.ffirst.size()) * 8) / 1e6);
    printf("summary: Elias-Fano %.1f ns per step (first run of c %.0f%%, predecessor search %.0f%%, from the runs %.0f%%); explicit %.1f ns (search %.0f%%)\n",
           t_all, 100 * t_b / (t_b + t_s + t_r), 100 * t_s / (t_b + t_s + t_r), 100 * t_r / (t_b + t_s + t_r), x_all, 100 * x_f / x_all);
    return (int)(sink & 0);
}
