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
    if (argc != 3) { fprintf(stderr, "usage: rlcsa_profile PREFIX STEPS\n"); return 1; }
    const std::string P = argv[1];
    rlcsa X, Y; X.load(P + ".csa"); Y.load(P + ".efcsa");
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
    printf("Elias-Fano RLCSA (%.1f MB)\n", Y.bytes() / 1e6);
    double t_all = rep("whole step (rank2)", [&](const st &x) { uint64_t a, b; E.rank2(x.c, x.sp, x.ep, a, b); return a + b; });
    double t_b = rep("first run of c: KS.count_less(c n)", [&](const st &x) { return E.KS.count_less((uint64_t)x.c * E.n); });
    double t_s = rep("runs of sp and ep: 2 KS.count_less (predecessor search)", [&](const st &x) { uint64_t base = (uint64_t)x.c * E.n; return E.KS.count_less(base + x.sp) + E.KS.count_less(base + x.ep); });
    std::vector<uint64_t> A1(N), A2(N), BB(N);
    for (uint64_t i = 0; i < N; ++i) { uint64_t base = (uint64_t)S[i].c * E.n; BB[i] = E.KS.count_less(base); A1[i] = S[i].sp ? E.KS.count_less(base + S[i].sp) : BB[i]; A2[i] = S[i].ep ? E.KS.count_less(base + S[i].ep) : BB[i]; }
    uint64_t k = 0;
    double t_r = rep("from the runs: 2 rank_from (KS.access, FP.access2, FP.access(b))", [&](const st &x) { uint64_t i = k++ % N; uint64_t r1 = A1[i] == BB[i] ? 0 : E.rank_from(BB[i], A1[i], x.c, x.sp), r2 = A2[i] == BB[i] ? 0 : E.rank_from(BB[i], A2[i], x.c, x.ep); return r1 + r2; });
    printf("  sum of the parts %.1f ns against %.1f for the whole step\n", t_b + t_s + t_r, t_all);
    // the predecessor search's own work: elements compared within the high bucket (get_low), per count_less
    uint64_t lows = 0, buck = 0;
    for (uint64_t i = 0; i < N; ++i) for (uint64_t x : {S[i].sp, S[i].ep}) {
        uint64_t v = (uint64_t)S[i].c * E.n + x, h = v >> E.KS.l;
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
    printf("summary: Elias-Fano %.1f ns per step (first run of c %.0f%%, predecessor search %.0f%%, from the runs %.0f%%); explicit %.1f ns (search %.0f%%)\n",
           t_all, 100 * t_b / (t_b + t_s + t_r), 100 * t_s / (t_b + t_s + t_r), 100 * t_r / (t_b + t_s + t_r), x_all, 100 * x_f / x_all);
    return (int)(sink & 0);
}
