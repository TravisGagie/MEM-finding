// Blocked binary interpolative coding of the RLCSA's two sorted lists (no coder: only the code lengths are summed).
//
//   blocked_interp PREFIX B...      (PREFIX.csa, an explicit RLCSA; B the block sizes, e.g. 64 128)
//
// Each character's run starts (universe [0, n)) and their F positions (universe [C[c], C[c+1])) are cut into blocks of B
// values. Each block's first value (its head) goes in an Elias-Fano list over the list's universe. The block's other B - 1
// values are coded by BIC between its head and the next head (or the universe's end), with minimal binary or centred
// minimal binary codes, in preorder (middle value, left half, right half). Reported, in bits per value:
//   heads     the heads' Elias-Fano (l = floor(log2(U / h)) low bits per head, and h + (U >> l) high bits)
//   BIC, cen  the codes inside the blocks
//   ptr       an Elias-Fano list of each block's starting bit offset (needed to jump to a block; not part of "heads")
// and the expected number of values decoded per lookup, in preorder, a left half being skipped only by decoding it:
//   access    access(j), j uniform over the list's values (a head costs nothing)
//   search    count_less(x) inside the block that the heads select, x uniform over the universe (so a gap is chosen
//             with probability proportional to its width): the values decoded until x's gap is known
#include "../src/rlcsa.hpp"
#include <cmath>
using namespace katka;
static inline uint64_t lg(uint64_t x) { return 63 - __builtin_clzll(x); }
static inline double ef_bits(uint64_t m, uint64_t U) {
    if (!m) return 0;
    const uint64_t l = U / m > 1 ? lg(U / m) : 0;
    return (double)m * l + m + (double)(U >> l);
}
struct Tot { double heads = 0, bic = 0, cen = 0, ptr = 0, acc = 0, srch = 0, width = 0; uint64_t m = 0; };
// BIC code lengths of x[i..j] in (lo, hi)
static void code(const std::vector<uint64_t> &x, int64_t i, int64_t j, int64_t lo, int64_t hi, double &bic, double &cen) {
    while (i <= j) {
        const int64_t mid = (i + j) / 2;
        const int64_t a = lo + 1 + (mid - i), b = hi - 1 - (j - mid);
        const uint64_t r = (uint64_t)(b - a + 1), v = (uint64_t)((int64_t)x[mid] - a);
        if (r > 1) {
            const uint64_t k = lg(r), u = (2ull << k) - r;
            bic += v < u ? k : k + 1;
            const uint64_t s = (r - u) / 2;
            cen += (v >= s && v < s + u) ? k : k + 1;
        }
        code(x, i, mid - 1, lo, (int64_t)x[mid], bic, cen);
        i = mid + 1; lo = (int64_t)x[mid];
    }
}
// access: for each position in [i, j], the number of values decoded before it is known (its preorder rank + 1); summed
static double acc_cost(int64_t i, int64_t j, double before) {
    if (i > j) return 0;
    const int64_t mid = (i + j) / 2, left = mid - i;
    return (before + 1) + acc_cost(i, mid - 1, before + 1) + acc_cost(mid + 1, j, before + 1 + left);
}
// search: x uniform in (lo, hi) among the block's values x[i..j]; returns the sum over x of the values decoded
static double srch_cost(const std::vector<uint64_t> &x, int64_t i, int64_t j, uint64_t lo, uint64_t hi, double before) {
    if (i > j) return (double)(hi - lo) * before;          // x's gap known: every x in [lo, hi) cost `before`
    const int64_t mid = (i + j) / 2, left = mid - i;
    return srch_cost(x, i, mid - 1, lo, x[mid], before + 1) + srch_cost(x, mid + 1, j, x[mid], hi, before + 1 + left);
}
static void list(const std::vector<uint64_t> &x, uint64_t ulo, uint64_t uhi, uint64_t B, Tot &t) {
    const uint64_t m = x.size();
    if (!m) return;
    const uint64_t nb = (m + B - 1) / B;
    t.heads += ef_bits(nb, uhi - ulo); t.m += m;
    std::vector<uint64_t> off; uint64_t tot = 0;
    for (uint64_t s = 0; s < m; s += B) {
        const uint64_t e = std::min(m, s + B), hi = e < m ? x[e] : uhi;
        double b = 0, c = 0;
        code(x, (int64_t)s + 1, (int64_t)e - 1, (int64_t)x[s], (int64_t)hi, b, c);
        t.bic += b; t.cen += c; off.push_back(tot); tot += (uint64_t)b;
        t.acc += acc_cost((int64_t)s + 1, (int64_t)e - 1, 0);
        t.srch += srch_cost(x, (int64_t)s + 1, (int64_t)e - 1, x[s], hi, 0); t.width += (double)(hi - x[s]);
    }
    t.ptr += ef_bits(nb, tot + 1);
}
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: blocked_interp PREFIX B...\n"); return 1; }
    rlcsa X; X.load(std::string(argv[1]) + ".csa");
    if (X.ef || X.wide) { fprintf(stderr, "needs an explicit 32-bit RLCSA\n"); return 1; }
    const auto &R = X.rr;
    printf("%-5s %-12s %10s %7s %7s %7s %7s %8s %8s %7s %8s\n", "B", "list", "values", "heads", "BIC", "centred", "ptr",
           "h+BIC", "h+BIC+p", "access", "search");
    for (int a = 2; a < argc; ++a) {
        const uint64_t B = strtoull(argv[a], nullptr, 10);
        Tot S, F;
        for (int c = 0; c < 256; ++c) {
            const uint64_t nr = R.off[c + 1] - R.off[c] - 1;
            if (!nr) continue;
            std::vector<uint64_t> st(nr), fp(nr);
            for (uint64_t k = 0; k < nr; ++k) { st[k] = R.runs[R.off[c] + k].start; fp[k] = X.C[c] + R.runs[R.off[c] + k].cum; }
            list(st, 0, X.n, B, S); list(fp, X.C[c], X.C[c + 1], B, F);
        }
        for (auto [name, t] : {std::pair<const char *, Tot &>{"run starts", S}, {"F positions", F}}) {
            const double m = t.m;
            printf("%-5lu %-12s %10lu %7.3f %7.3f %7.3f %7.3f %8.3f %8.3f %7.2f %8.2f   (MB h+BIC %.1f, h+BIC+p %.1f)\n", (unsigned long)B, name,
                   (unsigned long)t.m, t.heads / m, t.bic / m, t.cen / m, t.ptr / m, (t.heads + t.bic) / m, (t.heads + t.bic + t.ptr) / m,
                   t.acc / m, t.srch / t.width, (t.heads + t.bic) / 8e6, (t.heads + t.bic + t.ptr) / 8e6);
        }
    }
    return 0;
}
