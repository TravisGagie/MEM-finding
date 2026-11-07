// Bits per value for the RLCSA's sorted lists under several codes (no coder: only the code lengths are summed).
//
//   interp_bits PREFIX            (PREFIX.csa, an explicit RLCSA)
//
// Lists, per character c: its run starts in the BWT (universe [0, n)), and their positions in F (universe [C[c], C[c+1])), both
// strictly increasing.  Codes:
//   EF          Elias-Fano: l = floor(log2(U / m)) low bits per value, and m + (U >> l) high bits.
//   BIC         binary interpolative coding (Moffat and Stuiver): the middle value of each node in the range its neighbours and
//               positions allow, by a minimal binary code over that range.
//   BIC-centred the same with centred minimal binary codes (the short codewords for the middle of the range).
//   mid-gamma   the middle value's signed deviation from the linear interpolation of the node's two known neighbours (their
//               average, when the middle is the centre), zig-zagged and coded with Elias gamma.
//   mid-Rice    the same deviation with a Rice code whose parameter is 2^round(log2(R / (2 sqrt(2k + 3)))), R the neighbours'
//               difference and 2k + 1 the node's values.
// Both midpoint codes recurse over the same nodes as BIC.
#include "../src/rlcsa.hpp"
#include <cmath>
using namespace katka;
struct Tot { double ef = 0, bic = 0, cen = 0, gam = 0, rice = 0; uint64_t m = 0; };
static inline uint64_t lg(uint64_t x) { return 63 - __builtin_clzll(x); }
static void code(const std::vector<uint64_t> &x, int64_t i, int64_t j, int64_t lo, int64_t hi, Tot &t) {
    // x[i..j] lie in (lo, hi), lo and hi being the known neighbours (or the universe's ends, one past)
    while (i <= j) {
        const int64_t mid = (i + j) / 2, cnt = j - i + 1;
        const int64_t a = lo + 1 + (mid - i), b = hi - 1 - (j - mid);   // the range of x[mid]
        const uint64_t r = (uint64_t)(b - a + 1), v = (uint64_t)((int64_t)x[mid] - a);
        if (r > 1) {
            const uint64_t k = lg(r), u = (2ull << k) - r;           // u codewords of k bits, r - u of k + 1
            t.bic += v < u ? k : k + 1;
            const uint64_t s = (r - u) / 2;
            t.cen += (v >= s && v < s + u) ? k : k + 1;
        }
        const double pred = lo + (double)(hi - lo) * (double)(mid - i + 1) / (double)(cnt + 1);
        const int64_t d = (int64_t)x[mid] - (int64_t)llround(pred);
        const uint64_t z = d >= 0 ? 2 * (uint64_t)d : 2 * (uint64_t)(-d) - 1;
        t.gam += 2 * lg(z + 1) + 1;
        const double bp = (double)(hi - lo) / (2.0 * std::sqrt((double)cnt + 2.0));
        const int q = bp > 1 ? (int)std::lround(std::log2(bp)) : 0;
        t.rice += (double)(z >> q) + 1 + q;
        code(x, i, mid - 1, lo, (int64_t)x[mid], t);
        i = mid + 1; lo = (int64_t)x[mid];
    }
}
static void list(const std::vector<uint64_t> &x, uint64_t ulo, uint64_t uhi, Tot &t) {
    const uint64_t m = x.size(), U = uhi - ulo;
    if (!m) return;
    const uint64_t l = U / m > 1 ? lg(U / m) : 0;
    t.ef += (double)m * l + m + (double)(U >> l); t.m += m;
    code(x, 0, (int64_t)m - 1, (int64_t)ulo - 1, (int64_t)uhi, t);
}
int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: interp_bits PREFIX\n"); return 1; }
    rlcsa X; X.load(std::string(argv[1]) + ".csa");
    if (X.ef || X.wide) { fprintf(stderr, "needs an explicit 32-bit RLCSA\n"); return 1; }
    const auto &R = X.rr;
    Tot S, F;
    printf("%-6s %-10s %10s %8s %8s %8s %8s %8s\n", "char", "list", "values", "EF", "BIC", "centred", "mid-g", "mid-Rice");
    for (int c = 0; c < 256; ++c) {
        const uint64_t nr = R.off[c + 1] - R.off[c] - 1;
        if (!nr) continue;
        std::vector<uint64_t> st(nr), fp(nr);
        for (uint64_t k = 0; k < nr; ++k) { st[k] = R.runs[R.off[c] + k].start; fp[k] = X.C[c] + R.runs[R.off[c] + k].cum; }
        Tot s, f;
        list(st, 0, X.n, s); list(fp, X.C[c], X.C[c + 1], f);
        for (auto [name, t] : {std::pair<const char *, Tot &>{"starts", s}, {"F", f}})
            printf("%-6s %-10s %10lu %8.3f %8.3f %8.3f %8.3f %8.3f\n", c >= 33 && c < 127 ? std::string(1, (char)c).c_str() : std::to_string(c).c_str(), name, (unsigned long)t.m,
                   t.ef / t.m, t.bic / t.m, t.cen / t.m, t.gam / t.m, t.rice / t.m);
        S.ef += s.ef; S.bic += s.bic; S.cen += s.cen; S.gam += s.gam; S.rice += s.rice; S.m += s.m;
        F.ef += f.ef; F.bic += f.bic; F.cen += f.cen; F.gam += f.gam; F.rice += f.rice; F.m += f.m;
    }
    for (auto [name, t] : {std::pair<const char *, Tot &>{"run starts", S}, {"F positions", F}})
        printf("all    %-10s %10lu %8.3f %8.3f %8.3f %8.3f %8.3f   (bits per value; MB: %.1f %.1f %.1f %.1f %.1f)\n", name, (unsigned long)t.m, t.ef / t.m, t.bic / t.m, t.cen / t.m, t.gam / t.m, t.rice / t.m,
               t.ef / 8e6, t.bic / 8e6, t.cen / 8e6, t.gam / 8e6, t.rice / 8e6);
    return 0;
}
