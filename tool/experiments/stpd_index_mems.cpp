// stpd-index-mems: MEMs with the published STPD-index (github.com/regindex/STPD-index; Becker et al.,
// "Compressing Suffix Trees by Path Decompositions", ICALP 2026), for comparison with katka find.
//
// Their tools locate whole patterns; MEM-finding needs only the longest prefix of a string that occurs,
// which their Algorithm 1 computes on the way.  On their index of a text T holding both strands (katka
// stpd-prep), the longest suffix of X[lo..end) that occurs is the longest prefix of rc(X[lo..end)) that
// occurs, so BML and forward-backward search run on that one primitive, exactly as in katka (mems.hpp):
//   BML: window R[x..x+L) -> longest prefix of rc(window); then the left extension (longest prefix of
//   rc(R[0..x+L))); then forward from a (longest prefix of R[a..n)); report R[a..e); x = e - L + 1.
//
// The longest prefix of Y, with their components (text oracle O, sampled prefix array S): extend with
// O.LCP while the text agrees; at a mismatch after Y[0..i), find the sample that ends the colex-smallest
// text prefix suffixed by Y[0..i+1) (their tabulated binary search for i + 1 >= their table length, else a
// binary search over their samples with O.LCS_char), then check it with O.LCS: their binary search returns
// a candidate without checking it (they check once, at the end of a whole pattern), and a candidate that
// fails means Y[0..i+1) does not occur, so the longest prefix is i, exactly.
//
// usage: stpd-index-mems INDEX patterns.fa L [--fb] [--locate] [--pieces T.pieces] > mems.tsv
//   INDEX: built by stpd-build (with INDEX.param); patterns: FASTA, split at non-A/C/G/T characters;
//   --locate: also locate every occurrence of each MEM with their locate_all_occurrences (phi^-1 with an
//   exponential search and LCS checks), timed separately; --pieces: flag MEMs whose first occurrence (in
//   T) crosses a junction between A/C/G/T pieces (T has no separators, so such a match may be longer than
//   in katka's text, which has them).
// Output: name, start, end (0-based, half-open in the pattern), [occurrences], [junction flag].
#include <bits/stdc++.h>                   // everything else first: only their index and binary search are opened up
#include <malloc_count.h>                   // (in the order colex-index.hpp includes them)
#include <r-index_phi_inv.hpp>
#include <move-r_phi_inv.hpp>
#include <RLZ/RLZ_DNA_sux.hpp>
#include <bitpacked_text_oracle.hpp>
#define private public                      // (their index keeps its components private)
#include <tabulated_binary_search_DNA.hpp>
#include <colex-index.hpp>
#undef private
#include <chrono>
#include <fstream>
#include <iostream>

using u64 = uint64_t;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static char comp(char c) { switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return c; } }

struct counters {
    u64 searches = 0, jumps = 0, failed_checks = 0, chars = 0, checked_chars = 0, copied = 0;
    double t_lcp = 0, t_search = 0, t_check = 0;   // time in their LCP, their binary search, our check
};

template <class Index> struct engine {
    const Index &I;
    counters cnt;
    std::string Y;                                   // the string being matched (a copy, from its start)

    explicit engine(const Index &idx) : I(idx) {}

    // the sample for the colex-smallest prefix of T suffixed by Y[0..e), or -1 (exact, checked)
    long sample_for(u64 e) {
        const auto &S = I.S; const auto &O = I.O;
        long a;
        double t0 = now();
        if ((long)e >= (long)S.len) a = S.binary_search_with_tabluation(Y, 0, e);
        else {                                       // short: a binary search over the samples
            u64 lo = 0, hi = S.S;
            while (lo < hi) {
                u64 mid = (lo + hi) / 2, p = S.samples.get_value(mid) >> S.log_l;
                auto j = O.LCS_char(Y, e - 1, p);
                if (j.first >= e) { hi = mid; continue; }
                unsigned char tc = (unsigned char)j.second, pc = (unsigned char)Y[e - 1 - j.first];
                if (tc == 255 || tc < pc) lo = mid + 1; else hi = mid;
            }
            a = lo < S.S ? (long)(S.samples.get_value(lo) >> S.log_l) : -1;
        }
        double t1 = now(); cnt.t_search += t1 - t0;
        if (a < 0) return -1;
        u64 l = O.LCS(Y, e - 1, a);
        cnt.t_check += now() - t1; cnt.checked_chars += std::min<u64>(l + 1, e);
        if (l < e) { ++cnt.failed_checks; return -1; }
        return a;
    }
    // the length of the longest prefix of X[s..s+cap) that occurs in T (X a whole segment string)
    u64 longest_prefix(const std::string &X, u64 s, u64 cap) {
        ++cnt.searches;
        const auto &O = I.O;
        u64 have = std::min<u64>(cap, 4096);
        Y.assign(X, s, have);
        u64 i = 0, t = 0;                            // Y[0..i) occurs ending at text position t (exclusive)
        cnt.copied += have;
        for (;;) {
            double t0 = now();
            u64 f = O.LCP(Y, i, t);
            cnt.t_lcp += now() - t0;
            i += f; t += f; cnt.chars += f;
            if (i == have && have < cap) {           // the copy ran out: extend it and go on
                have = std::min<u64>(cap, 2 * have);
                Y.assign(X, s, have); cnt.copied += have;
                continue;
            }
            if (i >= have) return have;
            long a = sample_for(i + 1);
            ++cnt.jumps;
            if (a < 0) return i;
            t = (u64)a + 1; ++i;
        }
    }
};

struct mem { u64 a, e; };

template <class Index> static void bml(engine<Index> &E, const std::string &R, const std::string &RC, u64 L, std::vector<mem> &out) {
    const u64 n = R.size();
    u64 x = 0;
    while (x + L <= n) {
        u64 l = E.longest_prefix(RC, n - x - L, L);        // the longest suffix of R[x..x+L)
        if (l < L) { x = x + L - l; continue; }
        u64 a = x + L - E.longest_prefix(RC, n - x - L, x + L);
        u64 e = a + E.longest_prefix(R, a, n - a);
        out.push_back({a, e});
        x = e + 1 >= L ? e + 1 - L : 0;
    }
}
template <class Index> static void fb(engine<Index> &E, const std::string &R, const std::string &RC, u64 L, std::vector<mem> &out) {
    const u64 m = R.size();
    u64 s0 = 0;
    while (s0 < m) {
        u64 e = s0 + E.longest_prefix(R, s0, m - s0);
        if (e <= s0) { ++s0; continue; }
        if (e - s0 >= L) out.push_back({s0, e});
        if (e >= m) break;
        u64 s2 = e + 1 - E.longest_prefix(RC, m - e - 1, e + 1);
        s0 = s2 > e ? e + 1 : s2;
    }
}

template <class Index> static int run(const std::string &idxfile, const std::string &patfile, u64 L, bool use_fb, bool locate, const std::string &piecesfile) {
    Index I;
    {   // their load reports on standard output: send it to standard error, which keeps stdout for the MEMs
        std::streambuf *old = std::cout.rdbuf(std::cerr.rdbuf());
        I.load(idxfile, !locate && piecesfile.empty());
        std::cout.rdbuf(old);
    }
    std::vector<u64> pieces;
    if (!piecesfile.empty()) { std::ifstream in(piecesfile); u64 x; while (in >> x) pieces.push_back(x); }
    engine<Index> E(I);
    std::ifstream in(patfile);
    if (!in) { std::cerr << "cannot open " << patfile << "\n"; return 1; }
    std::string line, name, seq;
    std::vector<std::pair<std::string, std::string>> pats;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') { if (!name.empty()) pats.push_back({name, seq}); name = line.substr(1, line.find_first_of(" \t") - 1); seq.clear(); }
        else for (char c : line) seq.push_back((char)toupper((unsigned char)c));
    }
    if (!name.empty()) pats.push_back({name, seq});
    double tfind = 0, tloc = 0;
    u64 nm = 0, nocc = 0, njunction = 0;
    std::string out;
    for (auto &[nm_, P] : pats) {
        for (u64 i = 0; i < P.size();) {
            while (i < P.size() && !strchr("ACGT", P[i])) ++i;
            u64 r0 = i;
            while (i < P.size() && strchr("ACGT", P[i])) ++i;
            if (i - r0 < L) continue;
            std::string R = P.substr(r0, i - r0), RC(R.rbegin(), R.rend());
            for (auto &c : RC) c = comp(c);
            std::vector<mem> ms;
            double t0 = now();
            if (use_fb) fb(E, R, RC, L, ms); else bml(E, R, RC, L, ms);
            tfind += now() - t0;
            for (auto &M : ms) {
                out += nm_ + "\t" + std::to_string(r0 + M.a) + "\t" + std::to_string(r0 + M.e);
                if (locate) {
                    double t1 = now();
                    auto occ = I.locate_all_occurrences(R.substr(M.a, M.e - M.a));
                    tloc += now() - t1;
                    out += "\t" + std::to_string(occ.size()); nocc += occ.size();
                }
                if (!pieces.empty()) {               // does the first occurrence cross a junction?
                    long p = I.locate_primary_occurrence(R.substr(M.a, M.e - M.a));
                    bool cross = false;
                    if (p >= 0) {
                        u64 st = (u64)p, en = st + (M.e - M.a);
                        auto it = std::upper_bound(pieces.begin(), pieces.end(), st);
                        cross = it != pieces.end() && *it < en;
                    }
                    out += cross ? "\tjunction" : "\t-"; njunction += cross;
                }
                out += "\n"; ++nm;
            }
        }
    }
    std::cout << out;
    std::cerr << "stpd-index-mems: " << pats.size() << " patterns, " << nm << " MEMs (L = " << L << ", " << (use_fb ? "forward-backward" : "BML")
              << "); finding them " << tfind << " s (" << E.cnt.searches << " longest-prefix searches, " << E.cnt.jumps << " jumps, "
              << E.cnt.failed_checks << " candidates rejected by the check; time in their LCP " << E.cnt.t_lcp << " s, their binary search "
              << E.cnt.t_search << " s, the check " << E.cnt.t_check << " s over " << E.cnt.checked_chars << " characters; " << E.cnt.copied
              << " characters copied)";
    if (locate) std::cerr << "; locating " << nocc << " occurrences " << tloc << " s";
    if (!pieces.empty()) std::cerr << "; " << njunction << " MEMs whose first occurrence crosses a junction";
    std::cerr << "\n";
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) { std::cerr << "usage: stpd-index-mems INDEX patterns.fa L [--fb] [--locate] [--pieces T.pieces] > mems.tsv\n"; return 1; }
    std::string idx = argv[1], pat = argv[2], pieces;
    u64 L = std::stoull(argv[3]);
    bool use_fb = false, locate = false;
    for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--fb") use_fb = true; else if (a == "--locate") locate = true;
        else if (a == "--pieces" && i + 1 < argc) pieces = argv[++i];
        else { std::cerr << "unknown option " << a << "\n"; return 1; }
    }
    std::string phi, oracle;
    { std::ifstream in(idx + ".param"); in >> phi >> oracle; }
    using namespace stpd;
    if (oracle == "RLZ" && phi == "r-index")
        return run<colex_minus_index<tabulated_binary_search_DNA<RLZ_DNA_sux<>>, RLZ_DNA_sux<>, r_index_phi_inv>>(idx, pat, L, use_fb, locate, pieces);
    if (oracle == "bitpacked" && phi == "r-index")
        return run<colex_minus_index<tabulated_binary_search_DNA<bitpacked_text_oracle>, bitpacked_text_oracle, r_index_phi_inv>>(idx, pat, L, use_fb, locate, pieces);
    if (oracle == "RLZ" && phi == "move")
        return run<colex_minus_index<tabulated_binary_search_DNA<RLZ_DNA_sux<>>, RLZ_DNA_sux<>, move_r_phi_inv>>(idx, pat, L, use_fb, locate, pieces);
    std::cerr << "unsupported index variant " << oracle << " / " << phi << "\n";
    return 1;
}
