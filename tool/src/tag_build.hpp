// Building the tag array without a suffix array (Section 4 of the paper, "Construction"):
//   1. one pass of LF over each dataset, from its terminator to its first character, with "move"-style
//      arrays of O(r) words, gives the suffix-array samples at the ends of the BWT runs (as for the
//      r-index; KATKA's katka-build -A);
//   2. phi (phi(SA[x]) = SA[x-1]) from those samples streams the suffix array from the last row to the
//      first; each entry is mapped to its genome by a predecessor search over the genome starts, and
//      consecutive rows with the same tag are merged into runs (KATKA's katka-tagbuild);
//   3. the runs (cached in prefix.tagruns, or prefix.sp.tagruns with tags by column) give the tag array:
//      run starts, tags (all, sampled, or grammar-compressed) and the RMQ.
// The BWT may have one terminator (BWT(S$)) or one per dataset (from pfp-merge; katka prep -d).
#pragma once
#include "tag_array.hpp"
#include <deque>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <unistd.h>

namespace katka {

// prefix.tbl: one line per genome: index, tag, name, length, start in S
struct genome_table {
    std::vector<u64> start, len;
    std::vector<std::string> tagname;
    void load(const std::string &tbl) {
        std::ifstream in(tbl);
        if (!in) die("cannot open " + tbl);
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream ss(line); std::string f[5];
            for (int i = 0; i < 5; ++i) std::getline(ss, f[i], '\t');
            tagname.push_back(f[1]); len.push_back(std::stoull(f[3])); start.push_back(std::stoull(f[4]));
        }
        if (start.empty()) die("empty " + tbl);
    }
    // the datasets: maximal runs of consecutive genomes with the same tag (their starts in S)
    std::vector<u64> dataset_starts() const {
        std::vector<u64> d;
        for (size_t i = 0; i < start.size(); ++i) if (i == 0 || tagname[i] != tagname[i - 1]) d.push_back(start[i]);
        return d;
    }
    // the tag of each genome: its index, or (by_column) its tag name numbered by first appearance
    std::vector<uint32_t> tags(bool by_column, u64 *ntags = nullptr) const {
        std::map<std::string, uint32_t> id;
        std::vector<uint32_t> t(start.size());
        for (size_t i = 0; i < start.size(); ++i) {
            auto it = id.emplace(tagname[i], (uint32_t)id.size()).first;
            t[i] = by_column ? it->second : (uint32_t)i;
        }
        if (ntags) *ntags = by_column ? id.size() : start.size();
        return t;
    }
};

// LF by runs.  Rows below 2^40 and runs below 2^32: 14 bytes per run.
struct move_lf {
    struct __attribute__((packed)) run {
        uint32_t s0, l0, lf_run; uint8_t s1, l1;
        u64 start() const { return s0 | ((u64)s1 << 32); }
        u64 lf_start() const { return l0 | ((u64)l1 << 32); }
        void set_start(u64 v) { s0 = (uint32_t)v; s1 = (uint8_t)(v >> 32); }
        void set_lf_start(u64 v) { l0 = (uint32_t)v; l1 = (uint8_t)(v >> 32); }
    };
    std::vector<run> runs;                         // plus a sentinel starting at n
    std::vector<unsigned char> head;

    void build(const std::string &bwtfile, const rlbwt &B) {
        if (B.n >= (1ULL << 40)) die("BWT longer than 2^40");
        if (B.R >= (1ULL << 32) - 1) die("too many runs");
        runs.assign(B.R + 1, run{}); head.assign(B.R, 0);
        u64 before[256] = {0};
        for_each_bwt_run(bwtfile, [&](u64 k, u64 s, u64 len, unsigned char c) {
            runs[k].set_start(s); runs[k].set_lf_start(B.C[c] + before[c]); head[k] = c; before[c] += len;
        });
        runs[B.R].set_start(B.n);
        for (u64 k = 0; k < B.R; ++k) runs[k].lf_run = (uint32_t)run_containing(runs[k].lf_start());
    }
    u64 run_containing(u64 row) const {
        return std::upper_bound(runs.begin(), runs.end(), row, [](u64 v, const run &x) { return v < x.start(); }) - runs.begin() - 1;
    }
    void LF(u64 &row, u64 &k) const {
        row = runs[k].lf_start() + (row - runs[k].start());
        k = runs[k].lf_run;
        while (runs[k + 1].start() <= row) ++k;
    }
};

// With K > 1 terminators, rows 0..K-1 start with a terminator, in an order the BWT builder chose; match
// the strings read backwards from them against the datasets' tails.  Returns each dataset's terminator row.
inline std::vector<u64> terminator_rows(const move_lf &M, const std::vector<u64> &doff, u64 K, const std::string &Sfile) {
    if (K == 1) return {0};
    FILE *S = open_or_die(Sfile, "rb");
    std::vector<u64> drow(K, UINT64_MAX);
    std::map<std::string, u64> tails;
    u64 LEN = 64;
    for (;;) {
        tails.clear();
        bool unique = true;
        for (u64 j = 0; j < K; ++j) {
            u64 len = std::min(LEN, doff[j + 1] - doff[j]);
            std::string s(len, 0);
            fseeko(S, doff[j + 1] - len, SEEK_SET);
            if (fread(&s[0], 1, len, S) != len) die("cannot read S");
            std::reverse(s.begin(), s.end());
            if (!tails.emplace(s, j).second) unique = false;
        }
        if (unique) break;
        if (LEN > (1ULL << 26)) die("two datasets have identical tails");
        LEN *= 4;
    }
    fclose(S);
    for (u64 x = 0; x < K; ++x) {
        std::string s;
        u64 row = x, k = M.run_containing(x);
        while (s.size() < LEN && M.head[k] != 1) { s.push_back(M.head[k]); M.LF(row, k); }
        auto it = tails.find(s);
        if (it == tails.end()) die("terminator row " + std::to_string(x) + " matches no dataset");
        drow[it->second] = x;
    }
    return drow;
}

// phi from the samples at the run boundaries, in the coordinates of D = T_0 $ T_1 $ ... T_{K-1} $
struct phi_samples {
    u64 n = 0, R = 0;
    PackedInts esa, vals, extra;    // SA at the last row of each run; the index in esa (or R + index in extra) of SA[x-1]
    EliasFano keys;                 // sorted SA[x] for rows x >= 1 starting a run (and extra terminator rows)
    std::vector<u64> dstart;        // dataset starts in D
    u64 phi(u64 i) const {
        u64 t = keys.count_less(i + 1) - 1, j = keys.access(t), v = vals[t];
        return (v < R ? esa[v] : extra[v - R]) + (i - j);
    }
};

inline phi_samples sample_rindex(const std::string &Sfile, const std::string &bwtfile, const rlbwt &B, const genome_table &G) {
    const u64 N = file_size(Sfile), K = B.C[2] - B.C[1];
    std::vector<u64> doff;
    if (K == 1) doff = {0};
    else {
        doff = G.dataset_starts();
        if (doff.size() != K) die("the BWT has " + std::to_string(K) + " terminators but S has " + std::to_string(doff.size()) + " datasets");
    }
    doff.push_back(N);
    if (B.n != N + K) die("BWT length " + std::to_string(B.n) + " != |S| + " + std::to_string(K));
    move_lf M; M.build(bwtfile, B);
    std::vector<u64> drow = terminator_rows(M, doff, K, Sfile);
    const u64 R = B.R, n = B.n;
    // the samples are spilled to temporary files during the pass, to keep them out of its peak memory
    struct __attribute__((packed)) entry { uint32_t run; u64 sa; };
    const std::string tmp = Sfile + ".samples.tmp";
    FILE *fs = open_or_die(tmp + ".s", "wb"), *fe = open_or_die(tmp + ".e", "wb");
    std::vector<std::pair<u64, u64>> term;
    for (u64 j = K; j-- > 0;) {
        u64 len = doff[j + 1] - doff[j], dpos = doff[j] + j;
        u64 row = drow[j], k = M.run_containing(row);
        for (u64 p = len;; --p) {
            if (row == M.runs[k].start()) { entry x{(uint32_t)k, dpos + p}; fwrite(&x, sizeof x, 1, fs); }
            if (row + 1 == M.runs[k + 1].start()) { entry x{(uint32_t)k, dpos + p}; fwrite(&x, sizeof x, 1, fe); }
            if (p == 0) { term.push_back({row, dpos}); break; }
            M.LF(row, k);
        }
        if (M.head[k] != 1) die("the LF pass over dataset " + std::to_string(j) + " did not end at a terminator: BWT inconsistent");
    }
    fclose(fs); fclose(fe);
    std::sort(term.begin(), term.end());
    std::vector<std::pair<u64, u64>> xk;
    std::vector<u64> extra;
    for (u64 t = 0; t < term.size(); ++t) {
        u64 x = term[t].first;
        if (x == 0 || M.runs[M.run_containing(x)].start() == x) continue;
        if (t == 0 || term[t - 1].first != x - 1) die("internal error: terminator run");
        xk.push_back({term[t].second, R + extra.size()});
        extra.push_back(term[t - 1].second);
    }
    std::vector<move_lf::run>().swap(M.runs);
    auto load = [&](const std::string &f) {
        PackedInts a; a.init(R, PackedInts::width_for(n));
        std::vector<char> seen(R, 0);
        FILE *in = open_or_die(f, "rb");
        std::vector<entry> buf(1 << 16);
        size_t got;
        while ((got = fread(buf.data(), sizeof(entry), buf.size(), in)) > 0)
            for (size_t i = 0; i < got; ++i) { a.set(buf[i].run, buf[i].sa); seen[buf[i].run] = 1; }
        fclose(in); unlink(f.c_str());
        for (u64 k = 0; k < R; ++k) if (!seen[k]) die("internal error: a run without samples");
        return a;
    };
    phi_samples P; P.n = n; P.R = R;
    {
        PackedInts ssa = load(tmp + ".s");
        const u64 m = (R - 1) + xk.size();
        auto key = [&](u64 i) -> u64 { return i < R - 1 ? ssa[i + 1] : xk[i - (R - 1)].first; };
        auto val = [&](u64 i) -> u64 { return i < R - 1 ? i : xk[i - (R - 1)].second; };   // run i ends just before run i + 1
        std::vector<uint32_t> perm(m);
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b) { return key(a) < key(b); });
        std::vector<u64> keys(m);
        for (u64 i = 0; i < m; ++i) keys[i] = key(perm[i]);
        P.keys.build(keys, n);
        std::vector<u64>().swap(keys);
        P.vals.init(m, PackedInts::width_for(R + extra.size()));
        for (u64 i = 0; i < m; ++i) P.vals.set(i, val(perm[i]));
    }
    P.esa = load(tmp + ".e");
    P.extra.build(extra);
    for (u64 j = 0; j < K; ++j) P.dstart.push_back(doff[j] + j);
    return P;
}

// the tag runs, top-down: (first row, tag)
struct tag_runs {
    u64 n = 0;
    std::vector<u64> start;
    std::vector<uint32_t> tag;
    bool load(const std::string &f) {
        FILE *c = fopen(f.c_str(), "rb");
        if (!c) return false;
        char m[8];
        bool ok = fread(m, 1, 8, c) == 8 && !memcmp(m, "KTRUNS01", 8);
        if (ok) { rd(c, n); rvec(c, start); rvec(c, tag); }
        fclose(c);
        return ok;
    }
    void save(const std::string &f) const { FILE *o = open_or_die(f, "wb"); wmagic(o, "KTRUNS01"); wr(o, n); wvec(o, start); wvec(o, tag); fclose(o); }

    void compute(const std::string &Sfile, const std::string &bwtfile, const rlbwt &B, const genome_table &G, bool by_column) {
        phi_samples P = sample_rindex(Sfile, bwtfile, B, G);
        std::vector<uint32_t> gt = G.tags(by_column);
        n = P.n;
        const u64 K = P.dstart.size();
        std::vector<u64> ds(P.dstart); ds.push_back(n);
        std::vector<u64> rs; std::vector<uint32_t> rt;
        u64 d = P.esa[P.R - 1], j = K - 1;                  // SA[n-1], in D
        u64 g = 0, g0 = 0, g1 = G.start.size() > 1 ? G.start[1] : ~0ULL;
        int64_t cur = -1; u64 top = 0;
        for (u64 row = n; row-- > 0;) {
            u64 q;                                          // D position -> S position (a terminator: its dataset's last)
            if (K == 1) q = d < n - 1 ? d : n - 2;
            else {
                if (d < ds[j] || d >= ds[j + 1]) j = std::upper_bound(ds.begin(), ds.end(), d) - ds.begin() - 1;
                q = d - j;
                if (d + 1 == ds[j + 1]) --q;
            }
            if (q < g0 || q >= g1) {
                g = std::upper_bound(G.start.begin(), G.start.end(), q) - G.start.begin() - 1;
                g0 = G.start[g]; g1 = g + 1 < G.start.size() ? G.start[g + 1] : ~0ULL;
            }
            int64_t t = gt[g];
            if (t != cur) { if (cur >= 0) { rs.push_back(top); rt.push_back((uint32_t)cur); } cur = t; }
            top = row;
            if (row) d = P.phi(d);
        }
        rs.push_back(top); rt.push_back((uint32_t)cur);
        std::reverse(rs.begin(), rs.end()); std::reverse(rt.begin(), rt.end());
        start.swap(rs); tag.swap(rt);
        if (start.empty() || start[0] != 0) die("internal error: no tag run at row 0");
    }
};

// Sampling for s > 1 (KATKA's choose_samples): e[i] = the run containing LF of run i's first row; a run
// must be sampled if it is run 0 or e[i] has another tag; otherwise sample greedily from the leaves so
// that every run reaches a sampled one in fewer than s steps, and sample one run on every cycle.
inline std::vector<char> choose_samples(const tag_runs &TR, const rlbwt &B, u64 s) {
    const u64 rho = TR.start.size();
    std::vector<uint32_t> e(rho), indeg(rho, 0);
    std::vector<char> root(rho, 0), samp(rho, 0), done(rho, 0);
    std::vector<uint8_t> need(rho, 0);
    for (u64 i = 0; i < rho; ++i) {
        u64 lf = B.LF(TR.start[i]);
        e[i] = std::upper_bound(TR.start.begin(), TR.start.end(), lf) - TR.start.begin() - 1;
        if (i == 0 || TR.tag[e[i]] != TR.tag[i]) root[i] = 1;
    }
    for (u64 i = 0; i < rho; ++i) if (!root[i]) ++indeg[e[i]];
    auto process = [&](u64 v, bool force) {
        done[v] = 1;
        if (root[v] || force || need[v] + 1u >= s) { samp[v] = 1; return; }
        u64 w = e[v];
        if (need[v] + 1 > need[w]) need[w] = need[v] + 1;
    };
    std::deque<u64> Q;
    for (u64 i = 0; i < rho; ++i) if (!indeg[i]) Q.push_back(i);
    while (!Q.empty()) {
        u64 v = Q.front(); Q.pop_front();
        process(v, false);
        if (!root[v] && --indeg[e[v]] == 0) Q.push_back(e[v]);
    }
    for (u64 c = 0; c < rho; ++c) {
        if (done[c]) continue;
        process(c, true);
        for (u64 v = e[c]; v != c; v = e[v]) process(v, false);
    }
    return samp;
}

struct tag_build_options {
    u64 s = 1;
    bool rmq = true;
    bool grammar = false;
    std::string bigrepair;          // grammar from bigrepair -i output with this prefix
};

inline void build_tag_array(const tag_runs &TR, const rlbwt *B, const tag_build_options &o, tag_array &X) {
    const u64 rho = TR.start.size();
    X.n = TR.n; X.rho = rho; X.s = o.s;
    X.RB.build(TR.start, TR.n);
    if (o.rmq) {
        std::vector<uint64_t> last(1 << 16, 0);
        u64 i = 0;
        X.rmq.build(rho, [&]() {
            uint32_t t = TR.tag[i];
            if (t >= last.size()) last.resize(t + 1, 0);
            u64 c = last[t]; last[t] = ++i;
            return c;
        });
        X.has_rmq = true;
    }
    if (o.grammar || !o.bigrepair.empty()) {
        if (o.s != 1) die("grammar-compressed tags need s = 1");
        uint32_t alpha = 0;
        for (uint32_t t : TR.tag) alpha = std::max(alpha, t + 1);
        grammar_rules Gr = o.bigrepair.empty() ? repair(TR.tag, alpha) : read_bigrepair(o.bigrepair);
        X.G.reset(new tag_grammar);
        X.G->build(Gr);
        if (X.G->rho != rho) die("the grammar generates " + std::to_string(X.G->rho) + " runs, not " + std::to_string(rho));
        u64 i = 0, bad = 0;
        X.G->decode(0, rho, [&](uint32_t g) { if (g != TR.tag[i]) ++bad; ++i; });
        if (bad) die("the grammar gives " + std::to_string(bad) + " runs the wrong tag");
        fprintf(stderr, "grammar: %zu rules, top level %zu symbols, %.3f MB (tags as an array: %.3f MB)\n", Gr.rules.size(), Gr.C.size(),
                X.G->bytes() / 1e6, rho * PackedInts::width_for(alpha) / 8e6);
        return;
    }
    if (o.s == 1) { X.L.build_from(TR.tag); return; }
    std::vector<char> samp = choose_samples(TR, *B, o.s);
    X.SB.init(rho);
    std::vector<uint32_t> st;
    for (u64 i = 0; i < rho; ++i) if (samp[i]) { X.SB.set(i); st.push_back(TR.tag[i]); }
    X.SB.finish();
    X.L.build_from(st);
    X.bwt = B;
    u64 bad = 0, step = o.s > 4 ? 97 : 1;
    for (u64 i = 0; i < rho; i += step) if (X.tag(i) != TR.tag[i]) ++bad;
    if (bad) die(std::to_string(bad) + " runs get the wrong tag");
    fprintf(stderr, "s=%lu: %lu of %lu runs sampled\n", o.s, st.size(), rho);
}

}  // namespace katka
