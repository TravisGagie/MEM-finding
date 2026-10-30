// katka-tagbuild: build the run-length compressed tag array (tag_array.hpp) without a suffix array (Section 6,
// "Construction"; Appendix A, "Indexes").
//
// usage: katka-tagbuild prefix s out.tag [-n] [-S]
//   reads prefix.tbl (katka-prep), prefix.rix (katka-build -A) and, for s > 1, prefix.rlbwt (katka-build -G; for LF)
//   s   sampling parameter (1: store every run's tag, as in the paper)
//   -n  no RMQ: counting only (listing needs out.tag.rmq)
//   -S  a position's tag is the tag column of its genome in prefix.tbl (the species or genus), numbered in
//       order of first appearance; otherwise it is the genome's index
// writes out.tag and, unless -n, out.tag.rmq.
//
// Stage 1, the tag runs: phi (katka-build -A) gives SA[x-1] from SA[x], so starting from SA[n-1] we stream the
// suffix array from the last row to the first, map each entry to its genome by a predecessor search over the
// genome starts, and merge consecutive rows with the same tag into runs.  The memory is that of the r-index
// samples plus the runs; the runs are spilled to disk while the samples are in memory.  The runs are cached
// in prefix.tagruns (prefix.sp.tagruns with -S), so other values of s reuse them.
// Stage 2: the tag array's bitvector of run starts, its tags (or a sample of them, s > 1) and the RMQ.
#include "tag_array.hpp"
#include <deque>
#include <map>
#include <sstream>

using namespace katka;
using std::string;
using std::vector;

// The r-index's phi, from prefix.rix (see build.cpp): phi(SA[x]) = SA[x-1], in D coordinates
//     D = T_0 $ T_1 $ ... T_{K-1} $  (one string per dataset, each with its terminator).
struct phi_index {
    u64 n = 0, R = 0;
    sdsl::int_vector<> esa, vals, extra;
    sparse_bv keys, dstart;
    void load(const string &file) {
        std::ifstream in(file, std::ios::binary);
        if (!in) die("katka-tagbuild: " + file + " not found (build it with katka-build -A)");
        in.read((char *)&n, 8);
        { rlbwt B; B.load(in); R = B.R; }          // the RLBWT is in the file, but not needed here
        esa.load(in); keys.load(in); vals.load(in); extra.load(in); dstart.load(in);
        if (!in) die("cannot read " + file);
    }
    u64 last_row_sa() const { return esa[R - 1]; }  // SA[n-1], the last row of the last run
    u64 phi(u64 i) const {
        u64 t = keys.rank(i + 1), j = keys.select(t), v = vals[t - 1];   // j = the sampled SA[x] preceding i
        return (v < R ? esa[v] : extra[v - R]) + (i - j);
    }
};

// prefix.tbl: each genome's start in S and its tag
struct genome_tags {
    vector<u64> start;
    vector<uint16_t> tag;
    genome_tags(const string &tbl, bool by_column) {
        std::ifstream in(tbl);
        if (!in) die("cannot open " + tbl);
        std::map<string, int> id;
        string line;
        while (std::getline(in, line)) {
            std::istringstream ss(line); string f[5];
            for (int i = 0; i < 5; ++i) std::getline(ss, f[i], '\t');
            auto it = id.emplace(f[1], (int)id.size()).first;
            start.push_back(std::stoull(f[4]));
            tag.push_back((uint16_t)(by_column ? it->second : start.size() - 1));
        }
        if (start.empty() || (start.size() > 1 && start[1] == 0)) die("katka-tagbuild: bad " + tbl);
        if (by_column) fprintf(stderr, "tags: %zu over %zu genomes\n", id.size(), start.size());
    }
};

template <class T> struct tag_runs {
    u64 n = 0;
    vector<T> start;            // first row of each run
    vector<uint16_t> tag;

    bool load(const string &cache) {
        std::ifstream c(cache, std::ios::binary);
        if (!c) return false;
        u64 r; c.read((char *)&n, 8); c.read((char *)&r, 8);
        start.resize(r); tag.resize(r);
        c.read((char *)start.data(), sizeof(T) * r); c.read((char *)tag.data(), 2 * r);
        return (bool)c;
    }
    void save(const string &cache) const {
        std::ofstream o(cache, std::ios::binary); u64 r = start.size();
        o.write((char *)&n, 8); o.write((char *)&r, 8);
        o.write((char *)start.data(), sizeof(T) * r); o.write((char *)tag.data(), 2 * r);
    }

    void compute(const string &prefix, const genome_tags &G, const string &cache) {
        const string frows = cache + ".rows.tmp", ftags = cache + ".tags.tmp";
        u64 nruns = 0;
        {
            phi_index P; P.load(prefix + ".rix");
            n = P.n;
            const u64 K = P.dstart.ones();         // datasets, each followed by its terminator in D
            vector<u64> ds(K + 1);
            for (u64 j = 0; j < K; ++j) ds[j] = P.dstart.select(j + 1);
            ds[K] = n;
            fprintf(stderr, "r-index samples loaded; streaming the suffix array\n");
            FILE *fr = fopen(frows.c_str(), "wb"), *ft = fopen(ftags.c_str(), "wb");
            if (!fr || !ft) die("katka-tagbuild: cannot write temporary files next to " + cache);
            vector<T> brow; vector<uint16_t> btag;
            auto emit = [&](T row, uint16_t t) {   // a run, found bottom-up
                brow.push_back(row); btag.push_back(t); ++nruns;
                if (brow.size() == (1u << 20) ) { fwrite(brow.data(), sizeof(T), brow.size(), fr); fwrite(btag.data(), 2, btag.size(), ft); brow.clear(); btag.clear(); }
            };
            u64 d = P.last_row_sa(), j = K - 1;    // SA[row] in D; the dataset containing it
            u64 g = 0, g0 = 0, g1 = G.start.size() > 1 ? G.start[1] : ~0ULL;   // its genome [g0, g1) in S
            int cur = -1; T top = 0;               // the current run's tag and its topmost row so far
            for (u64 row = n; row-- > 0;) {
                // D position -> S position; a terminator gets its dataset's last character
                u64 q;
                if (K == 1) q = d < n - 1 ? d : n - 2;
                else {
                    if (d < ds[j] || d >= ds[j + 1]) j = std::upper_bound(ds.begin(), ds.end(), d) - ds.begin() - 1;
                    q = d - j;
                    if (d + 1 == ds[j + 1]) --q;
                }
                if (q < g0 || q >= g1) {           // S position -> genome
                    g = std::upper_bound(G.start.begin(), G.start.end(), q) - G.start.begin() - 1;
                    g0 = G.start[g]; g1 = g + 1 < G.start.size() ? G.start[g + 1] : ~0ULL;
                }
                int t = G.tag[g];
                if (t != cur) { if (cur >= 0) emit(top, (uint16_t)cur); cur = t; }
                top = (T)row;
                if (row) d = P.phi(d);
            }
            emit(top, (uint16_t)cur);
            fwrite(brow.data(), sizeof(T), brow.size(), fr); fwrite(btag.data(), 2, btag.size(), ft);
            fclose(fr); fclose(ft);
        }   // the r-index samples are freed here
        // read the runs back, reversing them into top-down order
        start.resize(nruns); tag.resize(nruns);
        FILE *fr = fopen(frows.c_str(), "rb"), *ft = fopen(ftags.c_str(), "rb");
        vector<T> brow(1 << 20); vector<uint16_t> btag(1 << 20);
        u64 i = nruns;
        for (size_t k; (k = fread(brow.data(), sizeof(T), brow.size(), fr)) > 0;) {
            if (fread(btag.data(), 2, k, ft) != k) die("katka-tagbuild: short read of temporary file");
            for (size_t q = 0; q < k; ++q) { --i; start[i] = brow[q]; tag[i] = btag[q]; }
        }
        fclose(fr); fclose(ft); std::remove(frows.c_str()); std::remove(ftags.c_str());
        if (i != 0) die("katka-tagbuild: temporary file has the wrong length");
        if (start.empty() || start[0] != 0) die("katka-tagbuild: internal error: no run at row 0");
    }
};

// Sampling for s > 1: e[i] = the run containing LF of run i's first row; a run is a root (must be sampled)
// if it is run 0 or e[i] has a different tag.  On the graph i -> e[i] (one out-edge, except roots), choose
// samples greedily from the leaves, so that every run reaches a sampled one in fewer than s steps, and
// sample one run on every cycle.
template <class T> static vector<char> choose_samples(const tag_runs<T> &TR, const rlbwt &B, u64 s, u64 &nroot, u64 &ncyc) {
    const u64 rho = TR.start.size();
    vector<uint32_t> e(rho), indeg(rho, 0);
    vector<char> root(rho, 0), samp(rho, 0), done(rho, 0);
    vector<uint8_t> need(rho, 0);                // steps from the farthest unsampled run that reaches this one
    nroot = ncyc = 0;
    for (u64 i = 0; i < rho; ++i) {
        u64 lf = B.LF(TR.start[i]);
        e[i] = std::upper_bound(TR.start.begin(), TR.start.end(), (T)lf) - TR.start.begin() - 1;
        if (i == 0 || TR.tag[e[i]] != TR.tag[i]) { root[i] = 1; ++nroot; }
    }
    for (u64 i = 0; i < rho; ++i) if (!root[i]) ++indeg[e[i]];
    auto process = [&](u64 v, bool force) {
        done[v] = 1;
        if (root[v] || force || need[v] + 1 >= s) { samp[v] = 1; return; }
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
    for (u64 c = 0; c < rho; ++c) {              // what is left lies on cycles
        if (done[c]) continue;
        ++ncyc; process(c, true);
        for (u64 v = e[c]; v != c; v = e[v]) process(v, false);
    }
    return samp;
}

template <class T> static int run(const string &prefix, u64 s, const string &out, bool no_rmq, bool by_column) {
    const string cache = prefix + (by_column ? ".sp.tagruns" : ".tagruns");
    double t0 = now();
    tag_runs<T> TR;
    if (!TR.load(cache)) {
        genome_tags G(prefix + ".tbl", by_column);
        TR.compute(prefix, G, cache);
        TR.save(cache);
    }
    const u64 rho = TR.start.size();
    fprintf(stderr, "n=%lu runs=%lu (stage 1 %.0f s)\n", TR.n, rho, now() - t0);

    tag_array X; X.n = TR.n; X.rho = rho; X.s = s;
    { sdsl::sd_vector_builder b(TR.n, rho); for (T r : TR.start) b.set(r); X.RB = sdsl::sd_vector<>(b); }
    sdsl::util::init_support(X.RBr, &X.RB); sdsl::util::init_support(X.RBs, &X.RB);
    if (!no_rmq) {                               // C[i] = 1 + the previous run with the same tag (0 if none)
        sdsl::int_vector<> C(rho, 0, sdsl::bits::hi(rho + 1) + 1);
        vector<uint32_t> last(1 << 16, 0);
        for (u64 i = 0; i < rho; ++i) { C[i] = last[TR.tag[i]]; last[TR.tag[i]] = i + 1; }
        X.rmq = sdsl::rmq_succinct_sct<>(&C); X.has_rmq = true;
    }
    u64 nroot = 0, ncyc = 0, nsamp = rho;
    if (s == 1) {
        X.L = sdsl::int_vector<>(rho, 0, 14);
        for (u64 i = 0; i < rho; ++i) X.L[i] = TR.tag[i];
    } else {
        rlbwt_file Z; Z.load(prefix + ".rlbwt");
        vector<char> samp = choose_samples(TR, Z.bwt, s, nroot, ncyc);
        X.SB = sdsl::bit_vector(rho, 0);
        nsamp = 0;
        for (u64 i = 0; i < rho; ++i) if (samp[i]) { X.SB[i] = 1; ++nsamp; }
        sdsl::util::init_support(X.SBr, &X.SB);
        X.L = sdsl::int_vector<>(nsamp, 0, 14);
        for (u64 i = 0, j = 0; i < rho; ++i) if (samp[i]) X.L[j++] = TR.tag[i];
        X.bwt = &Z.bwt;                          // check the tags (every run, or every 97th for s > 4)
        u64 bad = 0, step = s > 4 ? 97 : 1;
        for (u64 i = 0; i < rho; i += step) if (X.tag(i) != TR.tag[i]) ++bad;
        if (bad) die("ERROR: " + std::to_string(bad) + " runs get the wrong tag");
    }
    X.save(out);
    fprintf(stderr, "s=%lu: runs %lu, roots %lu, cycles %lu, sampled %lu (%.1f%%) | bytes %lu (%.2f bits/run) | %.0f s\n",
            s, rho, nroot, ncyc, nsamp, 100.0 * nsamp / rho, X.bytes(), 8.0 * X.bytes() / rho, now() - t0);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) die("usage: katka-tagbuild prefix s out.tag [-n] [-S]");
    bool no_rmq = false, by_column = false;
    for (int i = 4; i < argc; ++i) {
        string o = argv[i];
        if (o == "-n") no_rmq = true; else if (o == "-S") by_column = true; else die("unknown option " + o);
    }
    const string prefix = argv[1];
    // n decides the width of the run starts: 32 bits if n < 2^32
    u64 n = 0;
    { std::ifstream c(prefix + (by_column ? ".sp.tagruns" : ".tagruns"), std::ios::binary); if (c) c.read((char *)&n, 8); }
    if (!n) { std::ifstream r(prefix + ".rix", std::ios::binary); if (r) r.read((char *)&n, 8); }
    u64 s = std::stoull(argv[2]);
    if (s < 1 || s > 255) die("s must be between 1 and 255");
    return n >> 32 ? run<uint64_t>(prefix, s, argv[3], no_rmq, by_column) : run<uint32_t>(prefix, s, argv[3], no_rmq, by_column);
}
