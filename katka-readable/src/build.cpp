// katka-build: from the text S (katka-prep) and its BWT (Big-BWT, or pfp-merge for one string per dataset), build
//   -G  out.rlbwt  the RLBWT (rlbwt.hpp) and the strand starts: all that classification needs besides the tags
//   -A  out.rix    the r-index's suffix-array samples at the BWT run boundaries, with phi (Appendix A,
//                  "Indexes"), from which katka-tagbuild streams the suffix array to build the tag array
//
// usage: katka-build [-G] [-A] S bwt [- -] S.tbl out-prefix
//   (the two "-" stand for the LZ77 parses that the original code's rz-index needs; they are ignored)
//
// The BWT has one byte per row and either one terminator (BWT(S$)) or one per dataset, the datasets being
// the maximal runs of consecutive genomes with the same tag in S.tbl (katka-prep -d).  Terminators may be bytes
// 0, 1 or 2; the RLBWT reads them all as 1.
//
// No suffix array is built.  We walk LF over each dataset from its end to its start, which visits every BWT
// row once and tells us the text position of each row as we pass it; the rows at the start and end of BWT
// runs are the samples.  LF uses "move"-style arrays of O(r) words: for each run, its start, the LF of its
// start, and the run containing that.
#include "rlbwt.hpp"
#include <map>
#include <numeric>
#include <sstream>
#include <unistd.h>

using namespace katka;
using std::string;
using std::vector;

// S.tbl: the start in S of every strand (G_i and rc(G_i)) and of every dataset
struct text_layout {
    u64 N = 0;
    vector<u64> strand_start;
    vector<u64> dataset_start;
};
static text_layout read_table(const string &tbl, u64 N) {
    text_layout T; T.N = N;
    std::ifstream t(tbl);
    if (!t) die("cannot open " + tbl);
    string line, last;
    u64 pos = 0;
    while (std::getline(t, line)) {
        std::istringstream ss(line);
        string idx, tag, path; u64 len, start;
        std::getline(ss, idx, '\t'); std::getline(ss, tag, '\t'); std::getline(ss, path, '\t');
        ss >> len >> start;
        if (start != pos) die("table/offset mismatch");
        if (T.dataset_start.empty() || tag != last) { T.dataset_start.push_back(pos); last = tag; }
        T.strand_start.push_back(pos); pos += len + 1;      // G_i X
        T.strand_start.push_back(pos); pos += len + 1;      // rc(G_i) X
    }
    if (pos != N) die("table does not match |S| (" + std::to_string(pos) + " vs " + std::to_string(N) + ")");
    return T;
}

// LF by runs.  Rows are below 2^40 and runs below 2^32: 14 bytes per run.
struct move_lf {
    struct __attribute__((packed)) run {
        uint32_t s0, l0, lf_run; uint8_t s1, l1;   // start and LF(start), 40 bits each; the run containing LF(start)
        u64 start() const { return s0 | ((u64)s1 << 32); }
        u64 lf_start() const { return l0 | ((u64)l1 << 32); }
        void set_start(u64 v) { s0 = (uint32_t)v; s1 = (uint8_t)(v >> 32); }
        void set_lf_start(u64 v) { l0 = (uint32_t)v; l1 = (uint8_t)(v >> 32); }
    };
    vector<run> runs;                              // plus a sentinel starting at n
    vector<unsigned char> head;

    void build(const string &bwtfile, const rlbwt &B) {
        if (B.n >= (1ULL << 40)) die("BWT longer than 2^40");
        if (B.R >= (1ULL << 32) - 1) die("too many runs");
        runs.assign(B.R + 1, run{}); head.assign(B.R, 0);
        u64 before[256] = {0};
        rlbwt::for_each_run(bwtfile, [&](u64 k, u64 s, u64 len, unsigned char c) {
            runs[k].set_start(s);
            runs[k].set_lf_start(B.C[c] + before[c]);       // LF of the run's first row
            head[k] = c; before[c] += len;
        });
        runs[B.R].set_start(B.n);
        for (u64 k = 0; k < B.R; ++k) runs[k].lf_run = (uint32_t)run_containing(runs[k].lf_start());
    }
    u64 run_containing(u64 row) const {
        return std::upper_bound(runs.begin(), runs.end(), row, [](u64 v, const run &x) { return v < x.start(); }) - runs.begin() - 1;
    }
    // (row, run k containing it) -> (LF(row), the run containing that)
    void LF(u64 &row, u64 &k) const {
        row = runs[k].lf_start() + (row - runs[k].start());
        k = runs[k].lf_run;
        while (runs[k + 1].start() <= row) ++k;
    }
    void free() { vector<run>().swap(runs); vector<unsigned char>().swap(head); }
};

// With K > 1 terminators, rows 0..K-1 start with a terminator, in an order that the BWT builder chose.
// Reading the BWT backwards from such a row spells its dataset's text backwards; we match those strings
// against the datasets' tails.  Returns, for each dataset, its terminator row.
static vector<u64> terminator_rows(const move_lf &M, const vector<u64> &doff, u64 K, FILE *S) {
    if (K == 1) return {0};
    vector<u64> drow(K, UINT64_MAX);
    std::map<string, u64> tails;
    u64 LEN = 64;
    for (;;) {                                     // tail length at which all datasets differ
        tails.clear();
        bool unique = true;
        for (u64 j = 0; j < K; ++j) {
            u64 len = std::min(LEN, doff[j + 1] - doff[j]);
            string s(len, 0);
            fseeko(S, doff[j + 1] - len, SEEK_SET);
            if (fread(&s[0], 1, len, S) != len) die("cannot read S");
            std::reverse(s.begin(), s.end());
            if (!tails.emplace(s, j).second) unique = false;
        }
        if (unique) break;
        if (LEN > (1ULL << 26)) die("two datasets have identical tails");
        LEN *= 4;
    }
    for (u64 x = 0; x < K; ++x) {
        string s;
        u64 row = x, k = M.run_containing(x);
        while (s.size() < LEN && M.head[k] != 1) { s.push_back(M.head[k]); M.LF(row, k); }
        auto it = tails.find(s);
        if (it == tails.end()) die("terminator row " + std::to_string(x) + " matches no dataset");
        if (drow[it->second] != UINT64_MAX) die("two terminator rows match dataset " + std::to_string(it->second));
        drow[it->second] = x;
    }
    return drow;
}

// Suffix-array samples are spilled to a temporary file as the LF pass finds them, and read back into an
// array (indexed by run) only after the LF arrays are freed, to keep them out of the pass's peak memory.
struct sample_spill {
    struct __attribute__((packed)) entry { uint32_t run; u64 sa; };
    string file; FILE *f = nullptr;
    explicit sample_spill(const string &name) : file(name) {
        f = fopen(file.c_str(), "wb");
        if (!f) die("cannot write " + file);
        setvbuf(f, nullptr, _IOFBF, 1 << 24);
    }
    void add(u64 run, u64 sa) { entry e{(uint32_t)run, sa}; fwrite(&e, sizeof e, 1, f); }
    void close() { fclose(f); }
    sdsl::int_vector<> load(u64 R, u64 n) {
        sdsl::int_vector<> a(R, 0, sdsl::bits::hi(n) + 1);
        FILE *fp = fopen(file.c_str(), "rb");
        vector<entry> buf(1 << 20);
        size_t got, total = 0;
        while ((got = fread(buf.data(), sizeof(entry), buf.size(), fp)) > 0) {
            for (size_t i = 0; i < got; ++i) a[buf[i].run] = buf[i].sa;
            total += got;
        }
        fclose(fp); unlink(file.c_str());
        if (total != R) die("internal error: " + std::to_string(total) + " samples in " + file + ", expected " + std::to_string(R));
        return a;
    }
};

// The r-index's samples and phi (Gagie, Navarro and Prezza), in the coordinates of
//     D = T_0 $ T_1 $ ... T_{K-1} $   (the datasets, each with its terminator).
// phi(SA[x]) = SA[x-1].  For a row x >= 1 that starts a BWT run, row x-1 ends the previous run, so storing
// SA at the run ends (esa) and, sorted by SA[x], the pairs (SA[x], index of SA[x-1] in esa) (keys, vals)
// gives phi(i) = SA[x-1] + (i - SA[x]) for the predecessor SA[x] of i.  The few terminator rows that do not
// start a run are sampled too (their SA[x-1] in extra), so that phi never runs across a dataset boundary.
static void build_rindex(const string &out, const rlbwt &B, move_lf &M, const vector<u64> &doff, const vector<u64> &drow, u64 K) {
    const u64 R = B.R, n = B.n;
    sample_spill starts(out + ".ssa.tmp"), ends(out + ".esa.tmp");
    vector<std::pair<u64, u64>> term;            // (row, SA) of the rows starting a dataset
    double t0 = now();
    for (u64 j = K; j-- > 0;) {                  // LF over dataset j, from its terminator to its first character
        u64 len = doff[j + 1] - doff[j], dpos = doff[j] + j;   // dataset j starts at dpos in D
        u64 row = drow[j], k = M.run_containing(row);
        for (u64 p = len;; --p) {                // row is the row of D position dpos + p
            if (row == M.runs[k].start()) starts.add(k, dpos + p);
            if (row + 1 == M.runs[k + 1].start()) ends.add(k, dpos + p);
            if (p == 0) { term.push_back({row, dpos}); break; }
            M.LF(row, k);
        }
        if (M.head[k] != 1) die("the LF pass over dataset " + std::to_string(j) + " did not end at a terminator: BWT inconsistent");
    }
    starts.close(); ends.close();
    fprintf(stderr, "LF pass: %.0f s\n", now() - t0);

    // extra samples for terminator rows that do not start a run (their predecessor row is also one)
    std::sort(term.begin(), term.end());
    vector<std::pair<u64, u64>> xk;              // (SA[x], index in esa or R + index in extra)
    vector<u64> extra;
    for (u64 t = 0; t < term.size(); ++t) {
        u64 x = term[t].first;
        if (x == 0 || M.runs[M.run_containing(x)].start() == x) continue;
        if (t == 0 || term[t - 1].first != x - 1) die("internal error: terminator run");
        xk.push_back({term[t].second, R + extra.size()});
        extra.push_back(term[t - 1].second);
    }
    M.free();

    sdsl::int_vector<> ssa = starts.load(R, n);  // SA at the first row of each run
    const u64 m = (R - 1) + xk.size();
    auto key = [&](u64 i) -> u64 { return i < R - 1 ? (u64)ssa[i + 1] : xk[i - (R - 1)].first; };
    auto val = [&](u64 i) -> u64 { return i < R - 1 ? i : xk[i - (R - 1)].second; };   // run i ends just before run i + 1
    vector<uint32_t> perm(m);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b) { return key(a) < key(b); });
    sparse_bv keys;
    sdsl::int_vector<> vals(m, 0, sdsl::bits::hi(R + extra.size()) + 1);
    {
        sdsl::sd_vector_builder kb(n, m);
        for (u64 i = 0; i < m; ++i) { kb.set(key(perm[i])); vals[i] = val(perm[i]); }
        vector<uint32_t>().swap(perm);
        sdsl::int_vector<>().swap(ssa);
        keys.build_from(kb);
    }
    sdsl::int_vector<> esa = ends.load(R, n);    // SA at the last row of each run
    sdsl::int_vector<> ext(extra.size(), 0, sdsl::bits::hi(n) + 1);
    for (u64 i = 0; i < extra.size(); ++i) ext[i] = extra[i];
    vector<u64> ds(K);
    for (u64 j = 0; j < K; ++j) ds[j] = doff[j] + j;
    sparse_bv dstart; dstart.build(ds, n);

    // out.rix: n, the RLBWT, esa, keys, vals, extra, dataset starts in D (the original code's format)
    std::ofstream o(out + ".rix", std::ios::binary);
    o.write((char *)&n, 8);
    u64 bb = B.serialize(o);
    u64 sb = esa.serialize(o) + keys.serialize(o) + vals.serialize(o) + ext.serialize(o) + dstart.serialize(o);
    printf("r-index bytes=%lu rlbwt=%lu samples=%lu\n", 8 + bb + sb, bb, sb);
}

int main(int argc, char **argv) {
    bool rl = false, ri = false;
    int opt;
    const char *usage = "usage: katka-build [-G] [-A] S bwt [- -] S.tbl out-prefix\n";
    while ((opt = getopt(argc, argv, "GA")) != -1) {
        if (opt == 'G') rl = true; else if (opt == 'A') ri = true;
        else { fputs(usage, stderr); fputs("(the original code's LZ77 grids are not in this version)\n", stderr); return 1; }
    }
    int np = argc - optind;
    if ((np != 4 && np != 6) || (!rl && !ri)) { fputs(usage, stderr); return 1; }
    string Sfile = argv[optind], bwtfile = argv[optind + 1], tblfile = argv[optind + np - 2], out = argv[optind + np - 1];
    double t0 = now();

    FILE *S = fopen(Sfile.c_str(), "rb");
    if (!S) die("cannot open " + Sfile);
    fseeko(S, 0, SEEK_END);
    const u64 N = ftello(S);
    text_layout T = read_table(tblfile, N);

    rlbwt_file Z; Z.N = N; Z.nstr = T.strand_start.size();
    Z.strand_starts.build(T.strand_start, N);
    rlbwt &B = Z.bwt;
    B.build(bwtfile);
    const u64 K = B.C[2] - B.C[1];               // number of terminators = number of strings in the BWT
    vector<u64> doff;                            // the datasets' starts in S, then N
    if (K == 1) doff = {0};
    else if (K == T.dataset_start.size()) doff = T.dataset_start;
    else die("the BWT has " + std::to_string(K) + " terminators, but S has 1 text and " + std::to_string(T.dataset_start.size()) + " datasets");
    doff.push_back(N);
    if (B.n != N + K) die("BWT length " + std::to_string(B.n) + " != |S| + " + std::to_string(K));
    fprintf(stderr, "n=%lu r=%lu strings=%lu\n", B.n, B.R, K);

    if (rl) {
        std::ofstream o(out + ".rlbwt", std::ios::binary);
        u64 bytes = Z.serialize(o);
        printf("N=%lu genomes=%lu r=%lu .rlbwt bytes=%lu\n", N, Z.nstr / 2, B.R, bytes);
    }
    if (ri) {
        move_lf M; M.build(bwtfile, B);
        vector<u64> drow = terminator_rows(M, doff, K, S);
        build_rindex(out, B, M, doff, drow, K);
    }
    fclose(S);
    printf("build time %.1fs\n", now() - t0);
    return 0;
}
