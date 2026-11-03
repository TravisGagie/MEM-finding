// katka-test: checks the tool against brute force, on small data.
//
//   katka-test unit                                  bit structures, RMQ, grammar (random inputs)
//   katka-test mems PREFIX reads [options]           MEMs, intervals and tag answers against a suffix array
//       --tags FILE     also check the tag answers (count, list) of every MEM against brute force
//       --kebab k       also check KeBaB modes, with PREFIX.k<k>.kbf
//       --twolevel      also check two-level search (PREFIX.lvl), with and without block jumps
//       --hpc           homopolymer-compress the reads first (index built with katka prep -U)
//       --random n      also n random patterns made from S (mutated substrings, half reverse-complemented)
//
// The suffix array of S$ (libsais) gives the brute-force answers: R[a..e) occurs iff its SA interval is
// non-empty, MEMs are the (a, e(a)) with e(a) = the largest e such that R[a..e) occurs and e(a-1) < e(a),
// and the tags of a row are those of the genome containing SA[row].
#include "mems.hpp"
#include "tag_build.hpp"
#include "rmq.hpp"
#include "libsais.h"
#include <map>
#include <random>
#include <set>

using namespace katka;
using std::string;
using std::vector;

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++failures; if (failures < 20) { fprintf(stderr, "FAIL: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); } } } while (0)

// ------------------------------------------------------------------------------------------------
static int unit() {
    std::mt19937_64 rng(7);
    for (int it = 0; it < 200; ++it) {
        uint64_t n = 1 + rng() % 3000;
        vector<uint64_t> A(n);
        uint64_t mx = 1 + rng() % 50;
        for (auto &x : A) x = rng() % mx;
        uint64_t k = 0;
        succinct_rmq R; R.build(n, [&] { return A[k++]; });
        for (int q = 0; q < 300; ++q) {
            uint64_t l = rng() % n, r = rng() % n;
            if (l > r) std::swap(l, r);
            uint64_t m = R(l, r), mn = *std::min_element(A.begin() + l, A.begin() + r + 1);
            CHECK(m >= l && m <= r && A[m] == mn, "rmq n=%lu [%lu,%lu]", n, l, r);
        }
        vector<uint64_t> v;
        for (uint64_t i = 0; i < n; ++i) if (rng() % 3 == 0) v.push_back(i * 7 + rng() % 7);
        EliasFano E; E.build(v, n * 7 + 7);
        for (uint64_t i = 0; i < v.size(); ++i) CHECK(E.access(i) == v[i], "EF access");
        for (uint64_t x = 0; x < n * 7 + 7; x += 3) CHECK(E.count_less(x) == (uint64_t)(std::lower_bound(v.begin(), v.end(), x) - v.begin()), "EF count_less");
        BitRank B; B.init(n);
        vector<int> bb(n);
        for (uint64_t i = 0; i < n; ++i) if (rng() % 2) { B.set(i); bb[i] = 1; }
        B.finish();
        uint64_t c = 0, z = 0;
        for (uint64_t i = 0; i < n; ++i) {
            CHECK(B.rank1(i) == c, "rank1");
            if (bb[i]) { CHECK(B.select1(c) == i, "select1"); ++c; } else { CHECK(B.select0(z) == i, "select0"); ++z; }
        }
        vector<uint32_t> s(n);
        uint32_t sig = 1 + rng() % 7;
        for (auto &x : s) x = rng() % sig;
        WaveletMatrix W; W.build(s, std::max<uint32_t>(sig, 2));
        vector<uint64_t> cnt(sig, 0);
        for (uint64_t i = 0; i < n; ++i) {
            CHECK(W.access(i) == s[i], "wm access");
            for (uint32_t a = 0; a < sig; ++a) CHECK(W.rank(a, i) == cnt[a], "wm rank");
            cnt[s[i]]++;
        }
        EFRank F; F.build(s, sig);
        vector<uint64_t> cc(sig, 0);
        for (uint64_t i = 0; i <= n; ++i) {
            for (uint32_t a = 0; a < sig; ++a) { uint64_t r1, r2; F.rank2(a, i, i, r1, r2); CHECK(r1 == cc[a], "EFRank"); }
            if (i < n) cc[s[i]]++;
        }
        // grammar: repetitive sequences of tags
        vector<uint32_t> t;
        while (t.size() < n) {
            if (t.size() > 10 && rng() % 2) { uint64_t a = rng() % t.size(), l = 1 + rng() % std::min<uint64_t>(50, t.size() - a); for (uint64_t q = 0; q < l; ++q) t.push_back(t[a + q]); }
            else t.push_back(rng() % 20);
        }
        grammar_rules Gr = repair(t, 20);
        tag_grammar TG; TG.build(Gr);
        CHECK(TG.rho == t.size(), "grammar length");
        for (int q = 0; q < 50; ++q) {
            uint64_t a = rng() % t.size(), cnt2 = 1 + rng() % (t.size() - a), i = a;
            TG.decode(a, cnt2, [&](uint32_t g) { CHECK(g == t[i], "grammar decode at %lu", i); ++i; });
        }
    }
    printf("unit tests: %s\n", failures ? "FAILED" : "passed");
    return failures ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------
struct brute {
    string S;                       // with the terminator (byte 1) at the end
    vector<int32_t> SA;
    vector<uint32_t> rowtag;        // tag of each row
    void build(const string &prefix) {
        u64 n = file_size(prefix + ".S");
        S.resize(n);
        FILE *f = open_or_die(prefix + ".S", "rb");
        if (fread(&S[0], 1, n, f) != n) die("read S");
        fclose(f);
        SA.resize(n + 1);
        vector<int32_t> sa(n);
        libsais((const uint8_t *)S.data(), sa.data(), (int32_t)n, 0, nullptr);
        SA[0] = (int32_t)n;
        for (u64 i = 0; i < n; ++i) SA[i + 1] = sa[i];
        S.push_back(1);
    }
    // the SA interval of X (half-open)
    interval find(const char *X, u64 len) const {
        auto cmp_lt = [&](int32_t p, bool upper) {       // suffix p < X (or <= X as a prefix, for upper)
            int r = S.compare(p, len, X, len);
            if (r == 0 && (u64)(S.size() - p) < len) r = -1;
            return upper ? r <= 0 : r < 0;
        };
        u64 lo = std::partition_point(SA.begin(), SA.end(), [&](int32_t p) { return cmp_lt(p, false); }) - SA.begin();
        u64 hi = std::partition_point(SA.begin(), SA.end(), [&](int32_t p) { return cmp_lt(p, true); }) - SA.begin();
        return {lo, hi};
    }
    bool occurs(const char *X, u64 len) const { interval I = find(X, len); return I.sp < I.ep; }
    // all MEMs of R[0..n) (over A/C/G/T): (a, e)
    vector<std::pair<u64, u64>> mems(const char *R, u64 n) const {
        vector<std::pair<u64, u64>> out;
        u64 e = 0, prev = 0;
        for (u64 a = 0; a < n; ++a) {
            if (e < a) e = a;
            while (e < n && occurs(R + a, e + 1 - a)) ++e;
            if (e > a && (a == 0 || prev < e)) out.push_back({a, e});
            prev = e;
        }
        return out;
    }
};

struct config {
    string name;
    const rlcsa *csa = nullptr;
    const lookup_table *T = nullptr;
    const two_level *TL = nullptr;
    bool jumps = true;
    const stpd_index *SP = nullptr;
    const parse_stpd *PS = nullptr;
};

static vector<mem> run_find(const stepper &base, const config &c, const string &read, const mem_options &opt) {
    stepper S = base;
    if (c.csa) S.csa = c.csa;
    segment_search Q; Q.S = &S; Q.T = c.T; Q.TL = c.TL; Q.jumps = c.jumps; Q.SP = c.SP; Q.PS = c.PS;
    vector<mem> out;
    find_mems(Q, read, opt, [&](const mem &M) { out.push_back(M); });
    return out;
}

static bool same(const vector<mem> &x, const vector<mem> &y) {
    if (x.size() != y.size()) return false;
    for (size_t i = 0; i < x.size(); ++i) if (x[i].a != y[i].a || x[i].e != y[i].e || !(x[i].I == y[i].I)) return false;
    return true;
}

static std::set<std::pair<u64, u64>> as_set(const vector<mem> &v) { std::set<std::pair<u64, u64>> s; for (auto &M : v) s.insert({M.a, M.e}); return s; }

static int mems_test(int argc, char **argv) {
    if (argc < 4) die("usage: katka-test mems PREFIX reads [--tags FILE] [--kebab k] [--twolevel] [--hpc] [--random n]");
    const string prefix = argv[2], readfile = argv[3];
    string tagfile; uint32_t kk = 0; bool tl = false, hpc = false, digest = false, stpd = false; int nrand = 0;
    for (int i = 4; i < argc; ++i) {
        string a = argv[i];
        if (a == "--tags") tagfile = argv[++i];
        else if (a == "--kebab") kk = std::stoul(argv[++i]);
        else if (a == "--twolevel") tl = true;
        else if (a == "--hpc") hpc = true;
        else if (a == "--random") nrand = std::stoi(argv[++i]);
        else if (a == "--digest") digest = true;
        else if (a == "--stpd") stpd = true;
        else die("unknown option " + a);
    }
    brute BF; BF.build(prefix);
    rlbwt B;
    { FILE *f = open_or_die(prefix + ".rlbwt", "rb"); rmagic(f, "KTRLBW01", prefix + ".rlbwt"); B.load(f); fclose(f); }
    // the RLBWT against the suffix array
    CHECK(B.n == BF.SA.size(), "BWT length");
    {
        u64 cnt[256] = {0};
        for (u64 i = 0; i < B.n; ++i) {
            unsigned char c = BF.SA[i] == 0 ? 1 : (unsigned char)BF.S[BF.SA[i] - 1];
            if (i % 37 == 0 || i < 300) {
                CHECK(B.at(i) == c, "BWT[%lu]", i);
                for (unsigned char d : {(unsigned char)1, (unsigned char)'A', (unsigned char)'C', (unsigned char)'G', (unsigned char)'T', (unsigned char)'X'})
                    CHECK(B.rank(i, d) == cnt[d], "rank(%lu, %c)", i, d);
            }
            cnt[c]++;
        }
    }
    rlcsa X, Y;
    X.load(prefix + ".csa"); Y.load(prefix + ".efcsa");
    stepper base; base.bwt = &B;
    lookup_table LT;
    digest_map DM;
    if (digest) {               // a table over the digest symbols
        DM.load(prefix + ".map");
        std::vector<unsigned char> a;
        for (int c = 3; c < 256; ++c) if (c != 'X' && B.C[c + 1] > B.C[c]) a.push_back((unsigned char)c);
        LT.set_alphabet(a); LT.build(base, 2);
        nrand = 0;
    } else LT.build(base, 6);
    two_level L2;
    vector<config> cfgs = {{"rlbwt"}, {"csa", &X}, {"efcsa", &Y}, {"rlbwt+F6", nullptr, &LT}, {"csa+F6", &X, &LT}};
    if (tl) {
        L2.load(prefix + ".lvl");
        cfgs.push_back({"twolevel", &X, nullptr, &L2, false});
        cfgs.push_back({"twolevel+jumps", &X, nullptr, &L2, true});
        cfgs.push_back({"twolevel+jumps/rlbwt", nullptr, nullptr, &L2, true});
    }
    stpd_index SD, SDall, SDp;
    parse_text PT;
    if (stpd) {
        SD.load(prefix + ".stpd", prefix + ".S");
        SDall.build(prefix + ".S", prefix + ".bwt", B, false);
        // the same sample with the text read through a parse (short phrases: many boundary crossings)
        SDp.load(prefix + ".stpd", "");
        PT.build(prefix + ".S", 12, 3, prefix + ".ptxt.tmp");
        SDp.PT = &PT;
        {
            std::mt19937_64 rng(11);
            u64 bad = 0;
            for (u64 x = 0; x < SD.n; ++x) bad += PT.at(x) != SD.S[x];               // forward
            for (u64 x = SD.n; x-- > 0;) bad += PT.at(x) != SD.S[x];                 // backward
            for (int q = 0; q < 20000; ++q) {                                         // random, and lce
                u64 a = rng() % SD.n, b = rng() % SD.n, len = std::min(SD.n - std::max(a, b), (u64)(rng() % 300));
                bad += PT.at(a) != SD.S[a];
                bad += PT.lce(a, SD.S.data() + b, len) != lce_bytes(SD.S.data() + a, SD.S.data() + b, len);
            }
            CHECK(bad == 0, "parse oracle: %lu wrong reads", bad);
            printf("parse oracle (12,3): %lu phrases; every position read forward and backward, 20000 random reads and lce checked\n", PT.m);
        }
        cfgs.push_back({"stpd", &X, nullptr, nullptr, true, &SD});
        cfgs.push_back({"stpd/parse-oracle", &X, nullptr, nullptr, true, &SDp});
        cfgs.push_back({"stpd-all/rlbwt", nullptr, nullptr, nullptr, true, &SDall});
        if (tl) cfgs.push_back({"stpd+twolevel", &X, nullptr, &L2, true, &SD});
    }
    // STPD over the parse: the two-level index's parse as the oracle (and its own character-level copy of
    // the sample), on the RLCSA and on the RLBWT at the right end
    stpd_index SDq;
    parse_text PTq;
    parse_stpd PSq;
    if (stpd && tl) {
        SDq.load(prefix + ".stpd", "");
        PTq.build(prefix + ".S", L2.K, L2.SM, prefix + ".ptxt.tmp");
        SDq.PT = &PTq;
        PSq.build(PTq);
        cfgs.push_back({"stpd-parse", &X, nullptr, &L2, true, &SDq, &PSq});
        cfgs.push_back({"stpd-parse/rlbwt", nullptr, nullptr, &L2, true, &SDq, &PSq});
    }
    tag_array T;
    vector<uint32_t> rowtag;
    if (!tagfile.empty()) {
        T.load(tagfile, true); T.bwt = &B;
        // the tag of each row, from the genome table (by genome index, or by tag column if the array has fewer tags)
        genome_table G; G.load(prefix + ".tbl");
        u64 ntag;
        vector<uint32_t> byidx = G.tags(false), bycol = G.tags(true, &ntag);
        u64 mx = 0;
        for (u64 r = 0; r < T.rho; ++r) mx = std::max<u64>(mx, T.tag(r));
        const vector<uint32_t> &gt = mx < ntag && ntag < G.start.size() ? bycol : byidx;
        rowtag.resize(B.n);
        for (u64 i = 0; i < B.n; ++i) {
            u64 p = BF.SA[i] == (int32_t)(BF.S.size() - 1) ? BF.S.size() - 2 : BF.SA[i];
            u64 g = std::upper_bound(G.start.begin(), G.start.end(), p) - G.start.begin() - 1;
            rowtag[i] = gt[g];
        }
        for (u64 i = 0; i < B.n; i += 1 + i % 13) CHECK(T.tag(T.run_of(i)) == rowtag[i], "tag of row %lu", i);
    }
    kmer_filter KF;
    if (kk && !KF.load(prefix + ".k" + std::to_string(kk) + ".kbf")) die("cannot load the KeBaB filter");

    // patterns: the reads, then random ones
    vector<string> pats;
    {
        FILE *f = open_or_die(readfile, "rb"); fclose(f);
        std::ifstream in(readfile);
        string line, cur;
        bool fasta = in.peek() == '>';
        u64 ln = 0;
        while (std::getline(in, line)) {
            if (fasta) { if (!line.empty() && line[0] == '>') { if (!cur.empty()) pats.push_back(cur); cur.clear(); } else cur += line; }
            else if (ln++ % 4 == 1) pats.push_back(line);
        }
        if (!cur.empty()) pats.push_back(cur);
    }
    std::mt19937_64 rng(11);
    for (int i = 0; i < nrand; ++i) {
        u64 len = 50 + rng() % 400, a = rng() % (BF.S.size() - len - 1);
        string p = BF.S.substr(a, len);
        for (auto &c : p) { if (rng() % 60 == 0) c = "ACGT"[rng() % 4]; }
        if (rng() % 2) p = reverse_complement(p);
        pats.push_back(p);
    }
    for (auto &p : pats) {
        for (auto &c : p) c = (char)toupper((unsigned char)c);
        if (hpc) p.erase(std::unique(p.begin(), p.end()), p.end());
    }

    u64 nmems = 0, nchecked = 0, nlongest = 0, ntop = 0, ntop_differ = 0, kb_multi = 0, kb_fall = 0;
    u64 kc_runs = 0, kc_fall = 0, nfixed = 0;
    for (size_t pi = 0; pi < pats.size(); ++pi) {
        const string &P = pats[pi];
        if (digest) {
            // brute force (as KATKA's RZ_BML_CHECK): every digested MEM dg[a..b] of each run, with its span,
            // the spans of at least L bases, and of those the ones not inside another's span
            for (u64 L : {11ul, 15ul, 20ul, 30ul}) {
                std::set<std::pair<u64, u64>> want;
                for (u64 i = 0; i < P.size();) {
                    while (i < P.size() && base_code(P[i]) < 0) ++i;
                    u64 r0 = i;
                    while (i < P.size() && base_code(P[i]) >= 0) ++i;
                    u64 n = i - r0;
                    if (n < L || n < (u64)DM.w) continue;
                    read_digest D; D.build(DM, P.data() + r0, n);
                    const u64 nd = D.dg.size();
                    std::vector<u64> A(nd);
                    for (u64 b = 0; b < nd; ++b) {           // leftmost start of a match ending at b
                        u64 j = b + 1;
                        while (j > 0 && BF.occurs(D.dg.data() + j - 1, b + 2 - j)) --j;
                        A[b] = j;
                    }
                    std::vector<std::pair<u64, u64>> all;
                    for (u64 b = 0; b < nd; ++b) {
                        if (A[b] > b || (b + 1 < nd && A[b + 1] <= A[b])) continue;
                        u64 left, right, yR;
                        D.span(A[b], b, n, DM.w, left, right, yR);
                        if (right >= left && right - left + 1 >= L) all.push_back({left, right});
                    }
                    for (auto &p : all) {
                        bool in = false;
                        for (auto &q : all) if (q != p && q.first <= p.first && p.second <= q.second) { in = true; break; }
                        if (!in) want.insert({r0 + p.first, r0 + p.second + 1});
                    }
                }
                mem_options o; o.L = L; o.digest = &DM;
                vector<mem> ref;
                for (size_t c = 0; c < cfgs.size(); ++c) {
                    stepper S2 = base;
                    if (cfgs[c].csa) S2.csa = cfgs[c].csa;
                    segment_search Q; Q.S = &S2; Q.T = cfgs[c].T; Q.comp = DM.comp;
                    vector<mem> got;
                    find_mems(Q, P, o, [&](const mem &M) { got.push_back(M); });
                    if (c == 0) {
                        ref = got;
                        bool okd = as_set(got) == want && got.size() == want.size();
                        CHECK(okd, "digest L=%lu pattern %zu: %zu MEMs, brute %zu", L, pi, got.size(), want.size());
                        for (auto &M : got) CHECK(M.I.sp < M.I.ep, "digest: empty interval");
                        nmems += got.size(); ++nchecked;
                    } else CHECK(same(got, ref), "digest L=%lu pattern %zu: %s differs", L, pi, cfgs[c].name.c_str());
                }
            }
            continue;
        }
        // brute force: the MEMs of every A/C/G/T run
        vector<std::pair<u64, u64>> all;
        for (u64 i = 0; i < P.size();) {
            while (i < P.size() && base_code(P[i]) < 0) ++i;
            u64 r0 = i;
            while (i < P.size() && base_code(P[i]) >= 0) ++i;
            for (auto [a, e] : BF.mems(P.data() + r0, i - r0)) all.push_back({r0 + a, r0 + e});
        }
        auto check_intervals = [&](const vector<mem> &got, const char *what) {
            for (auto &M : got) {
                string rc = reverse_complement(P.substr(M.a, M.e - M.a));
                interval I = BF.find(rc.data(), rc.size());
                CHECK(I == M.I, "%s: pattern %zu, MEM [%lu,%lu): interval [%lu,%lu) != brute [%lu,%lu)", what, pi, M.a, M.e, M.I.sp, M.I.ep, I.sp, I.ep);
            }
        };
        for (u64 L : {1ul, 12ul, 25ul, 40ul, 100ul, 1000ul}) {
            std::set<std::pair<u64, u64>> want;
            for (auto [a, e] : all) if (e - a >= L) want.insert({a, e});
            for (int fb = 0; fb < 2; ++fb) {
                mem_options o; o.L = L; o.fb = fb;
                vector<mem> ref = run_find(base, cfgs[0], P, o);
                CHECK(as_set(ref) == want && ref.size() == want.size(), "%s L=%lu pattern %zu: %zu MEMs, brute %zu", fb ? "FB" : "BML", L, pi, ref.size(), want.size());
                check_intervals(ref, fb ? "FB" : "BML");
                nmems += ref.size(); ++nchecked;
                for (size_t c = 1; c < cfgs.size(); ++c)
                    CHECK(same(run_find(base, cfgs[c], P, o), ref), "%s L=%lu pattern %zu: %s differs from plain search", fb ? "FB" : "BML", L, pi, cfgs[c].name.c_str());
            }
            if (kk && L >= kk) {
                for (int mode = 0; mode < 2; ++mode) {
                    mem_options o; o.L = L; o.filter = &KF; (mode ? o.inline_check : o.pseudo) = true;
                    for (auto &cf : cfgs) {
                        vector<mem> got = run_find(base, cf, P, o);
                        CHECK(as_set(got) == want && got.size() == want.size(), "KeBaB %s L=%lu pattern %zu (%s): %zu MEMs, brute %zu", mode ? "inline" : "pseudo", L, pi, cf.name.c_str(), got.size(), want.size());
                    }
                }
            }
        }
        // adaptive L and top t
        u64 maxlen = 0;
        for (auto [a, e] : all) maxlen = std::max(maxlen, e - a);
        vector<u64> lens;
        for (auto [a, e] : all) lens.push_back(e - a);
        std::sort(lens.rbegin(), lens.rend());
        std::set<std::pair<u64, u64>> allset(all.begin(), all.end());
        for (int kb = 0; kb < (kk ? 2 : 1); ++kb) {
            mem_options o; o.L = 1; o.adaptive = true;
            if (kb) { o.filter = &KF; o.pseudo = true; o.L = kk; }
            vector<mem> got = run_find(base, cfgs[0], P, o);
            std::set<std::pair<u64, u64>> gs = as_set(got);
            for (auto &M : got) CHECK(allset.count({M.a, M.e}), "adaptive: [%lu,%lu) is not a MEM", M.a, M.e);
            if (maxlen >= o.L) for (auto [a, e] : all) if (e - a == maxlen) CHECK(gs.count({a, e}), "adaptive%s: longest MEM [%lu,%lu) missing (pattern %zu)", kb ? "+KeBaB" : "", a, e, pi);
            check_intervals(got, "adaptive");
            for (size_t c = 1; c < cfgs.size(); ++c) CHECK(same(run_find(base, cfgs[c], P, o), got), "adaptive: %s differs", cfgs[c].name.c_str());
            for (u64 t : {1ul, 3ul, 10ul}) {
                mem_options q = o; q.adaptive = false; q.top = t;
                vector<mem> top = run_find(base, cfgs[0], P, q);
                vector<u64> tl;
                for (auto &M : top) { tl.push_back(M.e - M.a); CHECK(allset.count({M.a, M.e}), "top: [%lu,%lu) is not a MEM", M.a, M.e); }
                std::sort(tl.rbegin(), tl.rend());
                vector<u64> want;
                for (u64 x : lens) if (x >= q.L && want.size() < t) want.push_back(x);
                CHECK(tl == want, "top %lu%s pattern %zu: lengths differ (%zu vs %zu)", t, kb ? "+KeBaB" : "", pi, tl.size(), want.size());
                check_intervals(top, "top");
                for (size_t c = 1; c < cfgs.size(); ++c) CHECK(same(run_find(base, cfgs[c], P, q), top), "top: %s differs", cfgs[c].name.c_str());
            }
        }
        // the t longest with ties (--longest), and the KeBaB rules: streaming and batch exact, top (KeBaB's
        // heuristic) only true MEMs
        for (int min3 : {0, 14})
            for (u64 t : {1ul, 3ul, 5ul}) {
                auto want_longest = [&](u64 L) {
                    vector<std::pair<u64, u64>> ok;
                    for (auto [a, e] : all) if (e - a >= L && (!min3 || distinct_3mers(P.data() + a, e - a) >= min3)) ok.push_back({a, e});
                    vector<u64> ls;
                    for (auto [a, e] : ok) ls.push_back(e - a);
                    std::sort(ls.rbegin(), ls.rend());
                    u64 lt = ls.size() >= t ? ls[t - 1] : 0;
                    std::set<std::pair<u64, u64>> w;
                    for (auto [a, e] : ok) if (e - a >= lt) w.insert({a, e});
                    return w;
                };
                auto run_check = [&](mem_options q, const std::set<std::pair<u64, u64>> &w, const char *what) {
                    vector<mem> got = run_find(base, cfgs[0], P, q);
                    CHECK(as_set(got) == w && got.size() == w.size(), "%s t=%lu min3=%d L=%lu pattern %zu: %zu MEMs, brute %zu", what, t, min3, q.L, pi, got.size(), w.size());
                    check_intervals(got, what);
                    for (size_t c = 1; c < cfgs.size(); ++c) CHECK(same(run_find(base, cfgs[c], P, q), got), "%s: %s differs", what, cfgs[c].name.c_str());
                    ++nlongest;
                    return got;
                };
                mem_options q; q.L = 1; q.longest = t; q.min3 = min3;
                run_check(q, want_longest(1), "longest");
                if (!kk) continue;
                q.filter = &KF; q.pseudo = true;
                for (u64 L : {(u64)kk, 1ul}) {
                    q.L = L;
                    const auto w = want_longest(L);
                    if (L >= kk) { q.kebab_rule = mem_options::STREAM; run_check(q, w, "KeBaB streaming"); }
                    q.kebab_rule = mem_options::BATCH;
                    segment_search Qs;                       // (for the batch and fallback counts)
                    {
                        stepper S2 = base; Qs.S = &S2;
                        find_mems(Qs, P, q, [](const mem &) {});
                        kb_multi += Qs.cnt.kb_multibatch; kb_fall += Qs.cnt.kb_fallbacks;
                    }
                    run_check(q, w, "KeBaB batch");
                    q.kebab_rule = mem_options::TOP;
                    vector<mem> got = run_find(base, cfgs[0], P, q);
                    for (auto &M : got) CHECK(allset.count({M.a, M.e}), "KeBaB top: [%lu,%lu) is not a MEM", M.a, M.e);
                    ++ntop; ntop_differ += as_set(got) != w;
                    // --kebab-check T t: strict gives exactly the t longest with ties; non-strict the same lengths
                    // of the t longest, and only MEMs among them
                    for (u64 T : {t, 2 * t, t + 2})
                        for (int strict = 0; strict < 2; ++strict) {
                            mem_options c = q; c.kebab_rule = mem_options::CHECK; c.check_T = T; c.check_strict = strict;
                            segment_search Qc; stepper S3 = base; Qc.S = &S3;
                            vector<mem> g;
                            find_mems(Qc, P, c, [&](const mem &M) { g.push_back(M); });
                            kc_fall += Qc.cnt.kb_fallbacks; ++kc_runs;
                            std::set<std::pair<u64, u64>> gs = as_set(g);
                            if (strict) CHECK(gs == w && g.size() == w.size(), "KeBaB check T=%lu t=%lu strict L=%lu pattern %zu: %zu MEMs, brute %zu", T, t, L, pi, g.size(), w.size());
                            else {
                                bool sub = true;
                                for (auto &x : gs) sub = sub && w.count(x);
                                vector<u64> lg, lw;
                                for (auto &x : gs) lg.push_back(x.second - x.first);
                                for (auto &x : w) lw.push_back(x.second - x.first);
                                std::sort(lg.rbegin(), lg.rend()); std::sort(lw.rbegin(), lw.rend());
                                if (lg.size() > t) lg.resize(t);
                                if (lw.size() > t) lw.resize(t);
                                CHECK(sub && lg == lw, "KeBaB check T=%lu t=%lu L=%lu pattern %zu: not the t longest", T, t, L, pi);
                            }
                            for (size_t cf = 1; cf < cfgs.size(); ++cf) CHECK(same(run_find(base, cfgs[cf], P, c), run_find(base, cfgs[0], P, c)), "KeBaB check: %s differs", cfgs[cf].name.c_str());
                        }
                    // --kebab-top t --fixed-L: every MEM >= max(L, k) inside the t longest pseudo-MEMs
                    if (!min3) {
                        mem_options f; f.L = L; f.filter = &KF; f.pseudo = true; f.fixed_top = t;
                        vector<std::pair<u64, u64>> sg;
                        segments(P, &KF, sg);
                        std::stable_sort(sg.begin(), sg.end(), [](const std::pair<u64, u64> &x, const std::pair<u64, u64> &y) { return x.second > y.second; });
                        if (sg.size() > t) sg.resize(t);
                        std::set<std::pair<u64, u64>> wf;
                        const u64 Lf = std::max<u64>(L, kk);
                        for (auto [a, e] : all)
                            if (e - a >= Lf) for (auto [r0, n] : sg) if (r0 <= a && e <= r0 + n) { wf.insert({a, e}); break; }
                        vector<mem> g = run_find(base, cfgs[0], P, f);
                        CHECK(as_set(g) == wf && g.size() == wf.size(), "KeBaB top fixed L=%lu t=%lu pattern %zu: %zu MEMs, brute %zu", L, t, pi, g.size(), wf.size());
                        check_intervals(g, "KeBaB top fixed");
                        ++nfixed;
                    }
                }
            }
        // tags of the MEMs of length >= 12
        if (!tagfile.empty()) {
            mem_options o; o.L = 12;
            for (auto &M : run_find(base, cfgs[1], P, o)) {
                std::map<u64, u64> want;
                for (u64 r = M.I.sp; r < M.I.ep; ++r) want[rowtag[r]]++;
                vector<std::pair<u64, u64>> got;
                T.count(M.I.sp, M.I.ep, got);
                bool okc = got == vector<std::pair<u64, u64>>(want.begin(), want.end());
                CHECK(okc, "count: pattern %zu MEM [%lu,%lu)", pi, M.a, M.e);
                if (T.has_rmq) {
                    vector<u64> lst;
                    T.list(M.I.sp, M.I.ep, lst);
                    std::sort(lst.begin(), lst.end());
                    vector<u64> wl;
                    for (auto &x : want) wl.push_back(x.first);
                    CHECK(lst == wl, "list: pattern %zu MEM [%lu,%lu)", pi, M.a, M.e);
                }
            }
        }
        // two-level search: every longest-suffix query against plain search, on both strands
        if (tl) {
            stepper S = base; S.csa = &X;
            segment_search Q0, Q1, Q2;
            Q0.S = Q1.S = Q2.S = &S;
            Q1.TL = Q2.TL = &L2; Q1.jumps = false;
            for (u64 i = 0; i < P.size();) {
                while (i < P.size() && base_code(P[i]) < 0) ++i;
                u64 r0 = i;
                while (i < P.size() && base_code(P[i]) >= 0) ++i;
                u64 n = i - r0;
                if (!n) continue;
                Q0.prepare(P.data() + r0, n); Q1.prepare(P.data() + r0, n); Q2.prepare(P.data() + r0, n);
                for (int q = 0; q < 60; ++q) {
                    u64 end = 1 + rng() % n, lo = q % 3 ? 0 : rng() % end;
                    int dir = q & 1;
                    interval I0, I1, I2;
                    u64 j0 = Q0.back(dir, lo, end, I0), j1 = Q1.back(dir, lo, end, I1), j2 = Q2.back(dir, lo, end, I2);
                    CHECK(j0 == j1 && I0 == I1, "two-level (parse) differs: pattern %zu dir %d [%lu,%lu)", pi, dir, lo, end);
                    CHECK(j0 == j2 && I0 == I2, "two-level (blocks) differs: pattern %zu dir %d [%lu,%lu)", pi, dir, lo, end);
                }
            }
        }
    }
    if (stpd) {   // every STPD search ends at the lexicographically smallest occurrence: t = SA[sp]
        stepper S = base; S.csa = &X;
        u64 nq = 0, nloc = 0;
        for (size_t pi = 0; pi < pats.size(); ++pi) {
            const string &P = pats[pi];
            for (u64 i = 0; i < P.size();) {
                while (i < P.size() && base_code(P[i]) < 0) ++i;
                u64 r0 = i;
                while (i < P.size() && base_code(P[i]) >= 0) ++i;
                u64 n = i - r0;
                if (!n) continue;
                for (const stpd_index *DD : {&SD, &SDp}) {
                segment_search Q0, Q1; Q0.S = Q1.S = &S; Q1.SP = DD;
                Q0.prepare(P.data() + r0, n); Q1.prepare(P.data() + r0, n);
                for (int q = 0; q < 20; ++q) {
                    u64 end = 1 + rng() % n, lo = q % 3 ? 0 : rng() % end;
                    int dir = q & 1;
                    interval I0, I1;
                    u64 j0 = Q0.back(dir, lo, end, I0), j1 = Q1.back(dir, lo, end, I1);
                    bool ok = j0 == j1 && (j0 == end || (u64)BF.SA[I0.sp] == Q1.sst[dir].t);
                    CHECK(ok, "STPD: pattern %zu dir %d [%lu,%lu): j %lu vs %lu", pi, dir, lo, end, j1, j0);
                    if (ok && j0 < end) {          // every occurrence by phi^-1, against the suffix array
                        vector<u64> occ; u64 checks = 0;
                        const char *X = dir ? Q1.rc.data() : Q1.R;
                        DD->locate((const uint8_t *)X + j0, end - j0, Q1.sst[dir].t, occ, checks);
                        vector<u64> want;
                        for (u64 x = I0.sp; x < I0.ep; ++x) want.push_back(BF.SA[x]);
                        CHECK(occ == want, "STPD locate: pattern %zu dir %d [%lu,%lu): %zu occurrences, SA range %zu", pi, dir, j0, end, occ.size(), want.size());
                        nloc += occ.size();
                    }
                    ++nq;
                }
                }
            }
        }
        if (tl) {      // STPD over the parse: longest suffixes and every occurrence (phi^-1 over the parse)
            u64 pq = 0, ploc = 0;
            for (size_t pi = 0; pi < pats.size(); ++pi) {
                const string &P = pats[pi];
                for (u64 i = 0; i < P.size();) {
                    while (i < P.size() && base_code(P[i]) < 0) ++i;
                    u64 r0 = i;
                    while (i < P.size() && base_code(P[i]) >= 0) ++i;
                    u64 n = i - r0;
                    if (!n) continue;
                    segment_search Q0, Q1; Q0.S = Q1.S = &S; Q1.SP = &SDq; Q1.PS = &PSq; Q1.TL = &L2; Q1.locate = true;
                    Q0.prepare(P.data() + r0, n); Q1.prepare(P.data() + r0, n);
                    for (int q = 0; q < 30; ++q) {
                        u64 end = 1 + rng() % n, lo = q % 3 ? 0 : rng() % end;
                        int dir = q & 1;
                        interval I0, I1;
                        u64 j0 = Q0.back(dir, lo, end, I0), j1 = Q1.back(dir, lo, end, I1);
                        CHECK(j0 == j1, "STPD over the parse: pattern %zu dir %d [%lu,%lu): j %lu vs %lu", pi, dir, lo, end, j1, j0);
                        if (dir == 1 && j0 == j1 && j0 < end) {
                            const u64 a = n - end, e = n - j0;
                            Q1.ps_locate(a, e);
                            vector<u64> got = Q1.occbuf, want;
                            for (u64 x = I0.sp; x < I0.ep; ++x) want.push_back(BF.SA[x]);
                            std::sort(got.begin(), got.end()); std::sort(want.begin(), want.end());
                            CHECK(got == want, "STPD over the parse, locate: pattern %zu [%lu,%lu): %zu occurrences, SA %zu", pi, j0, end, got.size(), want.size());
                            ploc += got.size();
                        }
                        ++pq;
                    }
                    // extend_left resumption: back to x, then extend to 0 (as BML does)
                    for (int q = 0; q < 10; ++q) {
                        u64 end = 1 + rng() % n, x = rng() % end;
                        interval I0, I1;
                        u64 j0 = Q0.back(0, x, end, I0), j1 = Q1.back(0, x, end, I1);
                        CHECK(j0 == j1, "STPD over the parse (BML window): pattern %zu", pi);
                        if (j0 == x && j1 == x) {
                            u64 a0 = Q0.extend_left(0, 0, x, end, I0), a1 = Q1.extend_left(0, 0, x, end, I1);
                            CHECK(a0 == a1, "STPD over the parse, extend_left: pattern %zu [%lu,%lu): %lu vs %lu", pi, x, end, a1, a0);
                        }
                        ++pq;
                    }
                }
            }
            printf("STPD over the parse: %lu longest-suffix searches checked against plain search, %lu occurrences located by phi^-1 over the parse checked against SA\n", pq, ploc);
        }
        printf("STPD (text and parse oracle): %lu longest-suffix searches checked (start and lex-smallest occurrence), %lu occurrences located by phi^-1 checked against SA\n", nq, nloc);
    }
    printf("longest t (t = 1, 3, 5; -D 0 and 14): %lu exact runs checked (no filter, KeBaB streaming, KeBaB batch: %lu second batches, %lu whole-read fallbacks); KeBaB top: %lu of %lu runs differ from the t longest\n",
           nlongest, kb_multi, kb_fall, ntop_differ, ntop);
    printf("KeBaB check (T = t, 2t, t + 2; strict and not): %lu runs checked, %lu fell back; KeBaB top with a fixed L: %lu runs checked\n", kc_runs, kc_fall, nfixed);
    printf("MEM tests on %zu patterns (%lu runs, %lu MEMs from plain BML/FB): %s\n", pats.size(), nchecked, nmems, failures ? "FAILED" : "passed");
    return failures ? 1 : 0;
}

// |st-colex-| of T = rev(S) $ by the paper's definition (Definition 7 with pi = IPA): the number of
// distinct i + LPF[i], where LPF[i] = max { lcp(T[i..], T[j..]) : pi(j) < pi(i) } and pi(i) = the colex rank
// of T[0..i), i.e. the lex rank of the suffix S[n-i..] of S$.  LPF[i] is the larger lcp with the nearest
// suffixes before and after T[i..] in SA order with smaller pi.  Compared with the size of PREFIX.stpd,
// which is built on S from the BWT runs (the mirror image).
static int stpd_size(const string &prefix) {
    brute BF; BF.build(prefix);                 // SA of S$ (BF.S ends with byte 1)
    const u64 n = BF.S.size() - 1, N = n + 1;
    vector<int32_t> isaS(N);
    for (u64 x = 0; x < N; ++x) isaS[BF.SA[x]] = (int32_t)x;
    string T(BF.S.rbegin() + 1, BF.S.rend());   // rev(S)
    vector<int32_t> sa(n), SA(N), rk(N), lcp(N, 0);
    libsais((const uint8_t *)T.data(), sa.data(), (int32_t)n, 0, nullptr);
    T.push_back(1);
    SA[0] = (int32_t)n; for (u64 i = 0; i < n; ++i) SA[i + 1] = sa[i];
    for (u64 x = 0; x < N; ++x) rk[SA[x]] = (int32_t)x;
    for (u64 i = 0, h = 0; i < N; ++i) {        // Kasai
        if (rk[i] == 0) { h = 0; continue; }
        u64 j = SA[rk[i] - 1];
        while (i + h < N && j + h < N && T[i + h] == T[j + h]) ++h;
        lcp[rk[i]] = (int32_t)h;
        if (h) --h;
    }
    auto pi = [&](u64 i) -> int64_t { return isaS[n - i]; };   // T[0..i) reversed is S[n-i..n)
    vector<int64_t> best(N, 0);
    for (int dir = 0; dir < 2; ++dir) {          // nearest smaller pi before / after, in SA order
        // stack of (pi, h): h = lcp with the entry above it, or, for the top, with the last position seen
        vector<std::pair<int64_t, int64_t>> stk;
        for (u64 q = 0; q < N; ++q) {
            u64 x = dir ? N - 1 - q : q;
            int64_t p = pi(SA[x]);
            if (!stk.empty()) {
                int64_t g = dir ? lcp[x + 1] : lcp[x];          // lcp with the previous position seen
                stk.back().second = std::min(stk.back().second, g);
                while (!stk.empty() && stk.back().first > p) {
                    int64_t h = stk.back().second;
                    stk.pop_back();
                    if (!stk.empty()) stk.back().second = std::min(stk.back().second, h);
                }
                if (!stk.empty()) best[SA[x]] = std::max(best[SA[x]], stk.back().second);
            }
            stk.push_back({p, INT64_MAX});
        }
    }
    std::set<u64> pda;
    for (u64 i = 0; i < N; ++i) pda.insert(i + best[i]);
    FILE *f = open_or_die(prefix + ".stpd", "rb"); rmagic(f, "KTSTPD02", prefix + ".stpd");
    u64 nn, r; rd(f, nn); rd(f, r); PackedInts pos; pos.load(f); fclose(f);
    printf("|st-colex-(rev(S)$)| by definition = %zu; PREFIX.stpd sample = %lu; r = %lu\n", pda.size(), pos.size(), r);
    return 0;
}

// PREFIX.stpd, built on S = rev(T), against the st-colex- sampling of T from the paper's code
// (STPD-index's stpd-sampling -c: 5-byte little-endian positions, 0-based ends of the sampled prefixes):
// our position b in S is their n - b - 1.  Their sampling also has 0 (the first prefix), which a search
// never needs.
static int stpd_vs(const string &prefix, const string &colex) {
    FILE *f = open_or_die(colex, "rb");
    std::set<u64> A;
    unsigned char b[5];
    while (fread(b, 1, 5, f) == 5) { u64 x = 0; for (int i = 4; i >= 0; --i) x = x << 8 | b[i]; A.insert(x); }
    fclose(f);
    stpd_index D; D.load(prefix + ".stpd", prefix + ".S");
    const u64 n = D.n - 1;
    std::set<u64> B;
    for (u64 i = 0; i < D.pos.size(); ++i) B.insert(n - D.pos[i] - 1);
    u64 common = 0;
    for (u64 x : B) common += A.count(x);
    std::set<u64> onlyA;
    for (u64 x : A) if (!B.count(x)) onlyA.insert(x);
    bool ok = common == B.size() && (onlyA.empty() || (onlyA.size() == 1 && *onlyA.begin() == 0));
    printf("STPD-index st-colex-: %zu positions; ours: %zu; common %lu; only theirs: %zu%s -> %s\n", A.size(), B.size(), common, onlyA.size(),
           onlyA.size() == 1 && *onlyA.begin() == 0 ? " (position 0)" : "", ok ? "the same" : "DIFFERENT");
    return ok ? 0 : 1;
}

// B_lambda by brute force (suffix array and LCP of S$): tag boundaries p (tag of row p-1 != tag of row p)
// with lcp(p) > lambda, tags by genome (bycol 0) or by tag column (1).  Prints the same lines as
// katka index --tag-lcp.
static int blam(const string &prefix, bool bycol, const string &ls) {
    brute BF; BF.build(prefix);
    genome_table G; G.load(prefix + ".tbl");
    vector<uint32_t> gt = G.tags(bycol);
    const u64 N = BF.S.size();
    vector<int32_t> rk(N);
    for (u64 x = 0; x < N; ++x) rk[BF.SA[x]] = (int32_t)x;
    vector<u64> lcp(N, 0);
    for (u64 i = 0, h = 0; i < N; ++i) {
        if (rk[i] == 0) { h = 0; continue; }
        u64 j = BF.SA[rk[i] - 1];
        while (i + h < N && j + h < N && BF.S[i + h] == BF.S[j + h]) ++h;
        lcp[rk[i]] = h;
        if (h) --h;
    }
    auto tag = [&](u64 x) { u64 p = BF.SA[x] == (int32_t)(N - 1) ? N - 2 : BF.SA[x];
        return gt[std::upper_bound(G.start.begin(), G.start.end(), p) - G.start.begin() - 1]; };
    vector<u64> lambdas; { string w; for (char c : ls + ",") { if (c == ',') { if (!w.empty()) lambdas.push_back(std::stoull(w)); w.clear(); } else w += c; } }
    u64 rho = 1;
    vector<u64> B(lambdas.size(), 0);
    for (u64 x = 1; x < N; ++x) if (tag(x) != tag(x - 1)) { ++rho; for (size_t i = 0; i < lambdas.size(); ++i) if (lcp[x] > lambdas[i]) ++B[i]; }
    printf("tag runs rho = %lu\n", rho);
    for (size_t i = 0; i < lambdas.size(); ++i) printf("  lambda = %lu: B_lambda = %lu\n", lambdas[i], B[i]);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && string(argv[1]) == "unit") return unit();
    if (argc == 5 && string(argv[1]) == "blam") return blam(argv[2], string(argv[3]) == "1", argv[4]);
    if (argc == 4 && string(argv[1]) == "stpd-vs") return stpd_vs(argv[2], argv[3]);
    if (argc == 3 && string(argv[1]) == "stpd-size") return stpd_size(argv[2]);
    if (argc >= 2 && string(argv[1]) == "mems") return mems_test(argc, argv);
    fputs("usage: katka-test unit | katka-test mems PREFIX reads [options]\n", stderr);
    return 1;
}
