// katka: one index and one search engine for the MEM-finding of "Three Examples of MEM-Finding".
//
//   katka prep  [-d] [-s] [-U] list.tsv PREFIX       the text S of both strands (prep.hpp)
//   katka bwt   PREFIX                               BWT(S$) with libsais (or bring your own: Big-BWT, pfp-merge)
//   katka index [options] PREFIX                     the RLBWT and whatever else is asked for
//   katka find  [options] PREFIX reads out           MEMs of each read, answered from the tags
//   katka info  PREFIX                               sizes of the index files
// Run a command with -h for its options; README.md lists them all, with the presets.
#include "mems.hpp"
#include "prep.hpp"
#include "tag_build.hpp"
#include "twolevel_build.hpp"
#include "stpd_export.hpp"
#include "parse_text.hpp"
#include "pstpd.hpp"
#include "tloc_build.hpp"
#include "libsais.h"
#include "libsais64.h"
#include <atomic>
#include <condition_variable>
#include <getopt.h>
#include <sys/resource.h>
#include <map>
#include <mutex>
#include <thread>

using namespace katka;
using std::string;

// ---------------------------------------------------------------------------------------------------
// presets

struct preset {
    const char *name, *use, *prep, *index, *find;
};
static const preset PRESETS[] = {
    {"zymo", "noisy long reads, few taxa",
     "-d -s -U", "--csa --tags --by-tag --no-rmq --kebab 20",
     "-U -X -L 20 -K 20 -D 14 --count --csa"},
    {"silva", "short reads, many taxa",
     "", "--csa --tags",
     "-L 30 --count --csa -F 10"},
    {"chr19", "long accurate matches",
     "", "--csa --twolevel",
     "-L 1000 --twolevel --format mems --csa"},
    {"novel", "a few long MEMs in mostly novel sequence",
     "", "--csa --kebab 20",
     "-L 1000 --inline-kebab 20 --format mems --csa"},
};

static const preset *find_preset(const string &name) {
    for (auto &p : PRESETS) if (name == p.name) return &p;
    die("unknown preset " + name + " (zymo, silva, chr19, novel)");
}

// argv with a preset's options inserted after the command (so that later options override them)
static std::vector<string> expand_preset(int argc, char **argv, const char *preset::*field) {
    std::vector<string> a{argv[0]};
    string pname;
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        if (s == "--preset" && i + 1 < argc) { pname = argv[++i]; continue; }
        if (s.rfind("--preset=", 0) == 0) { pname = s.substr(9); continue; }
    }
    if (!pname.empty()) {
        const preset *p = find_preset(pname);
        string opts = p->*field, w;
        for (char c : opts + " ") { if (c == ' ') { if (!w.empty()) a.push_back(w); w.clear(); } else w += c; }
        fprintf(stderr, "preset %s (%s): %s\n", p->name, p->use, (p->*field)[0] ? p->*field : "(no options)");
    }
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        if (s == "--preset") { ++i; continue; }
        if (s.rfind("--preset=", 0) == 0) continue;
        a.push_back(s);
    }
    return a;
}
struct argv_holder {
    std::vector<string> s; std::vector<char *> p;
    explicit argv_holder(std::vector<string> v) : s(std::move(v)) { for (auto &x : s) p.push_back(&x[0]); p.push_back(nullptr); }
    int argc() const { return (int)s.size(); }
    char **argv() { return p.data(); }
};

// ---------------------------------------------------------------------------------------------------
// katka prep

static int cmd_prep(int argc, char **argv) {
    argv_holder A(expand_preset(argc, argv, &preset::prep));
    prep_options o;
    int c;
    optind = 1;
    const char *usage = "usage: katka prep [--preset NAME] [-d] [-s] [-U] [--digest k,w] list.tsv PREFIX\n";
    static option longopts[] = {{"digest", 1, 0, 'B'}, {0, 0, 0, 0}};
    while ((c = getopt_long(A.argc(), A.argv(), "dsUhB:", longopts, nullptr)) != -1) {
        if (c == 'd') o.datasets = true; else if (c == 's') o.split = true; else if (c == 'U') o.hpc = true;
        else if (c == 'B') { if (sscanf(optarg, "%d,%d", &o.dk, &o.dw) != 2) die("--digest k,w"); }
        else { fputs(usage, stderr); return 1; }
    }
    if (A.argc() - optind != 2) { fputs(usage, stderr); return 1; }
    prep(A.argv()[optind], A.argv()[optind + 1], o);
    return 0;
}

// ---------------------------------------------------------------------------------------------------
// katka bwt: BWT(S$) with libsais, one byte per row, the terminator as byte 0

// With -d, the BWT of the datasets as separate strings, one terminator each (what pfp-merge builds): the
// suffix array of T_0 $_0 T_1 $_1 ..., with distinct terminators smaller than every base, gives the order
// of the rotations of the strings; each row's character is the one before it in its own string, cyclically.
static int bwt_datasets(const string &prefix) {
    genome_table G; G.load(prefix + ".tbl");
    std::vector<u64> ds = G.dataset_starts();
    const u64 N = file_size(prefix + ".S"), K = ds.size();
    if (K > 60) die("katka bwt -d: at most 60 datasets (use pfp-merge)");
    if (N + K >= (1ull << 31) - 16) die("katka bwt -d: text too long for this builder (use pfp-merge)");
    std::vector<uint8_t> T(N + K);
    std::vector<u64> where(N + K);            // the dataset of each position
    {
        FILE *f = open_or_die(prefix + ".S", "rb");
        std::vector<uint8_t> S(N);
        if (fread(S.data(), 1, N, f) != N) die("cannot read S");
        fclose(f);
        ds.push_back(N);
        for (u64 i = 0; i < N; ++i) if (S[i] <= K) die("katka bwt -d: S has symbols below " + std::to_string(K + 1) + " (a digest?); use pfp-merge");
        u64 p = 0;
        for (u64 j = 0; j < K; ++j) {
            for (u64 i = ds[j]; i < ds[j + 1]; ++i) { where[p] = j; T[p++] = S[i]; }
            where[p] = j; T[p++] = (uint8_t)(1 + j);
        }
    }
    std::vector<u64> dstart(K + 1, 0);        // the start of each string in T
    for (u64 j = 0; j < K; ++j) dstart[j + 1] = dstart[j] + (ds[j + 1] - ds[j]) + 1;
    std::vector<int32_t> SA(N + K);
    if (libsais(T.data(), SA.data(), (int32_t)(N + K), 0, nullptr) != 0) die("libsais failed");
    FILE *o = open_or_die(prefix + ".bwt", "wb");
    for (u64 r = 0; r < N + K; ++r) {
        u64 p = SA[r], j = where[p];
        uint8_t c = p == dstart[j] ? T[dstart[j + 1] - 1] : T[p - 1];
        fputc(c <= 60 ? 0 : c, o);
    }
    fclose(o);
    fprintf(stderr, "BWT of %lu datasets (%lu characters with terminators)\n", K, N + K);
    return 0;
}

static int cmd_bwt(int argc, char **argv) {
    if (argc == 3 && string(argv[1]) == "-d") return bwt_datasets(argv[2]);
    if (argc != 2) { fputs("usage: katka bwt [-d] PREFIX     (reads PREFIX.S, writes PREFIX.bwt; -d: one string per dataset)\n", stderr); return 1; }
    const string prefix = argv[1];
    double t0 = now();
    const u64 n = file_size(prefix + ".S");
    std::vector<uint8_t> T(n);
    { FILE *f = open_or_die(prefix + ".S", "rb"); if (fread(T.data(), 1, n, f) != n) die("cannot read " + prefix + ".S"); fclose(f); }
    for (u64 i = 0; i < n; ++i) if (T[i] <= 2) die("S contains a byte 0, 1 or 2");
    std::vector<uint8_t> U(n);
    int64_t pidx;
    if (n < (1ull << 31) - 16) {
        std::vector<int32_t> A(n + 1);
        pidx = libsais_bwt(T.data(), U.data(), A.data(), (int32_t)n, 0, nullptr);
    } else {
        std::vector<int64_t> A(n + 1);
        pidx = libsais64_bwt(T.data(), U.data(), A.data(), (int64_t)n, 0, nullptr);
    }
    if (pidx < 0) die("libsais failed");
    std::vector<uint8_t>().swap(T);
    // BWT(S$) = U[0 .. p) $ U[p .. n), where row 0 (the suffix $) has BWT character S[n-1] = U[0]
    FILE *o = open_or_die(prefix + ".bwt", "wb");
    fwrite(U.data(), 1, pidx, o); fputc(0, o); fwrite(U.data() + pidx, 1, n - pidx, o);
    fclose(o);
    fprintf(stderr, "BWT of %lu + 1 characters in %.1f s\n", n, now() - t0);
    return 0;
}

// peak resident memory of this process so far, in GB (10^9 bytes)
// Peak memory: VmHWM (which reset_step_peak resets, so each index step can report its own peak), with the
// peaks of earlier steps kept in g_peak; getrusage where /proc is missing.
static double g_peak = 0;
static double vmhwm_gb() {
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return -1;
    char line[256]; double kb = -1;
    while (fgets(line, sizeof line, f)) if (!strncmp(line, "VmHWM:", 6)) { kb = atof(line + 6); break; }
    fclose(f);
    return kb < 0 ? -1 : kb * 1024 / 1e9;
}
static double step_peak_gb() {
    double v = vmhwm_gb();
    if (v >= 0) return v;
    struct rusage u; getrusage(RUSAGE_SELF, &u); return u.ru_maxrss * 1024.0 / 1e9;
}
static double peak_gb() { return std::max(g_peak, step_peak_gb()); }
static void reset_step_peak() {
    g_peak = peak_gb();
    if (FILE *f = fopen("/proc/self/clear_refs", "w")) { fputs("5", f); fclose(f); }
}

// ---------------------------------------------------------------------------------------------------
// katka index

static const char *INDEX_USAGE =
    "usage: katka index [--preset NAME] [options] PREFIX\n"
    "  reads PREFIX.S, PREFIX.tbl and the BWT (PREFIX.bwt, or --bwt FILE); always writes PREFIX.rlbwt\n"
    "  --bwt FILE          the BWT of S (one byte per row; terminators 0, 1 or 2)\n"
    "  --csa               explicit RLCSA (PREFIX.csa)\n"
    "  --efcsa             Elias-Fano RLCSA (PREFIX.efcsa)\n"
    "  --tags              tag array (PREFIX.tag, PREFIX.tag.rmq)\n"
    "    --by-tag          tags are the tag column of PREFIX.tbl (numbered by first appearance), not genomes\n"
    "    --sample s        store one tag in about s runs, the others by LF (default 1: all)\n"
    "    --no-rmq          no RMQ: counting only, no listing\n"
    "    --grammar         grammar-compress the runs' tags (built-in approximate RePair)\n"
    "    --bigrepair P     grammar-compress them with bigrepair -i output P.R, P.C\n"
    "    --tag-out FILE    write the tag array to FILE instead of PREFIX.tag\n"
    "    --tag-lcp l1,l2,..  also report B_l, the tag boundaries whose suffixes share more than l characters,\n"
    "                      and the bound (r + sigma + 2d)(l + 1) + B_l + 1 on the number of tag runs\n"
    "    --blam-mem G      memory budget for B_lambda (default 6 GB): a text within it is mapped and compared\n"
    "                      at each boundary; a larger one goes through boundary records (on disk) and\n"
    "                      fingerprints, 8 bytes per tag boundary per lambda, lambdas grouped to fit\n"
    "    --blam-external   the boundary records and fingerprints even for a text within the budget\n"
    "  --tmp DIR           large temporary files (the phi samples, the tag runs, B_lambda's records) in DIR\n"
    "  --kebab k           KeBaB Bloom filter of the k-mers of S (PREFIX.k<k>.kbf)\n"
    "    --kebab-eps e     its false-positive rate (default 0.1)\n"
    "  --stpd              suffix-tree path decomposition sample (PREFIX.stpd; needs one terminator)\n"
    "    --stpd-all        keep every BWT run start (no pruning)\n"
    "  --twolevel          two-level index with blocks (PREFIX.lvl; needs one terminator)\n"
    "    --syncmer k,s     closed syncmers (default 51,12)\n"
    "    --no-blocks       parse only\n"
    "    --continuation smallest|frequent   successor that continues a block (default smallest)\n"
    "    --max-block B     a block also ends after B phrases (either rule; default no cap)\n"
    "  --stpd-parse        STPD over the parse (PREFIX.pstpd: parse sample, entry table, C, parse phi^-1), from\n"
    "                      PREFIX.ptxt (built first if missing or not of PREFIX.lvl's k, s); needs PREFIX.lvl;\n"
    "                      find --stpd-parse also needs PREFIX.stpd (index --stpd); alone, reads no BWT\n"
    "                      (--stpd --parse-oracle: the STPD sample is built reading the text through PREFIX.ptxt,\n"
    "                      not holding PREFIX.S: for collections whose text does not fit in memory)\n"
    "  --twolevel-locate   PREFIX.tloc: samples at the run starts of the BWT, the parse BWT and the block BWT,\n"
    "                      block -> phrase and phrase -> text maps, and phi^-1, for find --twolevel --locate\n"
    "                      (needs PREFIX.lvl with blocks, the BWT; builds PREFIX.ptxt if missing)\n"
    "  --rlzap[=J,D,M]     PREFIX.rptxt: the parse oracle with the parse as an RLZAP-style factorization against\n"
    "                      the first sequence's (both strands) parse (default 4,4,2), for find --rlzap-oracle;\n"
    "                      from PREFIX.ptxt (built first if missing); alone, reads no BWT\n"
    "  --parse-oracle      the text read through the two-level parse (PREFIX.ptxt: phrase starts, parse,\n"
    "                      dictionary), for find --stpd --parse-oracle; uses --syncmer; alone, it reads only\n"
    "                      PREFIX.S (no BWT)\n";

static void build_kebab(const string &Sfile, uint32_t k, double eps, const string &out) {
    kmer_filter F; F.k = k;
    double t0 = now();
    auto stream = [&](auto f) {
        FILE *in = open_or_die(Sfile, "rb");
        const size_t B = 1 << 24;
        std::vector<unsigned char> buf(B + 32);
        size_t keep = 0, got;
        while ((got = fread(buf.data() + keep, 1, B, in)) > 0) {
            size_t len = keep + got;
            F.for_each_kmer(buf.data(), len, f);
            keep = std::min<size_t>(k - 1, len);
            memmove(buf.data(), buf.data() + len - keep, keep);
        }
        fclose(in);
    };
    const int P = 16;                                   // HyperLogLog estimate of the distinct k-mers
    std::vector<uint8_t> reg(1u << P, 0);
    uint64_t total = 0;
    stream([&](uint64_t, uint64_t c) {
        uint64_t x = kmer_filter::mix(c ^ 0x5bd1e995ULL);
        uint32_t j = x >> (64 - P);
        uint8_t r = (uint8_t)(__builtin_clzll((x << P) | (1ULL << (P - 1))) + 1);
        if (r > reg[j]) reg[j] = r;
        ++total;
    });
    double m = 1u << P, z = 0; uint64_t zeros = 0;
    for (auto r : reg) { z += std::ldexp(1.0, -r); if (!r) ++zeros; }
    double est = 0.7213 / (1 + 1.079 / m) * m * m / z;
    if (est < 2.5 * m && zeros) est = m * std::log(m / zeros);
    double bpk = -std::log(eps) / (std::log(2.0) * std::log(2.0)) * 1.07;
    uint32_t h = std::max(1u, (uint32_t)std::lround(std::log(2.0) * bpk));
    F.init((uint64_t)est, bpk, h);
    stream([&](uint64_t, uint64_t c) { F.add(c); });
    if (!F.save(out)) die("cannot write " + out);
    fprintf(stderr, "KeBaB k=%u: %lu k-mers, ~%.0f distinct; %.2f bits/k-mer, h=%u, %.3f GB, %.0f s\n", k, total, est, bpk, h, F.bytes() / 1e9, now() - t0);
}

static void save_rlbwt(const rlbwt &B, const string &file) {
    FILE *f = open_or_die(file, "wb"); wmagic(f, "KTRLBW01"); B.save(f); fclose(f);
}
static void load_rlbwt(rlbwt &B, const string &file) {
    FILE *f = open_or_die(file, "rb"); rmagic(f, "KTRLBW01", file); B.load(f); fclose(f);
}

static int cmd_index(int argc, char **argv) {
    argv_holder A(expand_preset(argc, argv, &preset::index));
    string bwtfile, tagout, bigrep;
    bool csa = false, efcsa = false, tags = false, by_tag = false, no_rmq = false, grammar = false, tl = false;
    u64 sample = 1; uint32_t kebab = 0; double eps = 0.1;
    twolevel_options to;
    static option longopts[] = {
        {"bwt", 1, 0, 'b'}, {"csa", 0, 0, 'c'}, {"efcsa", 0, 0, 'e'}, {"tags", 0, 0, 't'}, {"by-tag", 0, 0, 'S'},
        {"sample", 1, 0, 's'}, {"no-rmq", 0, 0, 'n'}, {"grammar", 0, 0, 'g'}, {"bigrepair", 1, 0, 'G'}, {"tag-out", 1, 0, 'o'},
        {"kebab", 1, 0, 'k'}, {"kebab-eps", 1, 0, 'E'}, {"twolevel", 0, 0, '2'}, {"syncmer", 1, 0, 'y'}, {"no-blocks", 0, 0, 'B'},
        {"continuation", 1, 0, 'C'}, {"max-block", 1, 0, 1004}, {"tag-lcp", 1, 0, 'l'}, {"tmp", 1, 0, 1001}, {"blam-mem", 1, 0, 1002}, {"blam-external", 0, 0, 1003}, {"parse-oracle", 0, 0, 'O'}, {"stpd-parse", 0, 0, 'Y'}, {"rlzap", 2, 0, 'Z'}, {"twolevel-locate", 0, 0, 'W'}, {"stpd", 0, 0, 'P'}, {"stpd-all", 0, 0, 'Q'}, {"help", 0, 0, 'h'}, {0, 0, 0, 0}};
    int c;
    int stpd = 0;
    bool poracle = false, pstpd = false, rlzap = false, tlocate = false;
    u64 rz_J = 4, rz_M = 2; int64_t rz_D = 4;
    std::vector<u64> lambdas;
    u64 blam_mem = 6ull << 30;
    bool blam_ext = false;
    optind = 1;
    while ((c = getopt_long(A.argc(), A.argv(), "h", longopts, nullptr)) != -1) {
        switch (c) {
            case 'b': bwtfile = optarg; break;
            case 'c': csa = true; break;
            case 'e': efcsa = true; break;
            case 't': tags = true; break;
            case 'S': by_tag = true; break;
            case 's': sample = std::stoull(optarg); break;
            case 'n': no_rmq = true; break;
            case 'g': grammar = true; break;
            case 'G': bigrep = optarg; break;
            case 'o': tagout = optarg; break;
            case 'k': kebab = std::stoul(optarg); break;
            case 'E': eps = std::stod(optarg); break;
            case '2': tl = true; break;
            case 'y': if (sscanf(optarg, "%d,%d", &to.k, &to.s) != 2) die("--syncmer k,s"); break;
            case 'B': to.blocks = false; break;
            case 'P': stpd = std::max(stpd, 1); break;
            case 1001: tmp_dir() = optarg; break;
            case 1002: blam_mem = (u64)(std::stod(optarg) * (1ull << 30)); break;
            case 1003: blam_ext = true; break;
            case 'l': { std::string v = optarg, w; for (char ch : v + ",") { if (ch == ',') { if (!w.empty()) lambdas.push_back(std::stoull(w)); w.clear(); } else w += ch; } break; }
            case 'Q': stpd = 2; break;
            case 'O': poracle = true; break;
            case 'Y': pstpd = true; break;
            case 'W': tlocate = true; break;
            case 'Z':
                rlzap = true;
                if (optarg && sscanf(optarg, "%lu,%ld,%lu", &rz_J, &rz_D, &rz_M) != 3) die("--rlzap[=J,D,M]");
                break;
            case 'C': { string v = optarg; if (v == "frequent") to.frequent = true; else if (v != "smallest") die("--continuation smallest|frequent"); break; }
            case 1004: to.max_block = (uint32_t)std::stoul(optarg); break;
            default: fputs(INDEX_USAGE, stderr); return 1;
        }
    }
    if (A.argc() - optind != 1) { fputs(INDEX_USAGE, stderr); return 1; }
    const string prefix = A.argv()[optind];
    if (bwtfile.empty()) bwtfile = prefix + ".bwt";
    const string Sfile = prefix + ".S";
    double t0 = now();
    auto build_parse_oracle = [&]() {
        if (read_info_digest(prefix)) die("the parse oracle is for undigested texts");
        reset_step_peak();
        parse_text PT;
        PT.build(Sfile, to.k, to.s, prefix + ".ptxt.tmp");
        PT.save(prefix + ".ptxt");
        fprintf(stderr, "parse oracle: %.3f MB (peak memory of this step %.2f GB)\n", PT.bytes() / 1e6, step_peak_gb());
    };
    auto build_parse_stpd = [&]() {
        int32_t k = 0, sm = 0;
        {   // k and s of PREFIX.lvl
            FILE *f = open_or_die(prefix + ".lvl", "rb"); char mg[8]; if (fread(mg, 1, 8, f) != 8 || memcmp(mg, "KTLVL00", 7)) die(prefix + ".lvl: not a two-level index"); rd(f, k); rd(f, sm); fclose(f);
        }
        reset_step_peak();
        parse_text PT;
        bool ok = false;
        if (parse_text::current_format(prefix + ".ptxt")) { PT.load(prefix + ".ptxt"); ok = (int32_t)PT.K == k && (int32_t)PT.SM == sm && PT.n == file_size(Sfile) + 1; }
        if (!ok) {
            PT = parse_text();
            PT.build(Sfile, k, sm, prefix + ".ptxt.tmp");
            PT.save(prefix + ".ptxt");
        }
        parse_stpd PS;
        PS.build(PT);
        PS.save(prefix + ".pstpd");
        fprintf(stderr, "STPD over the parse: %.3f MB (+%.3f MB for locating); parse oracle %.3f MB (peak memory of this step %.2f GB)\n",
                (PS.pos.bytes() + PS.entry.bytes() + PS.C.bytes()) / 1e6, PS.bytes_locate() / 1e6, PT.bytes() / 1e6, step_peak_gb());
    };
    auto build_rlzap = [&]() {
        reset_step_peak();
        if (!parse_text::current_format(prefix + ".ptxt")) build_parse_oracle();
        parse_text PT; PT.load(prefix + ".ptxt");
        genome_table G; G.load(prefix + ".tbl");
        PT.make_rlz(2 * (G.len[0] + 1), rz_J, rz_D, rz_M);     // the reference: the first sequence, both strands
        PT.save(prefix + ".rptxt");
        fprintf(stderr, "RLZAP parse oracle: %.3f MB (peak memory of this step %.2f GB)\n", PT.bytes() / 1e6, step_peak_gb());
    };
    if ((poracle || pstpd || rlzap) && !csa && !efcsa && !tags && !stpd && !tl && !kebab) {   // no BWT needed
        if (poracle) build_parse_oracle();
        if (pstpd) build_parse_stpd();
        if (rlzap) build_rlzap();
        fprintf(stderr, "index built in %.1f s (peak memory %.2f GB)\n", now() - t0, peak_gb());
        return 0;
    }
    rlbwt B;
    B.build(bwtfile);
    save_rlbwt(B, prefix + ".rlbwt");
    fprintf(stderr, "RLBWT: n=%lu r=%lu, %.3f MB (%.1f s)\n", B.n, B.R, B.bytes() / 1e6, now() - t0);
    // Each step frees what it built before the next one starts: only the RLBWT is kept throughout (the tags,
    // STPD and two-level steps need it).  The two-level step reloads the RLCSA it uses for character steps.
    if (csa) {
        rlcsa X; X.build(bwtfile, false); X.save(prefix + ".csa");
        fprintf(stderr, "explicit RLCSA: %.3f MB (peak memory so far %.2f GB)\n", X.bytes() / 1e6, peak_gb());
    }
    if (efcsa) {
        rlcsa Y; Y.build(bwtfile, true); Y.save(prefix + ".efcsa");
        fprintf(stderr, "Elias-Fano RLCSA: %.3f MB (peak memory so far %.2f GB)\n", Y.bytes() / 1e6, peak_gb());
    }
    if (tags) {
        double t1 = now();
        genome_table G; G.load(prefix + ".tbl");
        tag_runs TR;
        const string cache = prefix + (by_tag ? ".sp.tagruns" : ".tagruns");
        if (!lambdas.empty() || !TR.load(cache) || TR.n != B.n) {
            std::vector<u64> blam;
            TR.compute(Sfile, bwtfile, B, G, by_tag, lambdas, &blam, blam_mem, blam_ext);
            TR.save(cache);
            if (!lambdas.empty()) {                    // Lemma "runs": rho <= (r + sigma + 2d)(lambda + 1) + B_lambda + 1
                u64 sigma = 0, dseg = 0;
                for (int ch = 0; ch < 256; ++ch) sigma += B.C[ch + 1] > B.C[ch];
                std::vector<uint32_t> gt = G.tags(by_tag);
                for (size_t i = 0; i < gt.size(); ++i) dseg += i == 0 || gt[i] != gt[i - 1];
                fprintf(stderr, "tag runs rho = %zu; r = %lu, sigma = %lu, d = %lu segments\n", TR.start.size(), B.R, sigma, dseg);
                for (size_t i = 0; i < lambdas.size(); ++i)
                    fprintf(stderr, "  lambda = %lu: B_lambda = %lu, bound %lu\n", lambdas[i], blam[i],
                            (B.R + sigma + 2 * dseg) * (lambdas[i] + 1) + blam[i] + 1);
            }
        }
        tag_build_options o; o.s = sample; o.rmq = !no_rmq; o.grammar = grammar; o.bigrepair = bigrep;
        tag_array T;
        build_tag_array(TR, &B, o, T);
        T.save(tagout.empty() ? prefix + ".tag" : tagout);
        fprintf(stderr, "tag array: n=%lu, %lu tag runs, %.3f MB (%.1f s; peak memory so far %.2f GB)\n", T.n, T.rho, T.bytes() / 1e6, now() - t1, peak_gb());
    }
    if (stpd && poracle) {                          // the STPD sample built with the parse oracle instead of the text
        build_parse_oracle();
        reset_step_peak();
        parse_text PT; PT.load(prefix + ".ptxt");
        stpd_index D;
        D.build(Sfile, bwtfile, B, stpd == 1, &PT);
        D.save(prefix + ".stpd");
        fprintf(stderr, "STPD sample built through the parse oracle (peak memory of this step %.2f GB)\n", step_peak_gb());
    } else {
        if (stpd) {
            stpd_index D;
            D.build(Sfile, bwtfile, B, stpd == 1);
            D.save(prefix + ".stpd");
        }
        if (poracle) build_parse_oracle();
    }
    if (read_info_digest(prefix) && (kebab || tl)) die("KeBaB and two-level indexing are for undigested texts");
    if (kebab) build_kebab(Sfile, kebab, eps, prefix + ".k" + std::to_string(kebab) + ".kbf");
    if (tl) {
        reset_step_peak();
        rlcsa X;                                    // for the character steps: the one just built, or build it
        if (file_exists(prefix + ".csa")) X.load(prefix + ".csa");
        if (X.n != B.n) X.build(bwtfile, false);   // (none, or one left from another BWT)
        stepper st; st.bwt = &B; st.csa = &X;
        two_level L;
        to.tmp = prefix + ".lvl.tmp";
        build_twolevel(Sfile, B, st, to, L, [&]() { X = rlcsa(); B = rlbwt(); });   // (nothing after this step needs them)
        L.save(prefix + ".lvl");
        fprintf(stderr, "two-level: parse %.3f MB, blocks %.3f MB (peak memory of this step %.2f GB)\n", L.bytes_parse() / 1e6,
                L.has_blocks ? L.bytes_blocks() / 1e6 : 0.0, step_peak_gb());
    }
    if (pstpd) build_parse_stpd();
    if (rlzap) build_rlzap();
    if (tlocate) {                                  // after the two-level step, which may have freed the RLBWT
        reset_step_peak();
        if (B.n == 0) load_rlbwt(B, prefix + ".rlbwt");
        two_level L; L.load(prefix + ".lvl");
        parse_text PT;
        if (parse_text::current_format(prefix + ".ptxt")) PT.load(prefix + ".ptxt");
        if (PT.rlz || (int32_t)PT.K != L.K || (int32_t)PT.SM != L.SM || PT.n != L.n) {
            PT = parse_text();
            PT.build(Sfile, L.K, L.SM, prefix + ".ptxt.tmp");
            PT.save(prefix + ".ptxt");
        }
        twolevel_locate T;
        build_twolevel_locate(bwtfile, B, L, PT, T);
        T.save(prefix + ".tloc");
        fprintf(stderr, "two-level locate: %.3f MB (peak memory of this step %.2f GB)\n", (T.bytes_char() + T.bytes_parse() + T.bytes_blocks()) / 1e6, step_peak_gb());
    }
    fprintf(stderr, "index built in %.1f s (peak memory %.2f GB)\n", now() - t0, peak_gb());
    return 0;
}

// ---------------------------------------------------------------------------------------------------
// katka find

static const char *FIND_USAGE =
    "usage: katka find [--preset NAME] [options] PREFIX reads out\n"
    "  reads: FASTQ or FASTA (several lines per record); out: '-' for standard output\n"
    "  matches\n"
    "    -L n              MEMs of at least n bases (with -X or --top: the starting value; default 1)\n"
    "    -X, --adaptive    adaptive L: keep L equal to the longest MEM found so far\n"
    "    --top t           report only the t longest MEMs (L = 1 + the shortest of the t kept)\n"
    "    --longest t       report the t longest MEMs with ties: every MEM at least as long as the t-th longest\n"
    "                      (with -K: pseudo-MEMs longest first, until the next is shorter than L; exact if\n"
    "                      L >= k; overrides -X)\n"
    "    --kebab-top t     with -K: --longest t over the t longest pseudo-MEMs only (KeBaB's heuristic)\n"
    "    --kebab-batch t   with -K: --longest t over batches of pseudo-MEMs (2t, then 4t, 8t, ...) until t\n"
    "                      MEMs are longer than the next pseudo-MEM; BML on the whole read if they run out\n"
    "                      with fewer than t MEMs of length >= k and L < k (exact)\n"
    "    --kebab-check T t with -K: --longest t over the T longest pseudo-MEMs (T >= t), accepted if there are at\n"
    "                      least T + 1 pseudo-MEMs and at least t MEMs are at least as long as the (T+1)st\n"
    "                      (--strict: longer); otherwise BML on the whole read (exact)\n"
    "    --strict          with --kebab-check: \"longer than\" (complete on ties)\n"
    "    --check-seed      with --kebab-check: the fallback starts with L = the t-th longest MEM found in the\n"
    "                      scanned pseudo-MEMs, if at least t were found (same output, fewer steps)\n"
    "    --fixed-L         with --kebab-top t: every MEM of length >= max(L, k) in the t longest pseudo-MEMs, as\n"
    "                      KeBaB's pipeline (no longest-t selection)\n"
    "    --read-stats FILE with --longest and -K: per read, name, pseudo-MEMs, searched, batches, fallback\n"
    "                      (0/1), backward steps, MEMs reported (TSV, in input order)\n"
    "    --fb              forward-backward search instead of BML\n"
    "    -D d              ignore MEMs with fewer than d distinct 3-mers\n"
    "    -U / --no-hpc     homopolymer-compress the reads (default: as the index was built)\n"
    "    -K k              KeBaB pseudo-MEMs, with PREFIX.k<k>.kbf (with -X or --top: longest first)\n"
    "    --inline-kebab k  KeBaB window check inside BML, with PREFIX.k<k>.kbf\n"
    "    --kbf FILE        the KeBaB filter file, instead of PREFIX.k<k>.kbf\n"
    "  search\n"
    "    --csa / --efcsa   backward steps on the explicit / Elias-Fano RLCSA (default: the RLBWT)\n"
    "    -F t              lookup table of the intervals of all t-mers (t <= 13)\n"
    "    --twolevel        two-level search (PREFIX.lvl)\n"
    "    --no-jumps        two-level search without block steps\n"
    "    --stpd            STPD search (PREFIX.stpd and the text PREFIX.S): binary searches over a sample of\n"
    "                      at most r positions; the intervals of reported MEMs by backward search\n"
    "    --no-intervals    with --stpd: do not compute intervals (mems format prints - -; no tags)\n"
    "    --stpd-parse      STPD over the parse (PREFIX.lvl's B and dictionary, PREFIX.pstpd, PREFIX.ptxt as the\n"
    "                      text, PREFIX.stpd at the left end): implies --stpd --parse-oracle; with --locate,\n"
    "                      phi^-1 over the parse, and the mems format prints occurrences and distinct sequences\n"
    "    --locate-out FILE with --locate: per MEM, name, start, end and its sorted positions (those of its reverse\n"
    "                      complement in S, for BML and FB), comma-separated, in output order\n"
    "    --rlzap-oracle    with --parse-oracle or --stpd-parse: the RLZAP parse oracle PREFIX.rptxt (index --rlzap)\n"
    "    --parse-oracle    with --stpd: read the text through the two-level parse (PREFIX.ptxt, from index\n"
    "                      --parse-oracle) instead of PREFIX.S\n"
    "    --locate          with --stpd: locate every occurrence of each MEM by phi^-1 steps checked against the\n"
    "                      text (doubling batches, then a binary search), and print their number (mems format);\n"
    "                      no intervals are computed unless tags are asked for (the mems format prints - -)\n"
    "                      (--stpd --twolevel: STPD finds the MEMs, two-level search their intervals)\n"
    "  answers (Cliffy's listing format, with counts)\n"
    "    -c, --count       each tag with its number of occurrences (proportional credit)\n"
    "    -l, --list        the distinct tags (even credit)\n"
    "    --range           the smallest and largest tags (range scheme, LCA) [default with tags]\n"
    "    -T FILE           the tag array (default PREFIX.tag)\n"
    "    --format F        listing (default), or mems: name, start, end, sp, ep per MEM, no tags\n"
    "    --classify bayes  per read, the genus by naive Bayes over its MEMs (flat prior): name, then the top documents\n"
    "                      (0-based tags, comma-separated if tied; - without MEMs); needs --nseq FILE (sequences per\n"
    "                      document, one per line) and counting tags; --epsilon E for a missing MEM (default 1e-6)\n"
    "    --mate2 FILE      --classify bayes: FILE holds each read's mate; a pair is classified over both mates' MEMs,\n"
    "                      under the first mate's name\n"
    "    --exclude-docs FILE  --classify bayes: 0-based documents that can't be called (score_silva.py --exclude-out)\n"
    "  -j p                p threads (the output is the same as with one)\n"
    "  --stats             search step counts\n";

enum answer_mode { RANGE, LIST, COUNT, NONE };

struct finder {
    stepper S;
    const lookup_table *table = nullptr;
    const two_level *TL = nullptr;
    bool jumps = true;
    const tag_array *tags = nullptr;
    mem_options opt;
    answer_mode mode = RANGE;
    bool hpc = false, mems_format = false;
    // --classify bayes: naive Bayes over the read's MEMs with a flat prior: genus (document) g's score is the sum
    // over the MEMs of log(c_g / (2 N_g)), or log epsilon where g lacks the MEM (as scripts/score_silva.py --mode
    // bayes --prior flat, in the same order of operations, so ties come out the same).  With --mate2, the MEMs of
    // both mates, mate 1's first; --exclude-docs: documents that can't be called (the scorer's unscored clades)
    bool classify = false;
    std::vector<u64> nseq;                    // N_g per document (--nseq)
    std::vector<bool> excluded;               // per document (--exclude-docs; empty: none)
    double logeps = 0;
    const unsigned char *comp = nullptr;      // digests: complement of each symbol
    const stpd_index *SP = nullptr;
    const parse_stpd *PS = nullptr;
    const twolevel_locate *TLC = nullptr;
    const std::vector<u64> *gstarts = nullptr;
    bool intervals = true, locate = false;

    void answer(const interval &I, string &out) const {
        char b[48];
        if (mode == NONE) return;
        if (mode == COUNT) {
            static thread_local std::vector<std::pair<u64, u64>> tc;
            tags->count(I.sp, I.ep, tc);
            for (u64 q = 0; q < tc.size(); ++q) { snprintf(b, sizeof b, q ? ",%lu:%lu" : "%lu:%lu", tc[q].first, tc[q].second); out += b; }
            return;
        }
        static thread_local std::vector<u64> tl;
        tl.clear(); tags->list(I.sp, I.ep, tl); std::sort(tl.begin(), tl.end());
        if (mode == LIST) for (u64 q = 0; q < tl.size(); ++q) { snprintf(b, sizeof b, q ? ",%lu" : "%lu", tl[q]); out += b; }
        else { snprintf(b, sizeof b, "%lu,%lu", tl.front(), tl.back()); out += b; }
    }

    bool read_stats = false;

    // the output of one read; returns its number of MEMs (with read_stats, its statistics line in *sl)
    u64 run(segment_search &Q, const string &name, const string &read, string &line, string *sl = nullptr, string *pl = nullptr, const string *mate2 = nullptr) const {
        if (pl) pl->clear();
        line.clear();
        const search_counters before = Q.cnt;
        static thread_local string h;
        static thread_local std::vector<u64> pos;
        if (hpc) {
            h.clear(); pos.clear();
            for (u64 i = 0; i < read.size(); ++i) if (i == 0 || read[i] != read[i - 1]) { h.push_back(read[i]); pos.push_back(i); }
            pos.push_back(read.size());
        }
        u64 nm = 0;
        char b[160];
        if (!mems_format && !classify) line = ">" + name + "\n";
        static thread_local std::vector<std::pair<u64, double>> ll;   // (document, sum of log p) in first-seen order
        static thread_local std::vector<u64> hits, slot;
        if (classify) { ll.clear(); hits.clear(); if (slot.size() < nseq.size()) slot.assign(nseq.size(), ~0ull); }
        auto emit = [&](const mem &M) {
            if (classify) {
                static thread_local std::vector<std::pair<u64, u64>> tc;
                tags->count(M.I.sp, M.I.ep, tc);
                for (auto &[doc, c] : tc) {
                    if (doc >= nseq.size()) die("--classify bayes: a tag beyond --nseq");
                    if (slot[doc] == ~0ull) { slot[doc] = ll.size(); ll.push_back({doc, 0.0}); hits.push_back(0); }
                    ll[slot[doc]].second += std::log((double)c / (double)(2 * nseq[doc]));
                    ++hits[slot[doc]];
                }
                ++nm;
                return;
            }
            u64 first = hpc ? pos[M.a] : M.a, last = hpc ? pos[M.e] - 1 : M.e - 1;
            if (pl) {                               // --locate-out: name, start, end, the sorted positions (of rc(MEM) in S)
                static thread_local std::vector<u64> ps;
                ps = Q.occbuf;
                std::sort(ps.begin(), ps.end());
                char hb[64];
                snprintf(hb, sizeof hb, "\t%lu\t%lu\t", first, last + 1);
                *pl += name; *pl += hb;
                for (size_t q = 0; q < ps.size(); ++q) { snprintf(hb, sizeof hb, q ? ",%lu" : "%lu", ps[q]); *pl += hb; }
                *pl += "\n";
            }
            if (mems_format) {
                if (intervals) snprintf(b, sizeof b, "\t%lu\t%lu\t%lu\t%lu", first, last + 1, M.I.sp, M.I.ep);
                else snprintf(b, sizeof b, "\t%lu\t%lu\t-\t-", first, last + 1);
                line += name; line += b;
                if (locate) { snprintf(b, sizeof b, gstarts ? "\t%lu\t%lu" : "\t%lu", M.occ, M.haps); line += b; }
                line += "\n";
            } else {
                snprintf(b, sizeof b, "[%lu,%lu] {", first, last);
                line += b; answer(M.I, line); line += "} ";
            }
            ++nm;
        };
        find_mems(Q, hpc ? h : read, opt, emit);
        if (mate2) find_mems(Q, *mate2, opt, emit);
        if (classify) {             // the top documents (0-based, sorted, comma-separated), or - for a read without MEMs
            double top = -INFINITY;
            static thread_local std::vector<double> sc;
            sc.resize(ll.size());
            auto scored = [&](size_t q) { return excluded.empty() || !excluded[ll[q].first]; };
            for (size_t q = 0; q < ll.size(); ++q) { sc[q] = 0.0 + ll[q].second + (double)(nm - hits[q]) * logeps; if (scored(q)) top = std::max(top, sc[q]); }
            static thread_local std::vector<u64> best;
            best.clear();
            for (size_t q = 0; q < ll.size(); ++q) if (scored(q) && sc[q] == top) best.push_back(ll[q].first);
            std::sort(best.begin(), best.end());
            line = name + "\t";
            if (best.empty()) line += "-";
            for (size_t q = 0; q < best.size(); ++q) { snprintf(b, sizeof b, q ? ",%lu" : "%lu", best[q]); line += b; }
            line += "\n";
            for (auto &x : ll) slot[x.first] = ~0ull;
        } else if (!mems_format) line += "\n";
        if (sl) {
            const search_counters &c = Q.cnt;
            const char *why = c.kb_fb_few > before.kb_fb_few ? "few_pseudo_mems" : c.kb_fb_short > before.kb_fb_short ? "few_long_mems"
                            : c.kb_fallbacks > before.kb_fallbacks ? "pseudo_mems_ran_out" : "-";
            snprintf(b, sizeof b, "\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%s\n", c.kb_pseudo - before.kb_pseudo, c.kb_searched - before.kb_searched,
                     c.kb_batches - before.kb_batches, c.kb_fallbacks - before.kb_fallbacks, c.char_steps - before.char_steps, nm, why);
            *sl = name + b;
        }
        return nm;
    }
};

struct read_parser {
    FILE *in = nullptr;
    gzFile gz = nullptr;
    bool fasta = false;
    string pending, buf;
    std::vector<char> raw;
    size_t rp = 0, rn = 0;
    explicit read_parser(const string &f) : raw(1 << 20) {
        gz = gzopen(f.c_str(), "rb");
        if (!gz) die("cannot open " + f);
        int c = gzgetc(gz);
        fasta = c == '>';
        if (c >= 0) gzungetc(c, gz);
    }
    ~read_parser() { if (gz) gzclose(gz); }
    bool getline(string &s) {
        s.clear();
        for (;;) {
            if (rp == rn) { int g = gzread(gz, raw.data(), raw.size()); if (g <= 0) return !s.empty(); rn = g; rp = 0; }
            char *p = raw.data() + rp, *e = raw.data() + rn, *q = (char *)memchr(p, '\n', e - p);
            if (q) { s.append(p, q); rp = q - raw.data() + 1; if (!s.empty() && s.back() == '\r') s.pop_back(); return true; }
            s.append(p, e); rp = rn;
        }
    }
    bool next(string &name, string &seq) {
        string h, line;
        if (!fasta) {
            string plus, qual;
            if (!(getline(h) && getline(seq) && getline(plus) && getline(qual))) return false;
        } else {
            if (pending.empty()) {
                while (getline(line)) if (!line.empty() && line[0] == '>') { pending = line; break; }
                if (pending.empty()) return false;
            }
            h.swap(pending); pending.clear(); seq.clear();
            while (getline(line)) {
                if (!line.empty() && line[0] == '>') { pending = line; break; }
                seq += line;
            }
        }
        for (auto &c : seq) c = (char)toupper((unsigned char)c);
        size_t blank = h.find_first_of(" \t");
        name = h.substr(1, blank == string::npos ? string::npos : blank - 1);
        return true;
    }
};

struct run_stats { u64 reads = 0, mems = 0, unclassified = 0; double read_time = 0; search_counters sc; };   // read_time: the sum of each read's own time

// threads take chunks of reads and write their outputs back in input order
static run_stats find_all(const finder &K, read_parser &P, FILE *out, int nthreads, FILE *sout = nullptr, FILE *lout = nullptr, read_parser *P2 = nullptr) {
    const size_t CHUNK = 256;
    std::mutex in_m, out_m, st_m;
    std::condition_variable out_cv;
    u64 next_in = 0, next_out = 0;
    bool eof = false;
    run_stats total;
    auto work = [&]() {
        std::vector<string> names, seqs, seqs2, lines, slines, plines;
        run_stats st;
        string name, seq, name2, seq2;
        segment_search Q;
        Q.S = &K.S; Q.T = K.table; Q.TL = K.TL; Q.jumps = K.jumps; Q.comp = K.comp; Q.SP = K.SP; Q.intervals = K.intervals; Q.locate = K.locate;
        Q.PS = K.PS; Q.gstarts = K.gstarts; Q.TLC = K.TLC;
        for (;;) {
            u64 my;
            {
                std::lock_guard<std::mutex> g(in_m);
                names.clear(); seqs.clear(); seqs2.clear();
                while (!eof && names.size() < CHUNK) {
                    if (!P.next(name, seq)) {
                        eof = true;
                        if (P2 && P2->next(name2, seq2)) die("--mate2: more reads than in the first file");
                        break;
                    }
                    names.push_back(name); seqs.push_back(seq);
                    if (P2) {
                        if (!P2->next(name2, seq2)) die("--mate2: fewer reads than in the first file");
                        seqs2.push_back(seq2);
                    }
                }
                my = next_in++;
            }
            lines.resize(names.size()); slines.resize(names.size()); plines.resize(names.size());
            for (size_t i = 0; i < names.size(); ++i) {
                const double t0r = now();
                u64 m = K.run(Q, names[i], seqs[i], lines[i], sout ? &slines[i] : nullptr, lout ? &plines[i] : nullptr, P2 ? &seqs2[i] : nullptr);
                st.read_time += now() - t0r;
                ++st.reads; st.mems += m; st.unclassified += m == 0;
            }
            {
                std::unique_lock<std::mutex> lk(out_m);
                out_cv.wait(lk, [&] { return next_out == my; });
                for (size_t i = 0; i < names.size(); ++i) fputs(lines[i].c_str(), out);
                if (sout) for (size_t i = 0; i < names.size(); ++i) fputs(slines[i].c_str(), sout);
                if (lout) for (size_t i = 0; i < names.size(); ++i) fputs(plines[i].c_str(), lout);
                ++next_out;
            }
            out_cv.notify_all();
            if (names.empty()) break;
        }
        std::lock_guard<std::mutex> g(st_m);
        total.reads += st.reads; total.mems += st.mems; total.unclassified += st.unclassified; total.read_time += st.read_time; total.sc.add(Q.cnt);
    };
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back(work);
    for (auto &x : th) x.join();
    return total;
}

static int cmd_find(int argc, char **argv) {
    argv_holder A(expand_preset(argc, argv, &preset::find));
    finder K;
    string tagfile, kbffile, format = "listing";
    int nthreads = 1, ft = 0, hpc = -1, csa = 0;
    uint32_t kebab_k = 0;
    bool twolevel = false, stats = false, use_stpd = false, poracle = false, kebab_rule_set = false, use_pstpd = false, rlzap_oracle = false;
    string readstats, locate_out;
    bool fixed_L = false;
    double epsilon = 1e-6;
    string mate2file, excludefile;
    int mode = -1;
    static option longopts[] = {
        {"adaptive", 0, 0, 'X'}, {"top", 1, 0, 1}, {"longest", 1, 0, 18}, {"kebab-top", 1, 0, 19}, {"kebab-batch", 1, 0, 20},
        {"read-stats", 1, 0, 21}, {"kebab-check", 1, 0, 23}, {"strict", 0, 0, 24}, {"check-seed", 0, 0, 26}, {"fixed-L", 0, 0, 25}, {"fb", 0, 0, 2}, {"no-hpc", 0, 0, 3}, {"inline-kebab", 1, 0, 4}, {"kbf", 1, 0, 5},
        {"csa", 0, 0, 6}, {"efcsa", 0, 0, 7}, {"twolevel", 0, 0, 8}, {"no-jumps", 0, 0, 9}, {"count", 0, 0, 'c'}, {"list", 0, 0, 'l'},
        {"range", 0, 0, 10}, {"format", 1, 0, 11}, {"stats", 0, 0, 12}, {"no-tags", 0, 0, 13}, {"stpd", 0, 0, 14},
        {"no-intervals", 0, 0, 15}, {"locate", 0, 0, 16}, {"parse-oracle", 0, 0, 17}, {"stpd-parse", 0, 0, 22}, {"rlzap-oracle", 0, 0, 27}, {"locate-out", 1, 0, 28}, {"classify", 1, 0, 29}, {"nseq", 1, 0, 30}, {"epsilon", 1, 0, 31}, {"mate2", 1, 0, 32}, {"exclude-docs", 1, 0, 33}, {"help", 0, 0, 'h'}, {0, 0, 0, 0}};
    int c;
    optind = 1;
    while ((c = getopt_long(A.argc(), A.argv(), "L:XD:UK:F:clT:j:h", longopts, nullptr)) != -1) {
        switch (c) {
            case 'L': K.opt.L = std::stoull(optarg); break;
            case 'X': K.opt.adaptive = true; break;
            case 1: K.opt.top = std::stoull(optarg); break;
            case 18: K.opt.longest = std::stoull(optarg); K.opt.kebab_rule = mem_options::STREAM; break;
            case 19: K.opt.longest = std::stoull(optarg); K.opt.kebab_rule = mem_options::TOP; kebab_rule_set = true; break;
            case 20: K.opt.longest = std::stoull(optarg); K.opt.kebab_rule = mem_options::BATCH; kebab_rule_set = true; break;
            case 21: readstats = optarg; break;
            case 23:                                         // --kebab-check T t
                if (optind >= A.argc()) die("--kebab-check T t");
                K.opt.check_T = std::stoull(optarg); K.opt.longest = std::stoull(A.argv()[optind++]);
                K.opt.kebab_rule = mem_options::CHECK; kebab_rule_set = true;
                if (K.opt.check_T < K.opt.longest || K.opt.longest == 0) die("--kebab-check T t needs T >= t >= 1");
                break;
            case 24: K.opt.check_strict = true; break;
            case 26: K.opt.check_seed = true; break;
            case 25: fixed_L = true; break;
            case 2: K.opt.fb = true; break;
            case 'D': K.opt.min3 = std::stoi(optarg); break;
            case 'U': hpc = 1; break;
            case 3: hpc = 0; break;
            case 'K': kebab_k = std::stoul(optarg); K.opt.pseudo = true; K.opt.inline_check = false; break;
            case 4: kebab_k = std::stoul(optarg); K.opt.inline_check = true; K.opt.pseudo = false; break;
            case 5: kbffile = optarg; break;
            case 6: csa = 1; break;
            case 7: csa = 2; break;
            case 8: twolevel = true; break;
            case 9: K.jumps = false; break;
            case 'c': mode = COUNT; break;
            case 'l': mode = LIST; break;
            case 10: mode = RANGE; break;
            case 'T': tagfile = optarg; break;
            case 11: format = optarg; break;
            case 'F': ft = std::stoi(optarg); break;
            case 'j': nthreads = std::max(1, std::stoi(optarg)); break;
            case 12: stats = true; break;
            case 29: if (string(optarg) != "bayes") die("--classify bayes"); K.classify = true; break;
            case 30: { std::ifstream in(optarg); if (!in) die(string("cannot open ") + optarg); u64 x; while (in >> x) K.nseq.push_back(x); break; }
            case 31: epsilon = std::stod(optarg); break;
            case 32: mate2file = optarg; break;
            case 33: excludefile = optarg; break;
            case 13: mode = NONE; break;
            case 14: use_stpd = true; break;
            case 17: poracle = true; break;
            case 22: use_stpd = poracle = twolevel = use_pstpd = true; break;
            case 27: rlzap_oracle = true; break;
            case 28: locate_out = optarg; break;
            case 15: K.intervals = false; break;
            case 16: K.locate = true; break;
            default: fputs(FIND_USAGE, stderr); return 1;
        }
    }
    if (A.argc() - optind != 3) { fputs(FIND_USAGE, stderr); return 1; }
    const string prefix = A.argv()[optind], readfile = A.argv()[optind + 1], outfile = A.argv()[optind + 2];
    if (K.opt.L == 0) K.opt.L = 1;
    if (ft < 0 || ft > 13) die("-F needs t <= 13");
    if (format != "listing" && format != "mems") die("--format listing|mems");
    K.mems_format = format == "mems";
    if (!(epsilon > 0)) die("--epsilon must be positive");
    K.logeps = std::log(epsilon);
    if ((!mate2file.empty() || !excludefile.empty()) && !K.classify) die("--mate2 and --exclude-docs need --classify bayes");
    if (K.classify) {
        if (K.nseq.empty()) die("--classify bayes needs --nseq FILE (sequences per document, one per line)");
        if (!excludefile.empty()) {
            std::ifstream in(excludefile);
            if (!in) die("cannot open " + excludefile);
            K.excluded.assign(K.nseq.size(), false);
            u64 x;
            while (in >> x) { if (x >= K.nseq.size()) die("--exclude-docs: a document beyond --nseq"); K.excluded[x] = true; }
        }
        if (K.mems_format) die("--classify bayes writes its own format (name, top documents)");
        mode = COUNT;
    }
    K.hpc = hpc < 0 ? read_info_hpc(prefix) : hpc == 1;
    if (hpc == 1 && !read_info_hpc(prefix) && file_exists(prefix + ".info")) fprintf(stderr, "warning: -U, but the index was built without homopolymer compression\n");
    if (K.opt.adaptive && K.opt.top) die("-X and --top are exclusive");
    if (fixed_L) {                                           // --kebab-top t --fixed-L: KeBaB's pipeline
        if (K.opt.kebab_rule != mem_options::TOP || !K.opt.longest) die("--fixed-L goes with --kebab-top t");
        K.opt.fixed_top = K.opt.longest; K.opt.longest = 0; K.opt.adaptive = false;
    }
    if ((K.opt.check_strict || K.opt.check_seed) && K.opt.kebab_rule != mem_options::CHECK) die("--strict and --check-seed go with --kebab-check");
    if (K.opt.longest && K.opt.top) die("--top and --longest/--kebab-top/--kebab-batch/--kebab-check are exclusive");
    if (K.opt.longest) K.opt.adaptive = false;
    if (kebab_rule_set && !K.opt.pseudo) die("--kebab-top and --kebab-batch need -K");
    if (!readstats.empty() && !((K.opt.longest || K.opt.fixed_top) && K.opt.pseudo)) die("--read-stats needs --longest, --kebab-top or --kebab-batch, with -K");

    double t0 = now();
    rlbwt B; load_rlbwt(B, prefix + ".rlbwt");
    rlcsa X;
    if (csa) X.load(prefix + (csa == 1 ? ".csa" : ".efcsa"));
    K.S.bwt = &B; K.S.csa = csa ? &X : nullptr;
    tag_array T;
    if (K.mems_format) mode = NONE;
    else {
        if (tagfile.empty()) tagfile = prefix + ".tag";
        if (mode == -1) mode = file_exists(tagfile) ? RANGE : NONE;
        if (mode != NONE) { T.load(tagfile, mode != COUNT); T.bwt = &B; K.tags = &T; }
    }
    K.mode = (answer_mode)mode;
    kmer_filter KF;
    if (kebab_k) {
        if (kbffile.empty()) kbffile = prefix + ".k" + std::to_string(kebab_k) + ".kbf";
        if (!KF.load(kbffile)) die("cannot load " + kbffile + " (katka index --kebab " + std::to_string(kebab_k) + ")");
        if (KF.k > K.opt.L && !K.opt.adaptive && !K.opt.top && !K.opt.longest && K.opt.pseudo) die("-K: the filter's k exceeds L");
        K.opt.filter = &KF;
    }
    digest_map DM;
    const bool digested = read_info_digest(prefix);
    if (digested) {
        DM.load(prefix + ".map");
        if (K.opt.adaptive || K.opt.top || K.opt.fb || kebab_k || twolevel || K.hpc)
            die("over digests, only BML with a fixed L (and -D) is supported: no -X, --top, --fb, KeBaB, two-level or -U");
        if (K.opt.L < (u64)DM.w) die("over digests, L must be at least w = " + std::to_string(DM.w));
        K.opt.digest = &DM;
        K.comp = DM.comp;
    }
    stpd_index SD;
    parse_text PT;
    if (poracle && !use_stpd) die("--parse-oracle needs --stpd");
    if (rlzap_oracle && !poracle) die("--rlzap-oracle goes with --stpd --parse-oracle or --stpd-parse");
    if (use_stpd) {
        SD.load(prefix + ".stpd", poracle ? string() : prefix + ".S");
        if (poracle) {
            PT.load(prefix + (rlzap_oracle ? ".rptxt" : ".ptxt"));
            if (rlzap_oracle && !PT.rlz) die(prefix + ".rptxt is not an RLZAP parse oracle (katka index --rlzap)");
            if (PT.n != SD.n) die(prefix + ".ptxt does not match " + prefix + ".stpd");
            SD.PT = &PT;
            fprintf(stderr, "parse oracle: %lu phrases (k = %u, s = %u), %u distinct: %.3f MB (%.2f bits per character) instead of the text (%.3f MB)\n",
                    PT.m, PT.K, PT.SM, PT.sigma, PT.bytes() / 1e6, 8.0 * PT.bytes() / PT.n, PT.n / 1e6);
        }
        K.SP = &SD;
        fprintf(stderr, "STPD: %lu sampled positions (r = %lu), %.3f MB on top of the text\n", SD.pos.size(), SD.runs, SD.bytes() / 1e6);
    }
    if (!K.intervals && (!use_stpd || K.mode != NONE)) die("--no-intervals needs --stpd and no tags (--format mems or --no-tags)");
    if (K.locate && !use_stpd && !twolevel) die("--locate needs --stpd, --stpd-parse or --twolevel");
    // locating lists the occurrences by phi^-1 from the first one and counts them, so the intervals (a
    // backward search over every MEM) are needed only for tags
    if (K.locate && K.mode == NONE && use_stpd) K.intervals = false;
    two_level L2;
    parse_stpd PS;
    genome_table GT;
    if (K.locate) {                     // the sequence starts, to count the distinct sequences of each MEM's occurrences
        GT.load(prefix + ".tbl");
        K.gstarts = &GT.start;
    }
    twolevel_locate TLc;
    if (twolevel) {
        L2.load(prefix + ".lvl", use_pstpd);
        K.TL = &L2;
        if (K.locate && !use_stpd) {
            TLc.load(prefix + ".tloc");
            if (TLc.n != L2.n || TLc.m != L2.m || TLc.nblocks != L2.nblocks || TLc.n != B.n) die(prefix + ".tloc does not match the index (katka index --twolevel-locate)");
            K.TLC = &TLc;
            fprintf(stderr, "two-level locate (MB): characters %.3f (run samples %.3f, phi^-1 %.3f), parse %.3f (run samples %.3f, phrase -> text %.3f), blocks %.3f (run samples %.3f, block -> phrase %.3f)\n",
                    TLc.bytes_char() / 1e6, TLc.csamp.bytes() / 1e6, (TLc.ikeys.bytes() + TLc.ivals.bytes()) / 1e6, TLc.bytes_parse() / 1e6,
                    (TLc.prow.bytes() + TLc.psamp.bytes()) / 1e6, TLc.pstarts.bytes() / 1e6, TLc.bytes_blocks() / 1e6,
                    (TLc.brow.bytes() + TLc.bsamp.bytes()) / 1e6, TLc.bstarts.bytes() / 1e6);
        }
        if (use_pstpd) {
            PS.load(prefix + ".pstpd");
            if (PS.m != PT.m || PS.K != PT.K || PS.SM != PT.SM || (int32_t)PT.K != L2.K || (int32_t)PT.SM != L2.SM || L2.n != SD.n)
                die(prefix + ".pstpd, .ptxt, .lvl and .stpd do not match (katka index --stpd-parse)");
            K.PS = &PS;
            const u64 rb = csa ? X.bytes() : B.bytes();
            if (PT.rlz)
                fprintf(stderr, "RLZAP parse (MB): reference ids %.3f + starts %.3f, factor starts %.3f (parse) + %.3f (text), literal counts %.3f, literals %.3f, pointers %.3f: %.3f in all (instead of the plain parse and boundaries)\n",
                        PT.Zr.ref.bytes() / 1e6, PT.Zr.refstart.bytes() / 1e6, PT.Zr.fpos.bytes() / 1e6, PT.Zr.ftext.bytes() / 1e6, PT.Zr.litcum.bytes() / 1e6,
                        PT.Zr.lits.bytes() / 1e6, (PT.Zr.chg.bytes() + PT.Zr.relv.bytes()) / 1e6, PT.Zr.bytes() / 1e6);
            fprintf(stderr, "STPD over the parse, components (MB): plain parse %.3f, boundaries %.3f, dictionary strings %.3f (2 bits per base, with exceptions; + offsets %.3f), "
                            "parse sample %.3f (%lu of r_Z = %lu), entry table %.3f (+ C %.3f), parse phi^-1 %.3f (locating), character-level sample %.3f "
                            "(+ phi^-1 %.3f, locating), right end: %s %.3f, B %.3f, dictionary hash %.3f\n",
                    PT.ids.bytes() / 1e6, PT.starts.bytes() / 1e6, PT.dict.bytes() / 1e6, PT.doff.bytes() / 1e6, PS.pos.bytes() / 1e6, PS.pos.size(), PS.runs,
                    PS.entry.bytes() / 1e6, PS.C.bytes() / 1e6, PS.bytes_locate() / 1e6, SD.bytes() / 1e6, SD.locate_bytes() / 1e6,
                    csa ? "RLCSA" : "RLBWT", rb / 1e6, L2.B.bytes() / 1e6, L2.dict.bytes() / 1e6);
        }
        if (ft) { fprintf(stderr, "note: -F is not used by two-level search\n"); ft = 0; }
    }
    lookup_table LT;
    if (ft) {
        if (digested) {          // the table over the digest symbols
            std::vector<unsigned char> a;
            for (int c = 3; c < 256; ++c) if (c != 'X' && B.C[c + 1] > B.C[c]) a.push_back((unsigned char)c);
            LT.set_alphabet(a);
        }
        LT.build(K.S, ft); K.table = &LT;
    }
    fprintf(stderr, "index loaded in %.1f s\n", now() - t0);

    if (!mate2file.empty() && K.hpc) die("--mate2 is for texts without homopolymer compression");
    read_parser P(readfile);
    std::unique_ptr<read_parser> P2;
    if (!mate2file.empty()) P2.reset(new read_parser(mate2file));
    FILE *out = outfile == "-" ? stdout : open_or_die(outfile, "w");
    double t1 = now();
    FILE *sout = nullptr;
    if (!readstats.empty()) { sout = open_or_die(readstats, "w"); fputs("read\tpseudo_mems\tsearched\tbatches\tfallback\tbackward_steps\tmems\tfallback_reason\n", sout); }
    if (!locate_out.empty() && !K.locate) die("--locate-out needs --locate");
    FILE *lout = locate_out.empty() ? nullptr : open_or_die(locate_out, "w");
    run_stats st = find_all(K, P, out, nthreads, sout, lout, P2.get());
    if (sout) fclose(sout);
    if (lout) fclose(lout);
    double wall = now() - t1;
    if (out != stdout) fclose(out);
    fprintf(stderr, "reads=%lu MEMs=%lu without MEMs=%lu | threads=%d wall %.2f s, %.1f us/read\n", st.reads, st.mems, st.unclassified, nthreads, wall,
            1e6 * wall * nthreads / std::max<u64>(1, st.reads));
    // the timing protocol's two numbers: each read timed alone (its search and output formatting, in its
    // thread), averaged; and the wall time of the whole search (reading, searching, writing; loading excluded)
    // divided by the number of reads
    fprintf(stderr, "timing: %lu reads, %d threads; mean per-read time %.2f us; total time / reads %.2f us (total %.3f s)\n",
            st.reads, nthreads, 1e6 * st.read_time / std::max<u64>(1, st.reads), 1e6 * wall / std::max<u64>(1, st.reads), wall);
    if (stats)
        fprintf(stderr, "searches=%lu char_steps=%lu parse_steps=%lu block_steps=%lu kmer_tests=%lu kmer_restarts=%lu stpd_jumps=%lu stpd_compares=%lu\n",
                st.sc.searches, st.sc.char_steps, st.sc.parse_steps, st.sc.block_steps, st.sc.kmer_tests, st.sc.kmer_restarts, st.sc.stpd_jumps, st.sc.stpd_compares);
    if (stats && st.sc.parse_steps)
            fprintf(stderr, "parse steps: entry=%lu inside=%lu end=%lu; block attempts=%lu, symbol absent=%lu, failed block steps=%lu, back-ups=%lu, exact exits=%lu\n",
                    st.sc.pe_entry, st.sc.pe_inside, st.sc.pe_end, st.sc.b_attempts, st.sc.b_nosym, st.sc.b_fail, st.sc.b_backups, st.sc.b_exact);
    if (stats && K.opt.fixed_top && K.opt.pseudo)
        fprintf(stderr, "kebab (top, fixed L = %lu, t = %lu): pseudo-MEMs %lu (%.2f per read), searched %lu (%.2f per read), backward steps %.1f per read\n",
                K.opt.L, K.opt.fixed_top, st.sc.kb_pseudo, (double)st.sc.kb_pseudo / std::max<u64>(1, st.reads), st.sc.kb_searched,
                (double)st.sc.kb_searched / std::max<u64>(1, st.reads), (double)st.sc.char_steps / std::max<u64>(1, st.reads));
    if (stats && K.opt.kebab_rule == mem_options::CHECK && K.opt.longest && K.opt.pseudo)
        fprintf(stderr, "kebab check (T = %lu, t = %lu%s%s): fallbacks %lu of %lu reads (%lu with at most T pseudo-MEMs, %lu with fewer than t long MEMs; %lu seeded)\n",
                K.opt.check_T, K.opt.longest, K.opt.check_strict ? ", strict" : "", K.opt.check_seed ? ", seeded" : "", st.sc.kb_fallbacks, st.reads,
                st.sc.kb_fb_few, st.sc.kb_fb_short, st.sc.kb_seeded);
    if (stats && K.opt.longest && K.opt.pseudo)
        fprintf(stderr, "kebab (%s, t = %lu): pseudo-MEMs %lu (%.2f per read), searched %lu (%.2f per read), batches %lu, reads with a second batch %lu, whole-read fallbacks %lu, backward steps %.1f per read\n",
                K.opt.kebab_rule == mem_options::TOP ? "top" : K.opt.kebab_rule == mem_options::BATCH ? "batch" : K.opt.kebab_rule == mem_options::CHECK ? "check" : "streaming", K.opt.longest,
                st.sc.kb_pseudo, (double)st.sc.kb_pseudo / std::max<u64>(1, st.reads), st.sc.kb_searched, (double)st.sc.kb_searched / std::max<u64>(1, st.reads),
                st.sc.kb_batches, st.sc.kb_multibatch, st.sc.kb_fallbacks, (double)st.sc.char_steps / std::max<u64>(1, st.reads));
    if (stats && use_pstpd)
        fprintf(stderr, "STPD over the parse: right-end character steps %lu, middle phrase comparisons %lu, parse-sample binary searches %lu (%lu comparisons), left-end character steps %lu (character-level STPD jumps %lu)\n",
                st.sc.ps_right, st.sc.ps_mid, st.sc.ps_jumps, st.sc.ps_compares, st.sc.ps_left, st.sc.stpd_jumps);
    if (stats && poracle) {
        const std::atomic<u64> *t = parse_text::totals();
        fprintf(stderr, "parse oracle seeks: same phrase %lu, next/previous phrase %lu, predecessor queries %lu\n", t[0].load(), t[1].load(), t[2].load());
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------------

static int cmd_info(int argc, char **argv) {
    if (argc != 2) { fputs("usage: katka info PREFIX\n", stderr); return 1; }
    const string p = argv[1];
    for (const char *s : {".S", ".bwt", ".rlbwt", ".csa", ".efcsa", ".tag", ".tag.rmq", ".lvl"})
        if (file_exists(p + s)) printf("%-12s %12.3f MB\n", s, file_size(p + s) / 1e6);
    printf("homopolymer-compressed: %s\n", read_info_hpc(p) ? "yes" : "no");
    if (read_info_digest(p)) { digest_map M; M.load(p + ".map"); printf("minimizer digest: k = %d, w = %d\n", M.k, M.w); }
    return 0;
}

// katka stpd-prep / stpd-export: input files for the published STPD-index (stpd_export.hpp)
static int cmd_stpd_prep(int argc, char **argv) {
    bool split = false; int a = 1;
    if (a < argc && std::string(argv[a]) == "-s") { split = true; ++a; }
    if (argc - a != 2) { fputs("usage: katka stpd-prep [-s] list.tsv OUT   (writes OUT.T, OUT.S = reverse(OUT.T), OUT.pieces, OUT.reflen)\n", stderr); return 1; }
    stpd_prep(argv[a], argv[a + 1], split);
    return 0;
}
static int cmd_stpd_export(int argc, char **argv) {
    if (argc != 2 && !(argc == 4 && std::string(argv[2]) == "--bwt")) { fputs("usage: katka stpd-export OUT [--bwt FILE]   (after katka bwt / Big-BWT on OUT.S and katka index --stpd OUT)\n", stderr); return 1; }
    stpd_export(argv[1], argc == 4 ? argv[3] : std::string(argv[1]) + ".bwt");
    return 0;
}

int main(int argc, char **argv) {
    const char *usage =
        "katka: MEM-finding with BWTs of both strands (Three Examples of MEM-Finding)\n"
        "  katka prep  [--preset NAME] [-d] [-s] [-U] [--digest k,w] list.tsv PREFIX\n"
        "  katka bwt   PREFIX\n"
        "  katka index [--preset NAME] [options] PREFIX\n"
        "  katka find  [--preset NAME] [options] PREFIX reads out\n"
        "  katka info  PREFIX\n"
        "  katka stpd-prep [-s] list.tsv OUT; katka stpd-export OUT   (input files for the published STPD-index)\n"
        "presets: zymo, silva, chr19, novel (katka presets lists their options)\n";
    if (argc < 2) { fputs(usage, stderr); return 1; }
    string cmd = argv[1];
    try {
        if (cmd == "prep") return cmd_prep(argc - 1, argv + 1);
        if (cmd == "bwt") return cmd_bwt(argc - 1, argv + 1);
        if (cmd == "index") return cmd_index(argc - 1, argv + 1);
        if (cmd == "find") return cmd_find(argc - 1, argv + 1);
        if (cmd == "info") return cmd_info(argc - 1, argv + 1);
        if (cmd == "stpd-prep") return cmd_stpd_prep(argc - 1, argv + 1);
        if (cmd == "stpd-export") return cmd_stpd_export(argc - 1, argv + 1);
        if (cmd == "presets") {
            for (auto &p : PRESETS) printf("%s: %s\n  prep:  %s\n  index: %s\n  find:  %s\n", p.name, p.use, p.prep, p.index, p.find);
            return 0;
        }
    } catch (const std::exception &e) { die(e.what()); }
    fputs(usage, stderr);
    return 1;
}
