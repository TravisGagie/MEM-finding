// katka-classify: find the long MEMs of each read with BML (bml.hpp) and answer each one from the tag array
// (tag_array.hpp).  The scoring (which credit scheme, which taxon wins) is done afterwards by the scripts in
// scripts/, from the listings this program writes.
//
// usage: katka-classify [options] -L bases -T index.tag index.rlbwt reads.fq out.listings
//   index.rlbwt   the RLBWT (katka-build -G); also needed for LF with a sampled tag array
//   reads         FASTQ or FASTA (several lines per record), plain text
//   -L bases      report the MEMs of at least that many bases (with -X: the starting value)
//   -T file       the tag array (katka-tagbuild)
//   -c            count: list each tag of each MEM with the number of occurrences of the MEM in it (for
//                 proportional credit)
//   -l            list the distinct tags of each MEM (even credit)
//                 Without -c or -l: the smallest and largest tags of each MEM (the range scheme)
//   -C file       backward search on the RLCSA (katka-csabuild) instead of the RLBWT
//   -F t          lookup table of all strings of t <= 13 bases (built at start-up)
//   -X            longest MEMs: adaptive L (bml.hpp)
//   -K file       KeBaB filter (katka-kbbuild); its k must not exceed L
//   -D d          ignore MEMs with fewer than d distinct 3-mers
//   -U            homopolymer compression: the index was built with katka-prep -U; each read's runs of equal
//                 bases are collapsed before BML (so L and -D apply to the compressed read), and the listings
//                 give each MEM's span in the original read
//   -j p          p threads, sharing the index; the output is the same as with one thread
// Output (Cliffy's listing format, extended with counts), one record per read:
//   >name
//   [a,b] {t:c,t:c,...} [a,b] {...} ...
// with [a,b] the MEM's span in the read (0-based, inclusive) and t the tags (0-based, in order of first
// appearance in the .tbl with katka-tagbuild -S, else the genome index).
#include "bml.hpp"
#include "tag_array.hpp"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <unistd.h>

using namespace katka;
using std::string;

enum answer_mode { RANGE, LIST, COUNT };

struct classifier {
    searcher S;
    tag_array *tags = nullptr;
    bml_options opt;
    answer_mode mode = RANGE;
    bool hpc = false;

    // the answer for a MEM whose reverse complement has BWT interval I, appended to out
    void answer(const interval &I, string &out) const {
        char b[48];
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

    // the listing line of one read; returns its number of MEMs
    u64 classify(const string &read, string &line) const {
        line.clear();
        static thread_local string h;                    // the homopolymer-compressed read
        static thread_local std::vector<u64> pos;        // pos[i] = position in the read of h[i]; pos[|h|] = |read|
        if (hpc) {
            h.clear(); pos.clear();
            for (u64 i = 0; i < read.size(); ++i) if (i == 0 || read[i] != read[i - 1]) { h.push_back(read[i]); pos.push_back(i); }
            pos.push_back(read.size());
        }
        u64 nmems = 0;
        char b[64];
        find_mems(S, hpc ? h : read, opt, [&](u64 a, u64 e, const interval &I) {
            u64 first = hpc ? pos[a] : a, last = hpc ? pos[e] - 1 : e - 1;   // the span of the MEM in the read
            snprintf(b, sizeof b, "[%lu,%lu] {", first, last);
            line += b; answer(I, line); line += "} ";
            ++nmems;
        });
        return nmems;
    }
};

// reads FASTQ (four lines per record) or FASTA (a header, then one or more sequence lines)
struct read_parser {
    std::ifstream in;
    bool fasta = false;
    string pending;                                     // FASTA: the next record's header, already read
    explicit read_parser(const string &f) : in(f) { fasta = in.peek() == '>'; }
    // the read's name (the header up to the first blank) and its sequence in upper case
    bool next(string &name, string &seq) {
        string h, line;
        if (!fasta) {
            string plus, qual;
            if (!(std::getline(in, h) && std::getline(in, seq) && std::getline(in, plus) && std::getline(in, qual))) return false;
        } else {
            if (pending.empty()) {
                while (std::getline(in, line)) if (!line.empty() && line[0] == '>') { pending = line; break; }
                if (pending.empty()) return false;
            }
            h.swap(pending); pending.clear(); seq.clear();
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '>') { pending = line; break; }
                if (!line.empty() && line.back() == '\r') line.pop_back();
                seq += line;
            }
        }
        for (auto &c : seq) c = (char)toupper((unsigned char)c);
        size_t blank = h.find_first_of(" \t");
        name = h.substr(1, blank == string::npos ? string::npos : blank - 1);
        return true;
    }
};

// Workers take chunks of reads from the input in turn and write their chunks in the same order.
struct run_stats { u64 reads = 0, mems = 0, unclassified = 0; };
static run_stats classify_all(const classifier &K, read_parser &P, FILE *out, int nthreads) {
    const size_t CHUNK = 256;
    std::mutex in_m, out_m, st_m;
    std::condition_variable out_cv;
    u64 next_in = 0, next_out = 0;
    bool eof = false;
    run_stats total;
    auto work = [&]() {
        std::vector<string> names, seqs, lines;
        run_stats st;
        string name, seq;
        for (;;) {
            u64 my;
            {
                std::lock_guard<std::mutex> g(in_m);
                names.clear(); seqs.clear();
                while (!eof && names.size() < CHUNK) {
                    if (!P.next(name, seq)) { eof = true; break; }
                    names.push_back(name); seqs.push_back(seq);
                }
                my = next_in++;
            }
            lines.resize(names.size());
            for (size_t i = 0; i < names.size(); ++i) {
                u64 m = K.classify(seqs[i], lines[i]);
                ++st.reads; st.mems += m; st.unclassified += m == 0;
            }
            {
                std::unique_lock<std::mutex> lk(out_m);
                out_cv.wait(lk, [&] { return next_out == my; });
                for (size_t i = 0; i < names.size(); ++i) fprintf(out, ">%s\n%s\n", names[i].c_str(), lines[i].c_str());
                ++next_out;
            }
            out_cv.notify_all();
            if (names.empty()) break;
        }
        std::lock_guard<std::mutex> g(st_m);
        total.reads += st.reads; total.mems += st.mems; total.unclassified += st.unclassified;
    };
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back(work);
    for (auto &x : th) x.join();
    return total;
}

int main(int argc, char **argv) {
    const char *usage = "usage: katka-classify [-c | -l] [-X] [-U] [-C csa] [-F t] [-K kbf] [-D d] [-j threads] -L bases -T tags index.rlbwt reads out.listings\n";
    string tagfile, csafile, kbffile;
    u64 L = 0; int nthreads = 1, ft = 0, opt; bool list = false, counts = false, longest = false, hpc = false;
    int min3 = 0;
    while ((opt = getopt(argc, argv, "clXUC:F:K:D:j:L:T:")) != -1) {
        switch (opt) {
            case 'c': counts = true; break;
            case 'l': list = true; break;
            case 'X': longest = true; break;
            case 'U': hpc = true; break;
            case 'C': csafile = optarg; break;
            case 'F': ft = std::stoi(optarg); break;
            case 'K': kbffile = optarg; break;
            case 'D': min3 = std::stoi(optarg); break;
            case 'j': nthreads = std::max(1, std::stoi(optarg)); break;
            case 'L': L = std::stoull(optarg); break;
            case 'T': tagfile = optarg; break;
            default: fputs(usage, stderr); fputs("(the original code's other options are not in this version)\n", stderr); return 1;
        }
    }
    if (argc - optind != 3 || L == 0 || tagfile.empty()) { fputs(usage, stderr); return 1; }
    if (ft < 0 || ft > 13) die("-F needs t <= 13");

    rlbwt_file Z; Z.load(argv[optind]);
    rlcsa X;
    if (!csafile.empty()) X.load(csafile);
    tag_array T;
    T.load(tagfile, !counts);                            // counting needs no RMQ
    T.bwt = &Z.bwt;
    kmer_filter KF;
    if (!kbffile.empty()) {
        if (!KF.load(kbffile)) die("cannot load " + kbffile);
        if (KF.k > L) die("-K: the filter's k = " + std::to_string(KF.k) + " exceeds L");
        fprintf(stderr, "KeBaB filter %s: k=%u, h=%u, %.3f GB\n", kbffile.c_str(), KF.k, KF.h, KF.bytes() / 1e9);
    }

    classifier K;
    K.S.bwt = &Z.bwt; K.S.csa = csafile.empty() ? nullptr : &X;
    K.tags = &T;
    K.mode = counts ? COUNT : list ? LIST : RANGE;
    K.hpc = hpc;
    K.opt.L = L; K.opt.longest = longest; K.opt.min3 = min3;
    K.opt.filter = kbffile.empty() ? nullptr : &KF;
    lookup_table LT;
    if (ft) {
        double t0 = now(); LT.build(K.S, ft);
        fprintf(stderr, "lookup table: t=%d, %.1f MB, built in %.1f s\n", ft, LT.bytes() / 1e6, now() - t0);
        K.opt.table = &LT;
    }

    read_parser P(argv[optind + 1]);
    FILE *out = fopen(argv[optind + 2], "w");
    if (!P.in || !out) die("cannot open the reads or the output");
    double t0 = now();
    run_stats st = classify_all(K, P, out, nthreads);
    double wall = now() - t0;
    fclose(out);
    fprintf(stderr, "reads=%lu MEMs=%lu unclassified=%lu | threads=%d wall %.2f s, %.0f reads/s, %.1f us/read\n",
            st.reads, st.mems, st.unclassified, nthreads, wall, wall > 0 ? st.reads / wall : 0.0, 1e6 * wall * nthreads / std::max<u64>(1, st.reads));
    return 0;
}
