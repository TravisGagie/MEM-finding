// katka-prep: build the text S from a list of FASTA files (Appendix A, "Reference text"):
//     S = G_1 X rc(G_1) X G_2 X rc(G_2) X ... G_k X rc(G_k) X
// where G_i is genome (or genus) i, its records joined by X, and every character other than A, C, G, T (in
// either case) becomes X.  Having both strands lets one BWT serve backward and forward search (Section 5).
//
// usage: katka-prep [-d] [-s] [-U] list.tsv out-prefix
//   list.tsv : one line per FASTA file (plain or gzipped), in the order the genomes should have in S
//              (tree order for SILVA):  <tag>\t<fasta>   (a line with one field is a path, tag "-")
//   -s : a FASTA file may hold many genomes; a genome is a maximal run of consecutive records with the same
//        sample ID (the header up to the first '.' or blank)
//   -d : also write one file per dataset (a maximal run of consecutive genomes with the same tag), as
//        out.ds<j>.S, listed in out.datasets, to build the BWT dataset by dataset and merge (pfp-merge)
//   -U : homopolymer compression: collapse every run of equal characters of each genome to one character
//        (before taking the reverse complement, with which it commutes); classify with katka-classify -U
// writes
//   out.S    the text (bytes A, C, G, T, X; no terminator)
//   out.tbl  one line per genome: index, tag, file (#sample ID with -s), length, start in S
#include "common.hpp"
#include <zlib.h>
#include <algorithm>
#include <fstream>
#include <vector>

using namespace katka;
using std::string;

static char clean(char c) {
    switch (c) { case 'A': case 'a': return 'A'; case 'C': case 'c': return 'C'; case 'G': case 'g': return 'G'; case 'T': case 't': return 'T'; default: return 'X'; }
}

// The genomes of a FASTA file, as (sample ID, sequence).  Without split, the whole file is one genome
// (with an empty ID); with it, a new genome starts whenever a header's sample ID changes.
static std::vector<std::pair<string, string>> read_genomes(const string &path, bool split) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) die("cannot open " + path);
    std::vector<std::pair<string, string>> out;
    std::vector<char> buf(1 << 20);
    bool header = false, bol = true;
    string hdr;
    int got;
    auto new_record = [&]() {                       // a header has just been read into hdr
        string id = split ? hdr.substr(0, hdr.find_first_of(". \t")) : "";
        if (out.empty() || (split && out.back().first != id)) out.push_back({id, string()});
        else if (!out.back().second.empty() || split) out.back().second.push_back('X');   // separate records
    };
    while ((got = gzread(f, buf.data(), buf.size())) > 0) {
        for (int i = 0; i < got; ++i) {
            char c = buf[i];
            if (bol && c == '>') { header = true; hdr.clear(); bol = false; continue; }
            if (c == '\n') { if (header) new_record(); header = false; bol = true; continue; }
            bol = false;
            if (header) { hdr.push_back(c); continue; }
            if (c == '\r' || c == ' ' || c == '\t') continue;
            if (out.empty()) out.push_back({split ? "-" : "", string()});
            out.back().second.push_back(clean(c));
        }
    }
    gzclose(f);
    if (out.empty() && !split) out.push_back({"", string()});
    return out;
}

static string reverse_complement(const string &g) {
    string rc(g.rbegin(), g.rend());
    for (auto &c : rc) c = (char)complement((unsigned char)c);
    return rc;
}

int main(int argc, char **argv) {
    bool datasets = false, split = false, hpc = false;
    int a = 1;
    for (; a < argc && argv[a][0] == '-'; ++a) {
        string o = argv[a];
        if (o == "-d") datasets = true; else if (o == "-s") split = true; else if (o == "-U") hpc = true;
        else die("unknown option " + o + " (digests are not in this version)");
    }
    if (argc - a != 2) die("usage: katka-prep [-d] [-s] [-U] list.tsv out-prefix");
    std::ifstream list(argv[a]);
    if (!list) die(string("cannot open ") + argv[a]);
    const string out = argv[a + 1];
    FILE *fs = fopen((out + ".S").c_str(), "wb");
    std::ofstream tbl(out + ".tbl"), dslist;
    if (datasets) dslist.open(out + ".datasets");
    FILE *fd = nullptr;                             // the current dataset file
    string last_tag, line;
    u64 pos = 0, idx = 0, raw = 0, nds = 0;
    while (std::getline(list, line)) {
        if (line.empty() || line[0] == '#') continue;
        string tag = "-", path = line;
        size_t tab = line.find('\t');
        if (tab != string::npos) { tag = line.substr(0, tab); path = line.substr(tab + 1); }
        for (auto &[id, g] : read_genomes(path, split)) {
            string name = split ? path + "#" + id : path;
            if (g.empty()) { fprintf(stderr, "warning: empty genome %s\n", name.c_str()); continue; }
            if (hpc) { raw += g.size(); g.erase(std::unique(g.begin(), g.end()), g.end()); }
            string rc = reverse_complement(g);
            tbl << idx << "\t" << tag << "\t" << name << "\t" << g.size() << "\t" << pos << "\n";
            fwrite(g.data(), 1, g.size(), fs); fputc('X', fs);
            fwrite(rc.data(), 1, rc.size(), fs); fputc('X', fs);
            if (datasets) {
                if (!fd || tag != last_tag) {
                    if (fd) fclose(fd);
                    string ds = out + ".ds" + std::to_string(nds++) + ".S";
                    fd = fopen(ds.c_str(), "wb");
                    dslist << ds << "\n";
                    last_tag = tag;
                }
                fwrite(g.data(), 1, g.size(), fd); fputc('X', fd);
                fwrite(rc.data(), 1, rc.size(), fd); fputc('X', fd);
            }
            pos += 2 * (g.size() + 1);
            ++idx;
        }
    }
    fclose(fs);
    if (fd) fclose(fd);
    fprintf(stderr, "genomes=%lu |S|=%lu", idx, pos);
    if (datasets) fprintf(stderr, " datasets=%lu", nds);
    if (hpc) fprintf(stderr, " homopolymer-compressed: %lu bases -> %lu per strand (%g)", raw, pos / 2, (double)pos / 2 / (raw ? raw : 1));
    fprintf(stderr, "\n");
    return 0;
}
