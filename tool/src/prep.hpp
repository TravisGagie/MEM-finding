// The text S (KATKA App. A, "Reference text"; Section 2 of the paper, "Both strands"):
//     S = G_1 X rc(G_1) X G_2 X rc(G_2) X ... G_k X rc(G_k) X
// where G_i is genome (or genus) i, its records joined by X, and every character other than A, C, G, T
// (in either case) becomes X.  Having both strands lets one BWT serve backward and forward search.
//
//   list.tsv : one line per FASTA file (plain or gzipped), in the order the genomes should have in S
//              (tree order for SILVA):  <tag>\t<fasta>   (a line with one field is a path, tag "-")
//   split    : a FASTA file may hold many genomes, a genome being a maximal run of consecutive records
//              with the same sample ID (the header up to the first '.' or blank)
//   datasets : also write one file per dataset (a maximal run of consecutive genomes with the same tag),
//              prefix.ds<j>.S, listed in prefix.datasets, to build the BWT dataset by dataset (pfp-merge)
//   hpc      : homopolymer compression: collapse every run of equal characters of each genome (before
//              taking the reverse complement, with which it commutes)
// writes prefix.S (no terminator), prefix.tbl (index, tag, name, length, start in S) and prefix.info.
// This is KATKA's katka-prep, unchanged except for prefix.info.
#pragma once
#include "common.hpp"
#include <algorithm>
#include <fstream>
#include <zlib.h>

namespace katka {

struct prep_options { bool datasets = false, split = false, hpc = false; };

inline char clean_base(char c) {
    switch (c) { case 'A': case 'a': return 'A'; case 'C': case 'c': return 'C'; case 'G': case 'g': return 'G'; case 'T': case 't': return 'T'; default: return 'X'; }
}

inline std::vector<std::pair<std::string, std::string>> read_genomes(const std::string &path, bool split) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) die("cannot open " + path);
    std::vector<std::pair<std::string, std::string>> out;
    std::vector<char> buf(1 << 20);
    bool header = false, bol = true;
    std::string hdr;
    int got;
    auto new_record = [&]() {
        std::string id = split ? hdr.substr(0, hdr.find_first_of(". \t")) : "";
        if (out.empty() || (split && out.back().first != id)) out.push_back({id, std::string()});
        else if (!out.back().second.empty() || split) out.back().second.push_back('X');
    };
    while ((got = gzread(f, buf.data(), buf.size())) > 0) {
        for (int i = 0; i < got; ++i) {
            char c = buf[i];
            if (bol && c == '>') { header = true; hdr.clear(); bol = false; continue; }
            if (c == '\n') { if (header) new_record(); header = false; bol = true; continue; }
            bol = false;
            if (header) { hdr.push_back(c); continue; }
            if (c == '\r' || c == ' ' || c == '\t') continue;
            if (out.empty()) out.push_back({split ? "-" : "", std::string()});
            out.back().second.push_back(clean_base(c));
        }
    }
    gzclose(f);
    if (out.empty() && !split) out.push_back({"", std::string()});
    return out;
}

inline void write_info(const std::string &prefix, bool hpc) {
    FILE *f = open_or_die(prefix + ".info", "w");
    fprintf(f, "hpc\t%d\n", hpc ? 1 : 0);
    fclose(f);
}
inline bool read_info_hpc(const std::string &prefix) {
    FILE *f = fopen((prefix + ".info").c_str(), "r");
    if (!f) return false;
    int h = 0;
    if (fscanf(f, "hpc\t%d", &h) != 1) h = 0;
    fclose(f);
    return h;
}

inline void prep(const std::string &listfile, const std::string &out, const prep_options &o) {
    std::ifstream list(listfile);
    if (!list) die("cannot open " + listfile);
    FILE *fs = open_or_die(out + ".S", "wb");
    std::ofstream tbl(out + ".tbl"), dslist;
    if (o.datasets) dslist.open(out + ".datasets");
    FILE *fd = nullptr;
    std::string last_tag, line;
    u64 pos = 0, idx = 0, raw = 0, nds = 0;
    while (std::getline(list, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::string tag = "-", path = line;
        size_t tab = line.find('\t');
        if (tab != std::string::npos) { tag = line.substr(0, tab); path = line.substr(tab + 1); }
        for (auto &[id, g] : read_genomes(path, o.split)) {
            std::string name = o.split ? path + "#" + id : path;
            if (g.empty()) { fprintf(stderr, "warning: empty genome %s\n", name.c_str()); continue; }
            if (o.hpc) { raw += g.size(); g.erase(std::unique(g.begin(), g.end()), g.end()); }
            std::string rc = reverse_complement(g);
            tbl << idx << "\t" << tag << "\t" << name << "\t" << g.size() << "\t" << pos << "\n";
            fwrite(g.data(), 1, g.size(), fs); fputc('X', fs);
            fwrite(rc.data(), 1, rc.size(), fs); fputc('X', fs);
            if (o.datasets) {
                if (!fd || tag != last_tag) {
                    if (fd) fclose(fd);
                    std::string ds = out + ".ds" + std::to_string(nds++) + ".S";
                    fd = open_or_die(ds, "wb");
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
    write_info(out, o.hpc);
    fprintf(stderr, "genomes=%lu |S|=%lu", idx, pos);
    if (o.datasets) fprintf(stderr, " datasets=%lu", nds);
    if (o.hpc) fprintf(stderr, " homopolymer-compressed: %lu bases -> %lu per strand (%g)", raw, pos / 2, (double)pos / 2 / (raw ? raw : 1));
    fprintf(stderr, "\n");
}

}  // namespace katka
