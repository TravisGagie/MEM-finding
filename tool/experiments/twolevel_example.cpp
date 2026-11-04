// Two-level search with blocks, used from code (tool/CHR19.md, "Using two-level search with blocks from
// code").  Loads an index built with `katka index --preset chr19 PREFIX` (and, for locating, `katka index
// --twolevel-locate PREFIX`), then for each pattern of a FASTA file:
//   1. one longest-suffix query: the longest suffix of the pattern that occurs, with its BWT interval
//      (search2, the primitive every MEM-finding algorithm here is built on);
//   2. its MEMs of length at least L by BML (find_mems over a segment_search with the two-level index),
//      each with the BWT interval of its reverse complement and, with the .tloc file, its occurrences.
//
// Build, from tool/:
//   g++ -std=c++17 -O3 -march=native -Isrc -Ithird_party/libsais/include experiments/twolevel_example.cpp \
//       build/libsais.o build/libsais64.o -o bin/twolevel-example -lz -pthread
// Run:
//   bin/twolevel-example PREFIX patterns.fa L [--locate]
#include "mems.hpp"
#include "tloc_build.hpp"
#include <fstream>
using namespace katka;

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: twolevel-example PREFIX patterns.fa L [--locate]\n"); return 1; }
    const std::string prefix = argv[1];
    const u64 L = std::stoull(argv[3]);
    const bool locate = argc > 4 && std::string(argv[4]) == "--locate";

    // the index: the RLBWT (needed for the toehold when locating), the RLCSA (character steps), the two-level
    // index with blocks
    rlbwt B;
    { FILE *f = open_or_die(prefix + ".rlbwt", "rb"); rmagic(f, "KTRLBW01", prefix + ".rlbwt"); B.load(f); fclose(f); }
    rlcsa X; X.load(prefix + ".csa");
    two_level TL; TL.load(prefix + ".lvl");
    twolevel_locate TLC;
    if (locate) TLC.load(prefix + ".tloc");
    stepper S; S.bwt = &B; S.csa = &X;

    segment_search Q;                    // one per thread: it holds the pattern's parse and scratch space
    Q.S = &S; Q.TL = &TL; Q.jumps = true;
    if (locate) { Q.TLC = &TLC; Q.locate = true; }

    std::ifstream in(argv[2]);
    std::string line, name, seq;
    auto run = [&]() {
        if (seq.empty()) return;
        for (auto &c : seq) c = (char)toupper((unsigned char)c);
        // 1. the longest suffix of the first A/C/G/T run that occurs (a two-level search with block steps)
        u64 r0 = 0; while (r0 < seq.size() && base_code(seq[r0]) < 0) ++r0;
        u64 r1 = r0; while (r1 < seq.size() && base_code(seq[r1]) >= 0) ++r1;
        if (r1 > r0) {
            Q.prepare(seq.data() + r0, r1 - r0);                  // parses the run (and its reverse complement)
            interval I;
            u64 j = Q.back(0, 0, r1 - r0, I);                     // X[j..end) is the longest suffix that occurs
            printf("%s: longest occurring suffix of its first run starts at %lu (%lu characters), %lu occurrences\n",
                   name.c_str(), r0 + j, r1 - r0 - j, I.ep - I.sp);
        }
        // 2. the MEMs of length >= L, by BML (find_mems splits the pattern at characters other than A/C/G/T)
        mem_options o; o.L = L;
        u64 nm = 0;
        find_mems(Q, seq, o, [&](const mem &M) {
            ++nm;
            printf("%s\t%lu\t%lu\t%lu\t%lu", name.c_str(), M.a, M.e, M.I.sp, M.I.ep);   // [a, e) and rc's interval
            if (locate) {                                         // Q.occbuf: the positions of rc(MEM) in S
                std::vector<u64> p = Q.occbuf;
                std::sort(p.begin(), p.end());
                printf("\t%zu occurrences, first at %lu", p.size(), p.empty() ? 0 : p[0]);
            }
            printf("\n");
        });
        fprintf(stderr, "%s: %lu MEMs of length >= %lu\n", name.c_str(), nm, L);
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') { run(); name = line.substr(1, line.find_first_of(" \t") - 1); seq.clear(); }
        else seq += line;
    }
    run();
    return 0;
}
