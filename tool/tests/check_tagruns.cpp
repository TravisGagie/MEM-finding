// Independent check of a KTRUNS01 tag-runs file against the BWT. Walk LF back from each dataset's terminator
// Build: g++ -std=c++17 -O3 -march=native -Isrc tests/check_tagruns.cpp -o bin/check_tagruns -lz -pthread   (from tool/)
// for STEPS steps (as sample_rindex does); each visited row holds a known text position, hence a known genome
// and tag, which must equal the tag of the run containing that row in the file. The move structure is freed
// before the file is read, and the file is read through mmap, so the peak is about the move structure's size.
//   check_tagruns PREFIX TAGRUNS STEPS [by_tag=1]
#include "tag_build.hpp"
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
using namespace katka;
static void load_rlbwt(rlbwt &B, const std::string &file) { FILE *f = open_or_die(file, "rb"); rmagic(f, "KTRLBW01", file); B.load(f); fclose(f); }
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: check_tagruns PREFIX TAGRUNS STEPS [by_tag]\n"); return 1; }
    const std::string prefix = argv[1], Sfile = prefix + ".S", bwtfile = prefix + ".bwt";
    const u64 STEPS = std::stoull(argv[3]); const bool by_tag = argc < 5 || std::atoi(argv[4]);
    genome_table G; G.load(prefix + ".tbl");
    std::vector<uint32_t> gt = G.tags(by_tag);
    std::vector<std::pair<u64, uint32_t>> want;           // (row, tag from the text)
    u64 n, K, R;
    {
        rlbwt B; load_rlbwt(B, prefix + ".rlbwt");
        n = B.n; R = B.R; K = B.C[2] - B.C[1];
        const u64 N = file_size(Sfile);
        std::vector<u64> doff = K == 1 ? std::vector<u64>{0} : G.dataset_starts(); doff.push_back(N);
        move_lf M; M.build(bwtfile, B);
        std::vector<u64> drow = terminator_rows(M, doff, K, Sfile);
        auto text_tag = [&](u64 q) { return gt[std::upper_bound(G.start.begin(), G.start.end(), q) - G.start.begin() - 1]; };
        for (u64 j = 0; j < K; ++j) {
            u64 len = doff[j + 1] - doff[j], row = drow[j], k = M.run_containing(row);
            u64 steps = std::min(STEPS, len + 1);
            for (u64 s = 0; s < steps; ++s) {             // row holds the suffix at S position doff[j] + (len - s)
                u64 p = len - s, q = p == len ? doff[j + 1] - 1 : doff[j] + p;   // the terminator: its dataset's last character
                want.push_back({row, text_tag(q)});
                if (p) M.LF(row, k);
            }
        }
    }
    int fd = open(argv[2], O_RDONLY); if (fd < 0) die(std::string("cannot open ") + argv[2]);
    const u64 fsz = file_size(argv[2]);
    const char *m = (const char *)mmap(nullptr, fsz, PROT_READ, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) die("mmap");
    if (memcmp(m, "KTRUNS01", 8)) die("not a KTRUNS01 file");
    const u64 fn = *(const u64 *)(m + 8), r = *(const u64 *)(m + 16);
    const u64 *st = (const u64 *)(m + 24);
    const u64 r2 = *(const u64 *)(m + 24 + 8 * r);
    const uint32_t *tg = (const uint32_t *)(m + 32 + 8 * r);
    printf("BWT n = %lu, R = %lu, datasets K = %lu; file n = %lu, runs = %lu\n", n, R, K, fn, r);
    if (fn != n) die("n differs");
    if (r2 != r || fsz != 32 + 12 * r) die("bad file layout");
    std::sort(want.begin(), want.end());
    u64 bad = 0, rows = 0, last = ~0ULL; std::vector<char> seen_run;   // distinct runs touched
    u64 runs_touched = 0, prev_run = ~0ULL;
    for (auto &[row, t] : want) {
        if (row == last) continue; last = row; ++rows;
        u64 i = std::upper_bound(st, st + r, row) - st - 1;
        if (i != prev_run) { ++runs_touched; prev_run = i; }
        if (tg[i] != t && bad++ < 10) printf("  MISMATCH row %lu (run %lu): file %u, text %u\n", row, i, tg[i], t);
    }
    printf("checked %lu distinct rows (%lu runs) from %lu datasets: %lu mismatches\n", rows, runs_touched, K, bad);
    return bad != 0;
}
