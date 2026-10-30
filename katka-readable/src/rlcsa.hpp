// The RLCSA used for backward search (Section 5, Appendix A "Indexes"): for each character c, the list of
// the runs of c in the BWT, each as (start row, number of c's before it), sorted by start.  A backward step
// by c is a predecessor search in c's list, so its cost depends on c's number of runs, not on the alphabet,
// and it gives exactly the same intervals as rlbwt::extend.
//
// Explicit lists (the default): plain arrays of 32-bit integers (64-bit if n >= 2^32), plus, for a character
// with many runs, a jump table indexed by the top bits of the row, so the binary search covers only a short
// slice (about 8 runs).  Compressed lists: Elias-Fano (elias_fano.hpp).
//
// File (katka-csabuild): "RZCSA001", n, r, C[257], one byte (0 explicit 32-bit, 1 Elias-Fano, 2 explicit
// 64-bit), then the lists.
#pragma once
#include "elias_fano.hpp"
#include "rlbwt.hpp"
#include <cstring>

namespace katka {

template <class W> struct run_lists {
    struct Run { W start, cum; };                  // first row of the run; number of c's before it
    static constexpr uint64_t JUMP_MIN_RUNS = 64;  // characters with fewer runs get no jump table

    uint64_t n = 0;
    uint32_t sigma = 0;
    std::vector<uint64_t> off;     // the runs of c are runs[off[c] .. off[c+1]), the last one a sentinel (n, total)
    std::vector<Run> runs;
    std::vector<uint64_t> joff;    // the jump table of c is jumps[joff[c] .. joff[c+1])
    std::vector<uint8_t> shift;    // row >> shift[c] indexes c's jump table
    std::vector<uint32_t> jumps;   // jumps[b] = number of runs of c starting before row b << shift[c]
    uint64_t nruns = 0;

    // number of c's in BWT[0, i)
    uint64_t rank(uint32_t c, uint64_t i) const {
        if (c >= sigma || i == 0) return 0;
        const Run *r = runs.data() + off[c];
        uint64_t R = off[c + 1] - off[c] - 1, lo = 0, hi = R;
        if (R == 0) return 0;
        if (joff[c + 1] != joff[c]) {                  // narrow [lo, hi] with the jump table
            const uint32_t *J = jumps.data() + joff[c];
            uint64_t b = i >> shift[c];
            lo = J[b]; hi = J[b + 1];
        }
        while (lo < hi) {                              // first run starting at or after i
            uint64_t mid = (lo + hi) >> 1;
            if (r[mid].start < i) lo = mid + 1; else hi = mid;
        }
        if (lo == 0) return 0;
        const Run &a = r[lo - 1];                      // the last run of c starting before i
        uint64_t len = r[lo].cum - a.cum, d = i - a.start;
        return a.cum + (d < len ? d : len);
    }
    uint64_t bytes() const { return off.size() * 8 + runs.size() * sizeof(Run) + joff.size() * 8 + shift.size() + jumps.size() * 4; }
    void save(FILE *f) const {
        fwrite(&n, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&nruns, 8, 1, f);
        write_vec(f, off); write_vec(f, runs); write_vec(f, joff); write_vec(f, shift); write_vec(f, jumps);
    }
    void load(FILE *f) {
        if (fread(&n, 8, 1, f) != 1 || fread(&sigma, 4, 1, f) != 1 || fread(&nruns, 8, 1, f) != 1) throw std::runtime_error("read");
        read_vec(f, off); read_vec(f, runs); read_vec(f, joff); read_vec(f, shift); read_vec(f, jumps);
    }
};

struct rlcsa {
    uint64_t n = 0, R = 0, C[257];
    bool ef = false, wide = false;      // Elias-Fano lists; 64-bit explicit lists
    run_lists<uint32_t> rr;
    run_lists<uint64_t> rr64;
    EFRank er;

    // from the BWT file: one pass to count the runs of each character, one to fill the lists
    void build(const std::string &bwtfile, bool use_ef) {
        ef = use_ef;
        uint64_t count[256] = {0}, nruns[256] = {0};
        rlbwt::for_each_run(bwtfile, [&](uint64_t k, uint64_t s, uint64_t len, unsigned char c) { count[c] += len; nruns[c]++; R = k + 1; n = s + len; });
        C[0] = 0;
        for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + count[c];
        if (ef) build_ef(bwtfile, nruns);
        else if ((wide = (n >> 32) != 0)) build_explicit(rr64, bwtfile, nruns);
        else build_explicit(rr, bwtfile, nruns);
    }
    template <class RL> void build_explicit(RL &rl, const std::string &bwtfile, const uint64_t *nruns) {
        typedef decltype(rl.runs[0].start) W;
        rl.n = n; rl.sigma = 256; rl.nruns = R;
        rl.off.assign(257, 0);
        for (int c = 0; c < 256; ++c) rl.off[c + 1] = rl.off[c] + nruns[c] + 1;
        rl.runs.assign(rl.off[256], typename RL::Run{0, 0});
        std::vector<uint64_t> next(rl.off.begin(), rl.off.end() - 1), seen(256, 0);
        rlbwt::for_each_run(bwtfile, [&](uint64_t, uint64_t s, uint64_t len, unsigned char c) {
            rl.runs[next[c]++] = typename RL::Run{(W)s, (W)seen[c]}; seen[c] += len; });
        for (int c = 0; c < 256; ++c) rl.runs[next[c]] = typename RL::Run{(W)n, (W)seen[c]};    // sentinels
        // jump tables: about Rc/8 buckets for a character with Rc runs
        rl.joff.assign(257, 0); rl.shift.assign(256, 0); rl.jumps.clear();
        for (int c = 0; c < 256; ++c) {
            rl.joff[c] = rl.jumps.size();
            uint64_t Rc = rl.off[c + 1] - rl.off[c] - 1;
            if (Rc < RL::JUMP_MIN_RUNS) continue;
            uint8_t sh = 0;
            while ((n >> sh) > Rc / 8) ++sh;
            rl.shift[c] = sh;
            const typename RL::Run *r = &rl.runs[rl.off[c]];
            for (uint64_t b = 0, k = 0; b < (n >> sh) + 2; ++b) {
                while (k < Rc && r[k].start < (b << sh)) ++k;
                rl.jumps.push_back((uint32_t)k);
            }
        }
        rl.joff[256] = rl.jumps.size();
    }
    void build_ef(const std::string &bwtfile, const uint64_t *nruns) {
        er.n = n; er.sigma = 256; er.nruns = R;
        std::vector<uint64_t> slot(257, 0), seen(256, 0), ks(R), fp(R);
        for (int c = 0; c < 256; ++c) slot[c + 1] = slot[c] + nruns[c];
        rlbwt::for_each_run(bwtfile, [&](uint64_t, uint64_t s, uint64_t len, unsigned char c) {
            uint64_t k = slot[c]++; ks[k] = (uint64_t)c * n + s; fp[k] = C[c] + seen[c]; seen[c] += len; });
        er.KS.build(ks, 256 * n);
        std::vector<uint64_t>().swap(ks);
        er.FP.build(fp, n);
    }

    // one backward step: [sp, ep) for P -> [sp, ep) for cP (empty if sp >= ep)
    void extend(uint64_t &sp, uint64_t &ep, unsigned char c) const {
        uint64_t a, b;
        if (ef) er.rank2(c, sp, ep, a, b);
        else if (wide) { a = rr64.rank(c, sp); b = rr64.rank(c, ep); }
        else { a = rr.rank(c, sp); b = rr.rank(c, ep); }
        sp = C[c] + a; ep = C[c] + b;
    }

    uint64_t bytes() const { return 8 * 259 + (ef ? er.bytes() : wide ? rr64.bytes() : rr.bytes()); }
    void save(const std::string &file) const {
        FILE *f = fopen(file.c_str(), "wb");
        if (!f) die("cannot write " + file);
        fwrite("RZCSA001", 1, 8, f); fwrite(&n, 8, 1, f); fwrite(&R, 8, 1, f); fwrite(C, 8, 257, f);
        uint8_t e = ef ? 1 : wide ? 2 : 0; fwrite(&e, 1, 1, f);
        if (ef) er.save(f); else if (wide) rr64.save(f); else rr.save(f);
        fclose(f);
    }
    void load(const std::string &file) {
        FILE *f = fopen(file.c_str(), "rb");
        char magic[8]; uint8_t e = 0;
        if (!f || fread(magic, 1, 8, f) != 8 || memcmp(magic, "RZCSA001", 8) || fread(&n, 8, 1, f) != 1 || fread(&R, 8, 1, f) != 1
            || fread(C, 8, 257, f) != 257 || fread(&e, 1, 1, f) != 1) die("bad RLCSA file " + file);
        ef = e == 1; wide = e == 2;
        if (ef) er.load(f); else if (wide) rr64.load(f); else rr.load(f);
        fclose(f);
    }
};

}  // namespace katka
