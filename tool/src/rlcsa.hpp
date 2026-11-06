// The RLCSA used for backward search (Section 2, "RLBWTs and RLCSAs"): for each character c, the list of
// the runs of c in the BWT, each as (start row, number of c's before it), sorted by start.  A backward step
// by c is a predecessor search in c's list, so its cost depends on c's number of runs, not on the alphabet,
// and it gives exactly the same intervals as rlbwt::extend.
//
//   explicit       plain arrays of 32-bit integers (64-bit if n >= 2^32), plus, for a character with many
//                  runs, a jump table indexed by the top bits of the row, so the binary search covers only
//                  about 8 runs (as in KATKA);
//   Elias-Fano     two Elias-Fano lists over all the runs, the run starts shifted by c * n so that one list
//                  serves every character (EFRank in bits.hpp, as in More RLBWT Tricks).
#pragma once
#include "rlbwt.hpp"

namespace katka {

// probes of the binary searches over the runs, and the steps, for a build with -DKATKA_COUNT_PROBES (single-threaded runs)
#ifdef KATKA_COUNT_PROBES
inline uint64_t g_probes = 0, g_steps = 0;
#define KATKA_PROBE() (++g_probes)
#define KATKA_STEP() (++g_steps)
#else
#define KATKA_PROBE() ((void)0)
#define KATKA_STEP() ((void)0)
#endif

template <class W> struct run_lists {
    struct Run { W start, cum; };                  // first row of the run; number of c's before it
    static constexpr uint64_t JUMP_MIN_RUNS = 64;  // characters with fewer runs get no jump table

    uint64_t n = 0;
    std::vector<uint64_t> off;     // the runs of c are runs[off[c] .. off[c+1]), the last one a sentinel (n, total)
    std::vector<Run> runs;
    std::vector<uint64_t> joff;    // the jump table of c is jumps[joff[c] .. joff[c+1])
    std::vector<uint8_t> shift;    // row >> shift[c] indexes c's jump table
    std::vector<uint32_t> jumps;   // jumps[b] = number of runs of c starting before row b << shift[c]

    // the first of c's runs starting at or after i, within runs lo..hi (the answer being known to lie there)
    uint64_t find_run(const Run *r, uint64_t i, uint64_t lo, uint64_t hi) const {
        while (lo < hi) {
            KATKA_PROBE();
            uint64_t mid = (lo + hi) >> 1;
            if (r[mid].start < i) lo = mid + 1; else hi = mid;
        }
        return lo;
    }
    // number of c's in BWT[0, i), from the first of c's runs starting at or after i
    uint64_t rank_from(const Run *r, uint64_t lo, uint64_t i) const {
        if (lo == 0) return 0;
        const Run &a = r[lo - 1];
        uint64_t len = r[lo].cum - a.cum, d = i - a.start;
        return a.cum + (d < len ? d : len);
    }
    // if row i of the BWT holds c: true, with the number of c's before it
    bool lf_in(uint32_t c, uint64_t i, uint64_t &o) const {
        const uint64_t lo = first_run(c, i + 1);
        if (lo == 0) return false;
        const Run *r = runs.data() + off[c];
        const Run &a = r[lo - 1];
        const uint64_t d = i - a.start;
        if (d >= r[lo].cum - a.cum) return false;
        o = a.cum + d;
        return true;
    }
    uint64_t first_run(uint32_t c, uint64_t i) const {
        const Run *r = runs.data() + off[c];
        uint64_t R = off[c + 1] - off[c] - 1, lo = 0, hi = R;
        if (R && joff[c + 1] != joff[c]) {
            const uint32_t *J = jumps.data() + joff[c];
            uint64_t b = i >> shift[c];
            lo = J[b]; hi = J[b + 1];
        }
        return find_run(r, i, lo, hi);
    }
    // number of c's in BWT[0, i) (written out in full: going through first_run and rank_from costs 3.5% on the default search)
    uint64_t rank(uint32_t c, uint64_t i) const {
        if (i == 0) return 0;
        const Run *r = runs.data() + off[c];
        uint64_t R = off[c + 1] - off[c] - 1, lo = 0, hi = R;
        if (R == 0) return 0;
        if (joff[c + 1] != joff[c]) {
            const uint32_t *J = jumps.data() + joff[c];
            uint64_t b = i >> shift[c];
            lo = J[b]; hi = J[b + 1];
        }
        while (lo < hi) {                              // first run starting at or after i
            KATKA_PROBE();
            uint64_t mid = (lo + hi) >> 1;
            if (r[mid].start < i) lo = mid + 1; else hi = mid;
        }
        if (lo == 0) return 0;
        const Run &a = r[lo - 1];
        uint64_t len = r[lo].cum - a.cum, d = i - a.start;
        return a.cum + (d < len ? d : len);
    }
    // the same, knowing that the first of c's runs at or after i is one of runs rl..rh; the jump table narrows it further
    // when that range is wider than the 8 or so runs of a jump-table entry
    uint64_t rank_in(uint32_t c, uint64_t i, uint64_t rl, uint64_t rh) const {
        if (i == 0 || off[c + 1] - off[c] == 1) return 0;
        if (rh - rl > 8 && joff[c + 1] != joff[c]) {
            const uint32_t *J = jumps.data() + joff[c];
            uint64_t b = i >> shift[c];
            if (J[b] > rl) rl = J[b];
            if (J[b + 1] < rh) rh = J[b + 1];
        }
        const Run *r = runs.data() + off[c];
        return rank_from(r, find_run(r, i, rl, rh), i);
    }
    uint64_t bytes() const { return off.size() * 8 + runs.size() * sizeof(Run) + joff.size() * 8 + shift.size() + jumps.size() * 4; }
    void save(FILE *f) const { wr(f, n); wvec(f, off); wvec(f, runs); wvec(f, joff); wvec(f, shift); wvec(f, jumps); }
    void load(FILE *f) { rd(f, n); rvec(f, off); rvec(f, runs); rvec(f, joff); rvec(f, shift); rvec(f, jumps); }
};

struct rlcsa {
    u64 n = 0, R = 0, C[257];
    bool ef = false, wide = false;      // Elias-Fano lists; 64-bit explicit lists
    run_lists<uint32_t> rr;
    run_lists<uint64_t> rr64;
    EFRank er;

    void build(const std::string &bwtfile, bool use_ef) {
        ef = use_ef;
        uint64_t count[256] = {0}, nruns[256] = {0};
        for_each_bwt_run(bwtfile, [&](uint64_t k, uint64_t s, uint64_t len, unsigned char c) { count[c] += len; nruns[c]++; R = k + 1; n = s + len; });
        C[0] = 0;
        for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + count[c];
        if (ef) build_ef(bwtfile, nruns);
        else if ((wide = (n >> 32) != 0)) build_explicit(rr64, bwtfile, nruns);
        else build_explicit(rr, bwtfile, nruns);
        prepare_lf();
    }
    template <class RL> void build_explicit(RL &rl, const std::string &bwtfile, const uint64_t *nruns) {
        typedef decltype(rl.runs[0].start) W;
        rl.n = n;
        rl.off.assign(257, 0);
        for (int c = 0; c < 256; ++c) rl.off[c + 1] = rl.off[c] + nruns[c] + 1;
        rl.runs.assign(rl.off[256], typename RL::Run{0, 0});
        std::vector<uint64_t> next(rl.off.begin(), rl.off.end() - 1), seen(256, 0);
        for_each_bwt_run(bwtfile, [&](uint64_t, uint64_t s, uint64_t len, unsigned char c) {
            rl.runs[next[c]++] = typename RL::Run{(W)s, (W)seen[c]}; seen[c] += len; });
        for (int c = 0; c < 256; ++c) rl.runs[next[c]] = typename RL::Run{(W)n, (W)seen[c]};
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
        for_each_bwt_run(bwtfile, [&](uint64_t, uint64_t s, uint64_t len, unsigned char c) {
            uint64_t k = slot[c]++; ks[k] = (uint64_t)c * n + s; fp[k] = C[c] + seen[c]; seen[c] += len; });
        er.KS.build(ks, 256 * n);
        std::vector<uint64_t>().swap(ks);
        er.FP.build(fp, n);
    }

    // the first of c's runs starting at or after row i, for extend_in
    uint64_t first_run(unsigned char c, uint64_t i) const {
        if (ef) return er.KS.count_less((uint64_t)c * n + i);
        return wide ? rr64.first_run(c, i) : rr.first_run(c, i);
    }
    // LF(i) = C[c] + rank_c(i), with no BWT: c is the character whose runs hold row i, tried in decreasing order of frequency (for
    // sampled tags: tag_array::tag)
    std::vector<uint8_t> lfc;
    void prepare_lf() {
        lfc.clear();
        for (int c = 0; c < 256; ++c) if (C[c + 1] > C[c]) lfc.push_back((uint8_t)c);
        std::stable_sort(lfc.begin(), lfc.end(), [&](uint8_t x, uint8_t y) { return C[x + 1] - C[x] > C[y + 1] - C[y]; });
    }
    uint64_t LF(uint64_t i) const {
        uint64_t a, b;
        for (uint8_t c : lfc) {
            if (ef) { er.rank2(c, i, i + 1, a, b); if (b > a) return C[c] + a; }
            else if (wide ? rr64.lf_in(c, i, a) : rr.lf_in(c, i, a)) return C[c] + a;
        }
        die("rlcsa::LF: no character at row " + std::to_string(i));
        return 0;
    }

    // where c's runs start among all the runs, as first_run counts them, and an upper bound on the sum
    uint64_t run_base(unsigned char c) const { return ef ? 0 : wide ? rr64.off[c] : rr.off[c]; }
    uint64_t run_universe() const { return (ef ? R : wide ? rr64.off[256] : rr.off[256]) + 2; }
    // extend, when the first of c's runs at or after sp and at or after ep are both within runs rl..rh
    // (first_run of the two ends of an interval containing [sp, ep): the binary searches start there)
    void extend_in(uint64_t &sp, uint64_t &ep, unsigned char c, uint64_t rl, uint64_t rh) const {
        if (C[c + 1] == C[c]) { sp = ep = 0; return; }
        KATKA_STEP();
        uint64_t a, b;
        if (ef) er.rank2_in(c, sp, ep, rl, rh, a, b);
        else if (wide) { a = rr64.rank_in(c, sp, rl, rh); b = rr64.rank_in(c, ep, rl, rh); }
        else { a = rr.rank_in(c, sp, rl, rh); b = rr.rank_in(c, ep, rl, rh); }
        sp = C[c] + a; ep = C[c] + b;
    }

    // one backward step: [sp, ep) for P -> [sp, ep) for cP (empty if sp >= ep)
    void extend(uint64_t &sp, uint64_t &ep, unsigned char c) const {
        if (C[c + 1] == C[c]) { sp = ep = 0; return; }
        KATKA_STEP();
        uint64_t a, b;
        if (ef) er.rank2(c, sp, ep, a, b);
        else if (wide) { a = rr64.rank(c, sp); b = rr64.rank(c, ep); }
        else { a = rr.rank(c, sp); b = rr.rank(c, ep); }
        sp = C[c] + a; ep = C[c] + b;
    }

    uint64_t bytes() const { return 8 * 259 + (ef ? er.bytes() : wide ? rr64.bytes() : rr.bytes()); }
    void save(const std::string &file) const {
        FILE *f = open_or_die(file, "wb");
        wmagic(f, "KTCSA001"); wr(f, n); wr(f, R); fwrite(C, 8, 257, f);
        uint8_t e = ef ? 1 : wide ? 2 : 0; wr(f, e);
        if (ef) er.save(f); else if (wide) rr64.save(f); else rr.save(f);
        fclose(f);
    }
    void load(const std::string &file) {
        FILE *f = open_or_die(file, "rb");
        rmagic(f, "KTCSA001", file);
        uint8_t e = 0;
        rd(f, n); rd(f, R);
        if (fread(C, 8, 257, f) != 257) die("cannot read " + file);
        rd(f, e);
        ef = e == 1; wide = e == 2;
        if (ef) er.load(f); else if (wide) rr64.load(f); else rr.load(f);
        fclose(f);
        prepare_lf();
    }
};

}  // namespace katka
