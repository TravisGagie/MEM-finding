// The run-length compressed tag array (Section 6).  The tag array T[0..n) gives, for each BWT row, the tag
// (genus or species) of the suffix starting there.  We store it by its tag runs: rho maximal runs of rows
// with the same tag.
//   RB : Elias-Fano bitvector marking the first row of every tag run
//   L  : the tag of every run (14-bit integers), or of the sampled runs only (s > 1, below)
//
// Counting (Section 6): the number of rows of [sp, ep) with each tag, by scanning the runs that overlap
// [sp, ep) and adding up their clipped lengths.  This is what proportional credit needs.
//
// Listing (equal credit, and the smallest and largest tags for the range scheme): the distinct tags in
// [sp, ep), by Muthukrishnan's method without the C array (Sadakane's trick, as in Olbrich and Ohlebusch):
// an RMQ over C[i] = 1 + the previous run with the same tag as run i (0 if none), plus marks for the tags
// already reported.  The RMQ is stored in its own file (out.tag.rmq), since counting does not need it.
//
// Sampling (s > 1; the paper's Table 3 uses s = 1): following LF from a run's first row stays in the same
// document except at a document's first character, so a run's tag is the tag of the run containing LF of its
// first row.  Only some runs' tags are stored (SB marks them), so that every run reaches a sampled run in
// fewer than s LF steps; tag(i) walks there.
//
// File (katka-tagbuild): "RZTAG002", n, rho, s, RB with its rank and select, [SB with its rank, if s > 1], L.
// Files of the older format "RZTAG001" carry the RMQ at the end instead of in out.tag.rmq.
#pragma once
#include "rlbwt.hpp"
#include <sdsl/rmq_support.hpp>
#include <algorithm>

namespace katka {

struct tag_array {
    u64 n = 0, rho = 0, s = 1;                      // rows, tag runs, sampling parameter
    sdsl::sd_vector<> RB; sdsl::rank_support_sd<> RBr; sdsl::select_support_sd<> RBs;
    sdsl::bit_vector SB; sdsl::rank_support_v5<> SBr;
    sdsl::int_vector<> L;
    sdsl::rmq_succinct_sct<> rmq; bool has_rmq = false;
    const rlbwt *bwt = nullptr;                     // for LF, when s > 1

    // scratch space, one copy per thread (katka-classify -j)
    static inline thread_local std::vector<uint32_t> cnt;      // counting: rows per tag
    static inline thread_local std::vector<char> mark;         // listing: tags already reported
    static inline thread_local std::vector<u64> touched;       // the tags whose cnt or mark was set
    static inline thread_local u64 lf_steps = 0;

    u64 run_of(u64 row) const { return RBr(row + 1) - 1; }
    u64 run_start(u64 i) const { return RBs(i + 1); }
    u64 run_end(u64 i) const { return i + 1 < rho ? (u64)RBs(i + 2) : n; }

    // the tag of run i
    u64 tag(u64 i) {
        if (s == 1) return L[i];
        while (!SB[i]) {                            // the run containing LF of the run's first row
            u64 b = run_start(i); unsigned char c = bwt->at(b);
            i = run_of(bwt->C[c] + bwt->rank(b, c)); ++lf_steps;
        }
        return L[SBr(i)];
    }

    // (tag, number of rows of [sp, ep) with that tag), sorted by tag
    void count(u64 sp, u64 ep, std::vector<std::pair<u64, u64>> &out) {
        out.clear();
        u64 a = run_of(sp), b = run_of(ep - 1);
        if (a == b) { out.push_back({tag(a), (uint32_t)(ep - sp)}); return; }   // inside one tag run (most intervals)
        if (cnt.empty()) cnt.assign(1 << 16, 0);
        u64 st = sp;                                    // each run's end is the next run's start: one select per run
        for (u64 i = a; i <= b; ++i) {
            u64 en = i < b ? run_start(i + 1) : ep, t = tag(i);
            if (t >= cnt.size()) cnt.resize(t + 1, 0);
            if (cnt[t] == 0) touched.push_back(t);
            cnt[t] += (uint32_t)(en - st);
            st = en;
        }
        for (u64 t : touched) { out.push_back({t, cnt[t]}); cnt[t] = 0; }
        touched.clear();
        std::sort(out.begin(), out.end());
    }

    // the distinct tags of rows [sp, ep), appended to out (unsorted)
    void list(u64 sp, u64 ep, std::vector<u64> &out) {
        if (mark.empty()) mark.assign(1 << 16, 0);
        u64 a = run_of(sp), b = run_of(ep - 1);
        report(tag(a), out);
        if (b > a) list_rec(a + 1, b, out);
        for (u64 t : touched) mark[t] = 0;
        touched.clear();
    }

    u64 bytes() const {
        return sdsl::size_in_bytes(RB) + sdsl::size_in_bytes(RBr) + sdsl::size_in_bytes(RBs) + sdsl::size_in_bytes(L)
               + (s > 1 ? sdsl::size_in_bytes(SB) + sdsl::size_in_bytes(SBr) : 0) + (has_rmq ? sdsl::size_in_bytes(rmq) : 0);
    }
    void save(const std::string &f) const {
        std::ofstream o(f, std::ios::binary);
        o.write("RZTAG002", 8); o.write((char *)&n, 8); o.write((char *)&rho, 8); o.write((char *)&s, 8);
        RB.serialize(o); RBr.serialize(o); RBs.serialize(o);
        if (s > 1) { SB.serialize(o); SBr.serialize(o); }
        L.serialize(o);
        if (has_rmq) { std::ofstream r(f + ".rmq", std::ios::binary); r.write("RZRMQ001", 8); rmq.serialize(r); }
    }
    // with_rmq: also load the RMQ, which listing needs
    void load(const std::string &f, bool with_rmq) {
        std::ifstream i(f, std::ios::binary); char mg[8];
        if (!i.read(mg, 8)) die("cannot read " + f);
        std::string magic(mg, 8);
        if (magic == "RZGTAG01") die(f + " is a grammar-compressed tag array (Section 6), which this version does not read");
        if (magic != "RZTAG001" && magic != "RZTAG002") die(f + " is not a tag array");
        i.read((char *)&n, 8); i.read((char *)&rho, 8); i.read((char *)&s, 8);
        RB.load(i); RBr.load(i, &RB); RBs.load(i, &RB);
        if (s > 1) { SB.load(i); SBr.load(i, &SB); }
        L.load(i);
        if (!i) die("cannot read " + f);
        if (!with_rmq) return;
        if (magic == "RZTAG001") {
            rmq.load(i);
            if (!i) die("cannot read the RMQ of " + f);
        } else {
            std::ifstream r(f + ".rmq", std::ios::binary);
            if (!r.read(mg, 8) || std::string(mg, 8) != "RZRMQ001") die(f + ".rmq missing: listing needs the RMQ (katka-tagbuild without -n)");
            rmq.load(r);
            if (!r) die("cannot read " + f + ".rmq");
        }
        has_rmq = true;
    }

private:
    void report(u64 t, std::vector<u64> &out) {
        if (t >= mark.size()) mark.resize(t + 1, 0);
        mark[t] = 1; touched.push_back(t); out.push_back(t);
    }
    // runs [l, r]: the run k with the smallest C[k] has a tag that does not occur in runs [a, k) iff its
    // tag is not marked; if it is marked, every tag in [l, r] already occurs earlier, so stop
    void list_rec(u64 l, u64 r, std::vector<u64> &out) {
        u64 k = rmq(l, r), t = tag(k);
        if (t < mark.size() && mark[t]) return;
        report(t, out);
        if (k > l) list_rec(l, k - 1, out);
        if (k < r) list_rec(k + 1, r, out);
    }
};

}  // namespace katka
