// The run-length compressed BWT (RLBWT) used for backward search (Section 5) and for LF in sampled tag
// arrays (Section 6), and the .rlbwt file that holds it.
//
// The text (Appendix A) is
//     S = G_1 X rc(G_1) X G_2 X rc(G_2) X ... G_k X rc(G_k) X
// so it is closed under reverse complement: a forward step on a pattern P is a backward step on rc(P)
// (Section 5).  The BWT file has one byte per row; the terminator(s) may be bytes 0, 1 or 2 (one per
// string, from Big-BWT or pfp-merge), and are all read as byte 1.
//
// The RLBWT is the Makinen-Navarro run-length FM-index:
//   starts : Elias-Fano bitvector over [0, n) marking the first row of every run   (r ones)
//   heads  : Huffman-shaped wavelet tree over the r run heads (the character of each run)
//   cum[c] : Elias-Fano bitvector over [0, number of c's), marking the last c of each run of c
// rank(i, c) = number of c's in BWT[0, i), from these three.
#pragma once
#include "common.hpp"
#include <sdsl/construct.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>

namespace katka {

// Elias-Fano bitvector with rank and select.  The serialized form (sd_vector, rank, select) is the
// original code's, so the files are interchangeable.
struct sparse_bv {
    sdsl::sd_vector<> bv;
    sdsl::rank_support_sd<1> rk;
    sdsl::select_support_sd<1> sl;

    sparse_bv() {}
    sparse_bv(const sparse_bv &) = delete;
    sparse_bv &operator=(const sparse_bv &) = delete;

    template <class V> void build(const V &ones, u64 len) {     // ones: sorted, distinct
        sdsl::sd_vector_builder b(len, ones.size());
        for (u64 x : ones) b.set(x);
        build_from(b);
    }
    void build_from(sdsl::sd_vector_builder &b) {
        bv = sdsl::sd_vector<>(b);
        sdsl::util::init_support(rk, &bv);
        sdsl::util::init_support(sl, &bv);
    }
    u64 rank(u64 i) const { return rk(i); }        // number of ones in [0, i)
    u64 select(u64 k) const { return sl(k); }      // position of the k-th one, k >= 1
    u64 ones() const { return rk(bv.size()); }
    u64 serialize(std::ostream &out) const { return bv.serialize(out) + rk.serialize(out) + sl.serialize(out); }
    void load(std::istream &in) { bv.load(in); rk.load(in, &bv); sl.load(in, &bv); }
};

struct rlbwt {
    u64 n = 0, R = 0;                  // BWT length, number of runs
    u64 C[257];                        // C[c] = number of characters smaller than c
    sparse_bv starts;
    sdsl::wt_huff<> heads;
    std::vector<sparse_bv> cum = std::vector<sparse_bv>(256);

    rlbwt() { std::fill(C, C + 257, 0); }
    rlbwt(const rlbwt &) = delete;
    rlbwt &operator=(const rlbwt &) = delete;

    // Calls f(run_index, start, length, character) for every run of the BWT file, in order.
    template <class F> static void for_each_run(const std::string &file, F f) {
        FILE *fp = fopen(file.c_str(), "rb");
        if (!fp) die("cannot open " + file);
        std::vector<unsigned char> buf(1 << 24);
        u64 pos = 0, k = 0, start = 0, len = 0;
        int last = -1;
        size_t got;
        while ((got = fread(buf.data(), 1, buf.size(), fp)) > 0) {
            for (size_t i = 0; i < got; ++i, ++pos) {
                int c = buf[i] > 2 ? buf[i] : 1;        // terminators 0, 1, 2 -> 1
                if (c != last) {
                    if (last >= 0) f(k++, start, len, (unsigned char)last);
                    last = c; start = pos; len = 0;
                }
                ++len;
            }
        }
        fclose(fp);
        if (last >= 0) f(k++, start, len, (unsigned char)last);
    }

    // two passes over the BWT file: count, then fill
    void build(const std::string &file) {
        u64 count[256] = {0}, nruns[256] = {0};
        for_each_run(file, [&](u64 k, u64 s, u64 len, unsigned char c) { count[c] += len; nruns[c]++; R = k + 1; n = s + len; });
        C[0] = 0;
        for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + count[c];

        sdsl::sd_vector_builder sb(n, R);
        std::vector<std::unique_ptr<sdsl::sd_vector_builder>> cb(256);
        for (int c = 0; c < 256; ++c) cb[c].reset(new sdsl::sd_vector_builder(count[c], nruns[c]));
        sdsl::int_vector<8> hv(R);
        u64 seen[256] = {0};
        for_each_run(file, [&](u64 k, u64 s, u64 len, unsigned char c) {
            sb.set(s);
            seen[c] += len;
            cb[c]->set(seen[c] - 1);
            hv[k] = c;
        });
        starts.build_from(sb);
        for (int c = 0; c < 256; ++c) cum[c].build_from(*cb[c]);
        sdsl::construct_im(heads, hv);
    }

    u64 run_of(u64 i) const { return starts.rank(i + 1) - 1; }     // the run containing row i
    unsigned char at(u64 i) const { return heads[run_of(i)]; }     // BWT[i]

    // number of c's in BWT[0, i)
    u64 rank(u64 i, unsigned char c) const {
        if (i == 0 || C[c + 1] == C[c]) return 0;
        if (i >= n) return C[c + 1] - C[c];
        u64 q = run_of(i - 1), s = starts.select(q + 1);
        auto is = heads.inverse_select(q);          // (number of heads[q] in heads[0, q), heads[q])
        u64 t = is.second == c ? is.first : heads.rank(q, c);   // runs of c before run q
        u64 r = t == 0 ? 0 : cum[c].select(t) + 1;              // c's in those runs
        if (is.second == c) r += i - s;                         // plus those of run q before row i
        return r;
    }
    // LF(i) = C[BWT[i]] + rank(i, BWT[i])
    u64 LF(u64 i) const { unsigned char c = at(i); return C[c] + rank(i, c); }

    // one backward step: [sp, ep) for P -> [sp, ep) for cP (empty if sp >= ep)
    void extend(u64 &sp, u64 &ep, unsigned char c) const {
        sp = C[c] + rank(sp, c);
        ep = C[c] + rank(ep, c);
    }

    u64 serialize(std::ostream &out) const {
        out.write((char *)&n, 8); out.write((char *)&R, 8); out.write((char *)C, 257 * 8);
        u64 w = 16 + 257 * 8 + starts.serialize(out) + heads.serialize(out);
        for (int c = 0; c < 256; ++c) w += cum[c].serialize(out);
        return w;
    }
    void load(std::istream &in) {
        in.read((char *)&n, 8); in.read((char *)&R, 8); in.read((char *)C, 257 * 8);
        starts.load(in); heads.load(in);
        for (int c = 0; c < 256; ++c) cum[c].load(in);
    }
};

// The .rlbwt file (katka-build -G): |S|, the number of strands (two per genome), the RLBWT and the strand starts
// in S.  The original code's .rz files have the same format (plus, if built with the LZ77
// grids of the rz-index of Appendix C, more after this, which we ignore), so they can be read too.
struct rlbwt_file {
    u64 N = 0, nstr = 0;
    rlbwt bwt;
    sparse_bv strand_starts;

    u64 serialize(std::ostream &out) const {
        out.write((char *)&N, 8); out.write((char *)&nstr, 8);
        return 16 + bwt.serialize(out) + strand_starts.serialize(out);
    }
    void load(const std::string &file) {
        std::ifstream in(file, std::ios::binary);
        if (!in) die("cannot open " + file);
        in.read((char *)&N, 8); in.read((char *)&nstr, 8);
        bwt.load(in); strand_starts.load(in);
        if (!in) die("cannot read " + file);
    }
};

}  // namespace katka
