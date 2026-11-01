// The RLBWT (Section 2 of the paper, "RLBWTs and RLCSAs"): the Makinen-Navarro run-length FM-index,
//   starts  Elias-Fano list of the first row of every run (r of them);
//   heads   a wavelet matrix over the run heads (the character of each run);
//   cum[c]  Elias-Fano list, for each character c, of the positions (among the c's) of the last c of each
//           run of c.
// rank(i, c) = number of c's in BWT[0, i) needs one predecessor search on starts, a rank on heads and a
// select on cum[c].  It gives LF, for sampled tags and the r-index samples, and a backward step.
//
// The BWT file has one byte per row; terminators may be bytes 0, 1 or 2 (one per string, from libsais,
// Big-BWT or pfp-merge) and are all read as byte 1.
#pragma once
#include "bits.hpp"
#include <memory>

namespace katka {

// Calls f(run_index, start, length, character) for every run of a BWT file, in order.
template <class F> void for_each_bwt_run(const std::string &file, F f) {
    FILE *fp = open_or_die(file, "rb");
    std::vector<unsigned char> buf(1 << 24);
    u64 pos = 0, k = 0, start = 0, len = 0;
    int last = -1;
    size_t got;
    while ((got = fread(buf.data(), 1, buf.size(), fp)) > 0) {
        for (size_t i = 0; i < got; ++i, ++pos) {
            int c = buf[i] > 2 ? buf[i] : 1;
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

struct rlbwt {
    u64 n = 0, R = 0;                   // rows, runs
    u64 C[257];                         // C[c] = number of characters smaller than c
    uint8_t code[256];                  // byte -> symbol of the wavelet matrix (255: absent)
    uint8_t byte_of[256];               // symbol -> byte
    EliasFano starts;
    WaveletMatrix heads;
    std::vector<EliasFano> cum = std::vector<EliasFano>(256);

    rlbwt() { std::fill(C, C + 257, 0); std::fill(code, code + 256, 255); }

    void build(const std::string &bwtfile) {
        u64 count[256] = {0}, nruns[256] = {0};
        for_each_bwt_run(bwtfile, [&](u64 k, u64 s, u64 len, unsigned char c) { count[c] += len; nruns[c]++; R = k + 1; n = s + len; });
        C[0] = 0;
        for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + count[c];
        uint32_t sig = 0;
        for (int c = 0; c < 256; ++c) code[c] = count[c] ? (uint8_t)sig++ : 255;
        set_byte_of();
        std::vector<u64> st; st.reserve(R);
        std::vector<uint8_t> hv; hv.reserve(R);
        std::vector<std::vector<u64>> ends(256);
        for (int c = 0; c < 256; ++c) ends[c].reserve(nruns[c]);
        u64 seen[256] = {0};
        for_each_bwt_run(bwtfile, [&](u64, u64 s, u64 len, unsigned char c) {
            st.push_back(s); hv.push_back(code[c]); seen[c] += len; ends[c].push_back(seen[c] - 1);
        });
        starts.build(st, n);
        std::vector<u64>().swap(st);
        for (int c = 0; c < 256; ++c) if (count[c]) { cum[c].build(ends[c], count[c]); std::vector<u64>().swap(ends[c]); }
        heads.build(hv, std::max<uint32_t>(sig, 2));
    }

    u64 run_of(u64 i) const { return starts.count_less(i + 1) - 1; }       // the run containing row i
    u64 run_start(u64 q) const { return starts.access(q); }
    unsigned char at(u64 i) const { return head(run_of(i)); }               // BWT[i]
    unsigned char head(u64 q) const { return byte_of[heads.access(q)]; }
    void set_byte_of() { std::fill(byte_of, byte_of + 256, 0); for (int c = 0; c < 256; ++c) if (code[c] != 255) byte_of[code[c]] = (uint8_t)c; }

    // number of c's in BWT[0, i)
    u64 rank(u64 i, unsigned char c) const {
        if (i == 0 || C[c + 1] == C[c]) return 0;
        if (i >= n) return C[c + 1] - C[c];
        u64 q = run_of(i - 1), s = starts.access(q);
        uint32_t cc = code[c];
        auto is = heads.inverse_select(q);                   // (runs of heads[q] before run q, heads[q])
        uint32_t hq = is.second;
        u64 t = hq == cc ? is.first : heads.rank(cc, q);     // runs of c before run q
        u64 r = t == 0 ? 0 : cum[c].access(t - 1) + 1;       // c's in those runs
        if (hq == cc) r += i - s;
        return r;
    }
    u64 LF(u64 i) const { unsigned char c = at(i); return C[c] + rank(i, c); }

    // one backward step: [sp, ep) for P -> [sp, ep) for cP (empty if sp >= ep)
    void extend(u64 &sp, u64 &ep, unsigned char c) const {
        if (C[c + 1] == C[c]) { sp = ep = 0; return; }
        sp = C[c] + rank(sp, c);
        ep = C[c] + rank(ep, c);
    }

    u64 bytes() const {
        u64 b = 16 + 257 * 8 + 256 + starts.bytes() + heads.bytes();
        for (int c = 0; c < 256; ++c) if (C[c + 1] > C[c]) b += cum[c].bytes();
        return b;
    }
    void save(FILE *f) const {
        wr(f, n); wr(f, R); fwrite(C, 8, 257, f); fwrite(code, 1, 256, f);
        starts.save(f); heads.save(f);
        for (int c = 0; c < 256; ++c) if (C[c + 1] > C[c]) cum[c].save(f);
    }
    void load(FILE *f) {
        rd(f, n); rd(f, R);
        if (fread(C, 8, 257, f) != 257 || fread(code, 1, 256, f) != 256) throw std::runtime_error("short read");
        set_byte_of();
        starts.load(f); heads.load(f);
        for (int c = 0; c < 256; ++c) if (C[c + 1] > C[c]) cum[c].load(f);
    }
};

}  // namespace katka
