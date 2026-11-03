// Input files for the published STPD-index (github.com/regindex/STPD-index), built in O(r) memory plus
// the text, without the suffix tree its stpd-sampling builds (about 68 bytes per character).
//
// Their index of a text T (over A/C/G/T only) is built by stpd-build from four files next to T, which it
// uses instead of running stpd-sampling when T.colex_m exists:
//   T.colex_m  the st-colex- sampling of T: 5-byte positions m (prefixes T[0..m]), in colex order of the
//              prefixes, without the smallest;
//   T.rbwt     the BWT of T_rev = reverse(T) followed by a 0 byte (one byte per row; the terminator 0);
//   T.pa       5 bytes per row of that BWT: N - 1 - SA[row], N = |T| + 1 (the prefix array);
//   T.lcs      5 bytes per row: the LCP of the suffixes of T_rev at rows row - 1 and row (0 for row 0).
// We index S = reverse(T) (so S$ is their T_rev), and:
//   - our STPD sample of S is their st-colex- sampling of T, b <-> m = |S| - 1 - b, except that theirs also
//     has m = 0, which we add (tests/regress_stpd.sh checks the two samplings against each other);
//   - our BWT file of S$ is T.rbwt (terminators rewritten as 0);
//   - SA of S$ is streamed from the last row to the first with phi (tag_build.hpp), and the LCP from
//     PLCP, computed in text order with phi (Kasai et al.'s argument: PLCP[i+1] >= PLCP[i] - 1) and kept as
//     a 2N-bit vector (bit PLCP[i] + 2i set).
// katka stpd-prep writes T (both strands, split into A/C/G/T pieces), S = reverse(T) and the piece starts.
#pragma once
#include "prep.hpp"
#include "stpd.hpp"

namespace katka {

// T = for each genome g, its A/C/G/T pieces, then those of rc(g); S = reverse(T)
inline void stpd_prep(const std::string &listfile, const std::string &out, bool split) {
    FILE *ft = open_or_die(out + ".T", "wb");
    std::ofstream pieces(out + ".pieces");
    u64 pos = 0, ng = 0, reflen = 0;
    auto put = [&](const std::string &g) {               // the A/C/G/T pieces of g
        u64 i = 0;
        while (i < g.size()) {
            while (i < g.size() && base_code(g[i]) < 0) ++i;
            u64 j = i;
            while (j < g.size() && base_code(g[j]) >= 0) ++j;
            if (j > i) { pieces << pos << "\n"; fwrite(g.data() + i, 1, j - i, ft); pos += j - i; }
            i = j;
        }
    };
    for (auto &[tag, path] : read_list(listfile))
        for_each_genome(path, split, [&](const std::string &, std::string &g) {
            put(g); put(reverse_complement(g));
            if (ng++ == 0) reflen = pos;
        });
    fclose(ft);
    // S = reverse(T), in chunks from the end
    FILE *in = open_or_die(out + ".T", "rb"), *fs = open_or_die(out + ".S", "wb");
    std::vector<char> buf(1 << 24);
    for (u64 end = pos; end > 0;) {
        u64 len = std::min<u64>(buf.size(), end), start = end - len;
        fseeko(in, start, SEEK_SET);
        if (fread(buf.data(), 1, len, in) != len) die("cannot read " + out + ".T");
        std::reverse(buf.begin(), buf.begin() + len);
        fwrite(buf.data(), 1, len, fs);
        end = start;
    }
    fclose(in); fclose(fs);
    { FILE *f = open_or_die(out + ".reflen", "w"); fprintf(f, "%lu\n", reflen); fclose(f); }
    write_info(out, false);
    fprintf(stderr, "stpd-prep: %lu genomes, |T| = %lu (both strands, A/C/G/T pieces); the first genome (both strands) has %lu characters (%s.reflen)\n",
            ng, pos, reflen, out.c_str());
}

inline void stpd_export(const std::string &prefix, const std::string &bwtfile) {
    double t0 = now();
    const std::string T = prefix + ".T";
    rlbwt B;
    { FILE *f = open_or_die(prefix + ".rlbwt", "rb"); rmagic(f, "KTRLBW01", prefix + ".rlbwt"); B.load(f); fclose(f); }
    if (B.C[2] - B.C[1] != 1) die("stpd-export needs a BWT with one terminator");
    const u64 N = B.n, n = N - 1;                       // N = |S$| = |T| + 1
    if (file_size(prefix + ".S") != n || file_size(T) != n) die(prefix + ".S, " + T + " and the BWT do not match");
    // T.rbwt: the BWT of S$, terminator as 0
    {
        FILE *in = open_or_die(bwtfile, "rb"), *o = open_or_die(T + ".rbwt", "wb");
        std::vector<unsigned char> buf(1 << 24);
        size_t got;
        while ((got = fread(buf.data(), 1, buf.size(), in)) > 0) { for (size_t i = 0; i < got; ++i) if (buf[i] <= 2) buf[i] = 0; fwrite(buf.data(), 1, got, o); }
        fclose(in); fclose(o);
    }
    // T.colex_m: our sample, b -> n - 1 - b, in row order, with m = 0 (b = n - 1) inserted at its rank
    {
        stpd_index D; D.load(prefix + ".stpd", prefix + ".S");
        const u64 bz = n - 1;
        u64 at = 0, cnt = D.pos.size();
        bool have = false;
        {   // the rank of the suffix S[n-1..]$ among the sampled suffixes
            u64 lo = 0, hi = cnt;
            while (lo < hi) {
                u64 mid = (lo + hi) / 2, b = D.pos[mid];
                if (b == bz) { have = true; lo = mid; break; }
                u64 m = std::min(N - b, N - bz), l = lce_bytes(D.S.data() + b, D.S.data() + bz, m);
                bool less = D.S[b + l] < D.S[bz + l];
                if (less) lo = mid + 1; else hi = mid;
            }
            at = lo;
        }
        FILE *o = open_or_die(T + ".colex_m", "wb");
        auto w5 = [&](u64 x) { fwrite(&x, 5, 1, o); };
        for (u64 i = 0; i < cnt; ++i) {
            if (i == at && !have) w5(0);
            w5(n - 1 - D.pos[i]);
        }
        if (at == cnt && !have) w5(0);
        fclose(o);
        fprintf(stderr, "stpd-export: %lu sampled positions (+ position 0%s)\n", cnt, have ? ", already there" : "");
    }
    // PLCP of S$ in text order with phi, as a 2N-bit vector
    std::vector<uint8_t> S(N);
    { FILE *f = open_or_die(prefix + ".S", "rb"); if (fread(S.data(), 1, n, f) != n) die("cannot read S"); fclose(f); S[n] = 1; }
    genome_table G;                                      // (one terminator: not used)
    phi_samples P = sample_rindex(prefix + ".S", bwtfile, B, G);
    const u64 sa0 = n;                                   // SA[0]: the suffix $
    BitRank plcp; plcp.init(2 * N + 2);
    {
        u64 h = 0;
        for (u64 i = 0; i < N; ++i) {
            if (i == sa0) h = 0;
            else {
                u64 j = P.phi(i);
                u64 m = N - std::max(i, j);
                if (h < m) h += lce_bytes(S.data() + i + h, S.data() + j + h, m - h);
            }
            plcp.set(h + 2 * i);
            if (h) --h;
        }
        plcp.finish();
    }
    fprintf(stderr, "stpd-export: PLCP done (%.1f s)\n", now() - t0);
    std::vector<uint8_t>().swap(S);
    // T.pa and T.lcs, streaming SA from the last row to the first
    {
        FILE *fp = open_or_die(T + ".pa", "wb"), *fl = open_or_die(T + ".lcs", "wb");
        const u64 CH = 1 << 20;
        std::vector<uint8_t> bp(5 * CH), bl(5 * CH);
        u64 d = P.esa[P.R - 1];                          // SA[N-1]
        u64 row = N, filled = 0;                         // rows [row, row + filled) are buffered
        auto flush = [&]() {
            if (!filled) return;
            fseeko(fp, 5 * row, SEEK_SET); fwrite(bp.data() + 5 * (CH - filled), 1, 5 * filled, fp);
            fseeko(fl, 5 * row, SEEK_SET); fwrite(bl.data() + 5 * (CH - filled), 1, 5 * filled, fl);
            filled = 0;
        };
        while (row > 0) {
            --row;
            u64 slot = CH - 1 - filled;
            u64 pa = N - 1 - d, lcs = row == 0 ? 0 : plcp.select1(d) - 2 * d;
            memcpy(&bp[5 * slot], &pa, 5); memcpy(&bl[5 * slot], &lcs, 5);
            ++filled;
            if (filled == CH) flush();
            if (row) d = P.phi(d);
        }
        flush();
        fclose(fp); fclose(fl);
    }
    fprintf(stderr, "stpd-export: wrote %s.{colex_m,rbwt,pa,lcs} for %lu rows (%.1f s)\n", T.c_str(), N, now() - t0);
}

}  // namespace katka
