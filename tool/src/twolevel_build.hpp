// Building the two-level index (twolevel.hpp) from the text S and its RLBWT (one terminator), in a
// streaming pass over S: only the distinct phrases and the parse are kept in memory.  Ported from More
// RLBWT Tricks' build_twolevel.cpp, for one text of both strands instead of T and T'.
#pragma once
#include "libsais.h"
#include "twolevel.hpp"
#include <unordered_map>

namespace katka {

struct twolevel_options {
    int k = 51, s = 12;
    bool blocks = true;
    bool frequent = false;          // continuation: the most frequent successor instead of the smallest
};

inline void build_twolevel(const std::string &Sfile, const rlbwt &bwt, const stepper &st, const twolevel_options &o, two_level &L) {
    const int K = o.k;
    if (o.s < 1 || o.s >= o.k || o.s > 32) die("two-level: need 1 <= s < k and s <= 32");
    if (bwt.C[2] - bwt.C[1] != 1) die("two-level indexing needs a BWT with one terminator (BWT(S$))");
    double t0 = now();
    L.K = o.k; L.SM = o.s;
    std::vector<uint32_t> pz;                // the parse, as indices of distinct phrases (later, ids)
    std::vector<uint8_t> dict;               // the distinct phrases, concatenated
    std::vector<uint64_t> doff, rhash;
    std::unordered_map<uint64_t, uint32_t> seen;
    bool zero_is_syncmer = false;
    uint64_t n = 0;
    {
        FILE *f = open_or_die(Sfile, "rb");
        std::vector<uint8_t> buf, raw(1 << 24);
        uint64_t base = 0, open = 0, scanned = 0;
        auto add_phrase = [&](uint64_t s, uint64_t e) {
            const uint8_t *p = buf.data() + (s - base);
            uint64_t len = e - s, hv = phrase_hash(p, len);
            auto it = seen.find(hv);
            if (it == seen.end()) {
                uint32_t idx = (uint32_t)rhash.size();
                seen.emplace(hv, idx); rhash.push_back(hv); doff.push_back(dict.size());
                dict.insert(dict.end(), p, p + len);
                pz.push_back(idx);
            } else {
                uint32_t idx = it->second;
                uint64_t a = doff[idx], b2 = idx + 1 < doff.size() ? doff[idx + 1] : dict.size();
                if (b2 - a != len || memcmp(dict.data() + a, p, len)) throw std::runtime_error("phrase hash collision");
                pz.push_back(idx);
            }
        };
        bool more = true;
        while (more) {
            size_t got = fread(raw.data(), 1, raw.size(), f);
            buf.insert(buf.end(), raw.begin(), raw.begin() + got);
            if (got == 0) { buf.push_back(1); more = false; }      // the terminator $
            const uint64_t end = base + buf.size();
            if (end >= (uint64_t)K && end - K + 1 > scanned) {
                const uint64_t from = scanned;
                for_each_syncmer(buf.data() + (from - base), end - from, K, o.s, [&](uint64_t q) {
                    uint64_t a = from + q;
                    if (a == 0) { zero_is_syncmer = true; return; }
                    add_phrase(open, a + K);
                    open = a;
                });
                scanned = end - K + 1;
            }
            if (!more) { n = end; add_phrase(open, n); break; }
            uint64_t keep = std::min(open, scanned);
            if (keep > base) { buf.erase(buf.begin(), buf.begin() + (keep - base)); base = keep; }
        }
        fclose(f);
    }
    if (n != bwt.n) die("two-level: |S$| = " + std::to_string(n) + " but the BWT has " + std::to_string(bwt.n) + " rows");
    doff.push_back(dict.size());
    seen = std::unordered_map<uint64_t, uint32_t>();
    const uint64_t m = pz.size();
    if (m + 1 >= (1ull << 31)) die("too many phrases for a 32-bit parse suffix array");
    const uint32_t sigma = (uint32_t)rhash.size();
    auto ptr = [&](uint32_t a) { return dict.data() + doff[a]; };
    auto dlen = [&](uint32_t a) { return doff[a + 1] - doff[a]; };
    fprintf(stderr, "two-level: %lu characters; %lu phrases, %u distinct (%.1f MB)\n", n, m, sigma, dict.size() / 1e6);

    // lexicographic ids 1..sigma (byte order; the dictionary is prefix-free)
    std::vector<uint32_t> order(sigma);
    for (uint32_t a = 0; a < sigma; ++a) order[a] = a;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t c) {
        uint64_t la = dlen(a), lc = dlen(c);
        int r = memcmp(ptr(a), ptr(c), std::min(la, lc));
        if (r == 0 && la != lc) throw std::runtime_error("dictionary is not prefix-free");
        return r < 0;
    });
    std::vector<uint32_t> idof(sigma);
    for (uint32_t r = 0; r < sigma; ++r) idof[order[r]] = r + 1;
    for (auto &x : pz) x = idof[x];
    L.n = n; L.m = m; L.sigma = sigma;

    // the dictionary: phrases over A/C/G/T, with the id of their reverse complements
    {
        std::unordered_map<uint64_t, uint32_t> byhash;
        byhash.reserve(sigma * 2);
        for (uint32_t a = 0; a < sigma; ++a) byhash.emplace(rhash[a], idof[a]);
        std::vector<std::pair<uint64_t, std::pair<uint32_t, uint32_t>>> v;
        std::vector<uint8_t> tmp;
        uint64_t missing = 0;
        for (uint32_t a = 0; a < sigma; ++a) {
            const uint8_t *p = ptr(a);
            uint64_t len = dlen(a);
            bool acgt = true;
            for (uint64_t i = 0; i < len && acgt; ++i) acgt = base_code(p[i]) >= 0;
            if (!acgt) continue;
            tmp.resize(len);
            for (uint64_t i = 0; i < len; ++i) tmp[i] = complement(p[len - 1 - i]);
            auto it = byhash.find(phrase_hash(tmp.data(), len));
            uint32_t rid = 0;
            if (it != byhash.end() && dlen(order[it->second - 1]) == len && !memcmp(ptr(order[it->second - 1]), tmp.data(), len)) rid = it->second;
            else ++missing;
            v.push_back({rhash[a], {idof[a], rid}});
        }
        L.dict.build(v);
        fprintf(stderr, "dictionary: %zu phrases over A/C/G/T, %lu without a reverse complement (the first phrase of S can be one)\n", v.size(), missing);
        if (missing > 2) die("two-level: phrases without reverse complements; is S closed under reverse complement (katka prep)?");
    }

    // distinct syncmers: phrases with equal first K characters are adjacent in id order
    std::vector<uint8_t> syn;
    {
        uint32_t skip = zero_is_syncmer ? 0 : pz[0];
        const uint8_t *prev = nullptr;
        for (uint32_t r = 0; r < sigma; ++r) {
            if (r + 1 == skip) continue;
            const uint8_t *p = ptr(order[r]);
            if (dlen(order[r]) < (uint64_t)K) continue;
            if (!prev || memcmp(prev, p, K)) { syn.insert(syn.end(), p, p + K); prev = p; }
        }
    }
    const size_t nsyn = syn.size() / K;
    std::vector<uint32_t>().swap(order); std::vector<uint32_t>().swap(idof);
    std::vector<uint64_t>().swap(rhash); std::vector<uint8_t>().swap(dict); std::vector<uint64_t>().swap(doff);

    // B: the union of the syncmers' intervals, and the row of position 0 (the row whose BWT character is $)
    {
        std::vector<std::pair<uint64_t, uint64_t>> iv;
        iv.reserve(nsyn + 1);
        for (size_t q = 0; q < nsyn; ++q) {
            const uint8_t *s = syn.data() + q * K;
            interval I = st.full();
            for (int c = K - 1; c >= 0; --c) if (!st.step(I, s[c])) die("two-level: a syncmer does not occur in the BWT");
            iv.push_back({I.sp, I.ep - I.sp});
        }
        if (!zero_is_syncmer) {
            u64 row = UINT64_MAX;
            for (u64 q = 0; q < bwt.R; ++q) if (bwt.head(q) == 1) { row = bwt.run_start(q); break; }
            iv.push_back({row, 1});
        }
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<uint64_t, uint64_t>> runs;
        uint64_t total = 0;
        for (auto &x : iv) {
            if (!runs.empty() && runs.back().first + runs.back().second > x.first) die("two-level: overlapping syncmer intervals");
            if (!runs.empty() && runs.back().first + runs.back().second == x.first) runs.back().second += x.second;
            else runs.push_back(x);
            total += x.second;
        }
        if (total != m) die("two-level: B has " + std::to_string(total) + " ones, not m = " + std::to_string(m) + "; is the BWT that of S?");
        L.B.build(runs, n);
        fprintf(stderr, "B: %zu runs, %.2f MB\n", runs.size(), L.B.bytes() / 1e6);
    }
    std::vector<uint8_t>().swap(syn);

    // the parse BWT
    std::vector<int32_t> SA(m + 1);
    {
        std::vector<int32_t> Tz(m + 1);
        for (uint64_t i = 0; i < m; ++i) Tz[i] = (int32_t)pz[i];
        Tz[m] = 0;
        if (libsais_int(Tz.data(), SA.data(), (int32_t)(m + 1), (int32_t)sigma + 1, 0) != 0) die("libsais_int failed");
        std::vector<uint32_t> bw(m + 1);
        for (uint64_t r = 0; r <= m; ++r) bw[r] = SA[r] == 0 ? 0 : pz[SA[r] - 1];
        L.pbwt.build(bw, sigma + 1);
    }
    fprintf(stderr, "parse BWT: %lu runs, %.2f MB\n", L.pbwt.nruns, L.pbwt.bytes() / 1e6);
    L.has_blocks = o.blocks;
    if (!o.blocks) { fprintf(stderr, "two-level index built in %.1f s\n", now() - t0); return; }

    // successors and blocks
    std::vector<uint32_t> succ(sigma + 1, 0);
    if (!o.frequent) {
        for (uint64_t i = 0; i + 1 < m; ++i) { uint32_t &s = succ[pz[i]]; if (!s || pz[i + 1] < s) s = pz[i + 1]; }
    } else {
        std::vector<uint64_t> pr;
        pr.reserve(m ? m - 1 : 0);
        for (uint64_t i = 0; i + 1 < m; ++i) pr.push_back((uint64_t)pz[i] << 32 | pz[i + 1]);
        std::sort(pr.begin(), pr.end());
        std::vector<uint64_t> best(sigma + 1, 0);
        for (uint64_t a = 0; a < pr.size();) {
            uint64_t b = a;
            while (b < pr.size() && pr[b] == pr[a]) ++b;
            uint32_t x = (uint32_t)(pr[a] >> 32);
            if (b - a > best[x]) { best[x] = b - a; succ[x] = (uint32_t)pr[a]; }
            a = b;
        }
    }
    std::vector<uint8_t> cont(m + 1, 0);
    L.ccont.assign(sigma / 64 + 1, 0);
    for (uint64_t i = 1; i < m; ++i) {
        cont[i] = succ[pz[i - 1]] == pz[i];
        if (cont[i]) L.ccont[pz[i] >> 6] |= 1ull << (pz[i] & 63);
    }
    // symbols: keys d << 33 | sym_key(j, side), sorted in the symbol order (the last block gets +, which is
    // safe since its last phrase occurs only once)
    std::vector<uint64_t> pairs;
    {
        uint64_t start = 0;
        for (uint64_t i = 1; i <= m; ++i)
            if (i == m || !cont[i]) {
                uint64_t j = i - 1 - start;
                if (j >= (1ull << 32) - 1) die("block too long");
                bool plus = !o.frequent || (i < m && pz[i] > succ[pz[i - 1]]);
                pairs.push_back(((uint64_t)pz[start] << 33) | two_level::sym_key(j, plus));
                start = i;
            }
    }
    L.nblocks = pairs.size();
    std::vector<uint64_t> syms(pairs);
    std::sort(syms.begin(), syms.end());
    syms.erase(std::unique(syms.begin(), syms.end()), syms.end());
    L.nsym = syms.size();
    {
        std::vector<uint64_t> ds(sigma + 2, 0), jvv(syms.size());
        std::vector<uint8_t> plus(syms.size());
        for (size_t q = 0; q < syms.size(); ++q) {
            ds[(syms[q] >> 33) + 1]++;
            plus[q] = (syms[q] >> 32) & 1;
            uint64_t low = syms[q] & 0xffffffffull;
            jvv[q] = plus[q] ? 0xffffffffull - low : low;
        }
        for (uint32_t d = 0; d <= sigma; ++d) ds[d + 1] += ds[d];
        L.dstart.build(ds, L.nsym + 1);
        L.jv.build(jvv);
        L.sided = o.frequent;
        L.sidev.assign(o.frequent ? jvv.size() / 64 + 1 : 0, 0);
        if (o.frequent) for (size_t q = 0; q < jvv.size(); ++q) if (plus[q]) L.sidev[q >> 6] |= 1ull << (q & 63);
    }
    L.succp.build_from(succ);
    std::vector<uint32_t>().swap(succ);
    {
        std::vector<uint32_t> prevsym(m + 1, 0);    // at block starts p > 0 (and m): the symbol of the block ending at p-1
        uint64_t q = 0;
        for (uint64_t i = 1; i <= m; ++i)
            if (i == m || !cont[i]) prevsym[i] = (uint32_t)(std::lower_bound(syms.begin(), syms.end(), pairs[q++]) - syms.begin()) + 1;
        std::vector<uint64_t> rows;
        std::vector<uint32_t> bw;
        rows.reserve(L.nblocks + 1); bw.reserve(L.nblocks + 1);
        for (uint64_t r = 0; r <= m; ++r) {
            uint64_t p = (uint64_t)SA[r];
            if (p == m || !cont[p]) { rows.push_back(r); bw.push_back(p == 0 ? 0 : prevsym[p]); }
        }
        if (rows.size() != L.nblocks + 1) die("internal error: block row count");
        L.BS.build(rows, m + 2);
        L.bbwt.build(bw, (uint32_t)L.nsym + 1);
    }
    fprintf(stderr, "blocks: %lu (%.2f phrases each), %lu symbols; block BWT %lu runs; built in %.1f s\n", L.nblocks,
            (double)m / L.nblocks, L.nsym, L.bbwt.nruns, now() - t0);
}

}  // namespace katka
