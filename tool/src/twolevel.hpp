// Two-level indexing with blocks (Section 5 of the paper, "Two-level indexing with a syncmer parse" and
// "Blocks"), on the text of both strands.
//
// The text S (with its terminator $) is parsed at the starts of closed syncmers and at position 0.  A
// k-mer over A/C/G/T is a closed syncmer if its smallest s-mer under the order g is its first or its last
// (ties included).  Unlike More RLBWT Tricks, which indexed T and the reversed collection T', g is
// canonical under reverse complement: g(gamma) is a hash of the smaller of the 2-bit codes of gamma and of
// rc(gamma).  The s-mers of rc(kappa) are the reverse complements of those of kappa, in reverse order, so
// kappa is a closed syncmer if and only if rc(kappa) is.  Consequences (see NOTES.md):
//   - whether a k-mer is a closed syncmer still does not depend on its context, so the parse is still
//     prefix-free and the order lemma (Lemma 6) holds unchanged: its proof uses nothing about g;
//   - the syncmer starts of rc(P) are m - k - q for the starts q of P, in reverse order, and the phrases
//     of rc(P) are the reverse complements of P's, so one parse of the pattern serves the backward search
//     on P and the forward search (backward on rc(P)), and one dictionary, storing for each phrase of S
//     its id and the id of its reverse complement, serves both strands.
//
// Phrase i is S[b_i .. b_{i+1} + k), the last one S[b_last .. n) (it ends with $); phrases are numbered
// 1..sigma in byte order ($ = 1 < A < C < G < T < X), and the parse is their sequence followed by 0.
//   B     marks the rows of the BWT of S whose suffixes start at boundaries; the u-th 1 (from 0) is parse
//         row u + 1.  It is the union of the intervals of the distinct syncmers (plus the row of position
//         0), stored as runs.
//   pbwt  the BWT of the parse, stored like an RLCSA (EFRank).
// Blocks: succ(x) is the smallest phrase that follows x in the parse (or, with --continuation frequent,
// the most frequent one).  A block is a phrase d that is not a continuation followed by its j >= 0
// continuations; its symbol is (d, j, side).  Symbols are ordered by d, then - before +, -'s by j
// ascending and +'s by j descending, so suffixes starting at blocks are in the same order in the block
// sequence as in the parse (Lemma "blockorder" and "Other continuations").
//   BS    marks the parse rows that start blocks (and row 0);  bbwt  the BWT of the block sequence.
#pragma once
#include "search.hpp"

namespace katka {

inline uint64_t mix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27; x *= 0x94d049bb133111ebull;
    x ^= x >> 31;
    return x;
}

// Calls f(q) for each start q of a closed syncmer of x[0 .. len), in order.  The order g is a hash of the
// smaller of the codes of an s-mer and of its reverse complement (0 marks s-mers with other characters).
// Window minima by van Herk / Gil-Werman, in chunks.
template <class F> void for_each_syncmer(const uint8_t *x, uint64_t len, int K, int SM, F f) {
    static int8_t code[256];
    static bool init = false;
    if (!init) {
        for (int c = 0; c < 256; ++c) code[c] = -1;
        code['A'] = 0; code['C'] = 1; code['G'] = 2; code['T'] = 3;
        init = true;
    }
    if (len < (uint64_t)K) return;
    const uint64_t W = (uint64_t)(K - SM + 1);
    const uint64_t nq = len - K + 1;
    const uint64_t CH = W * 256;
    const uint64_t mask = SM == 32 ? ~0ull : (1ull << (2 * SM)) - 1;
    static thread_local std::vector<uint64_t> g, pre, suf;
    g.resize(CH + W); pre.resize(CH + W); suf.resize(CH + W);
    for (uint64_t c0 = 0; c0 < nq; c0 += CH) {
        const uint64_t cq = std::min(CH, nq - c0), ns = cq + W - 1;
        uint64_t smer = 0, rsmer = 0;
        int64_t bad = -1;
        for (int64_t i = 0; i < SM - 1; ++i) {
            int cc = code[x[c0 + i]];
            if (cc < 0) { bad = i; cc = 0; }
            smer = ((smer << 2) | (uint64_t)cc) & mask;
            rsmer = (rsmer >> 2) | ((uint64_t)(3 - cc) << (2 * (SM - 1)));
        }
        for (uint64_t t = 0; t < ns; ++t) {
            int64_t i = (int64_t)t + SM - 1;
            int cc = code[x[c0 + i]];
            if (cc < 0) { bad = i; cc = 0; }
            smer = ((smer << 2) | (uint64_t)cc) & mask;
            rsmer = (rsmer >> 2) | ((uint64_t)(3 - cc) << (2 * (SM - 1)));
            uint64_t h = mix64(smer < rsmer ? smer : rsmer);
            g[t] = bad >= (int64_t)t ? 0 : (h >> 1) + 1;
        }
        for (uint64_t b0 = 0; b0 < ns; b0 += W) {
            uint64_t b1 = std::min(ns, b0 + W), mn = ~0ull;
            for (uint64_t t = b0; t < b1; ++t) { mn = g[t] < mn ? g[t] : mn; pre[t] = mn; }
            mn = ~0ull;
            for (uint64_t t = b1; t-- > b0;) { mn = g[t] < mn ? g[t] : mn; suf[t] = mn; }
        }
        for (uint64_t q = 0; q < cq; ++q) {
            uint64_t a = suf[q], bb = pre[q + W - 1], mn = a < bb ? a : bb;
            if (mn && (mn == g[q] || mn == g[q + W - 1])) f(c0 + q);
        }
    }
}

// hash of a byte string (never 0)
inline uint64_t phrase_hash(const uint8_t *p, uint64_t len) {
    uint64_t hv = mix64(len ^ 0x9e3779b97f4a7c15ull), i = 0;
    for (; i + 8 <= len; i += 8) { uint64_t w; memcpy(&w, p + i, 8); hv = mix64(hv ^ w) + i; }
    uint64_t w = 0;
    memcpy(&w, p + i, len - i);
    hv = mix64(hv ^ w ^ (len << 56));
    return hv ? hv : 1;
}

// The phrase dictionary: the 64-bit hashes of the phrases of S over A/C/G/T, sorted and split into
// buckets by their top hb bits (about 8 per bucket); each slot stores the other 64 - hb bits (so the whole
// hash is checked), the phrase's id and the id of its reverse complement.
struct PhraseDict {
    uint32_t hb = 1;
    uint64_t size = 0;
    std::vector<uint32_t> bstart;
    PackedInts low, idf, idr;
    void build(std::vector<std::pair<uint64_t, std::pair<uint32_t, uint32_t>>> v) {   // (hash, (id, rc id))
        std::sort(v.begin(), v.end());
        for (size_t k = 1; k < v.size(); ++k)
            if (v[k].first == v[k - 1].first) throw std::runtime_error("two phrases with the same hash");
        size = v.size(); hb = 1;
        while (hb < 32 && (8ull << hb) < size) ++hb;
        uint64_t nb = 1ull << hb, k = 0;
        bstart.assign(nb + 1, 0);
        for (uint64_t bk = 0; bk < nb; ++bk) {
            bstart[bk] = (uint32_t)k;
            while (k < v.size() && (v[k].first >> (64 - hb)) == bk) ++k;
        }
        bstart[nb] = (uint32_t)k;
        std::vector<uint64_t> lw(size), f(size), r(size);
        for (size_t t = 0; t < size; ++t) { lw[t] = v[t].first & ((1ull << (64 - hb)) - 1); f[t] = v[t].second.first; r[t] = v[t].second.second; }
        low.build(lw); idf.build(f); idr.build(r);
    }
    inline void prefetch(uint64_t hv) const { __builtin_prefetch(&low.d[((uint64_t)bstart[hv >> (64 - hb)] * low.w) >> 6]); }
    inline uint64_t slot(uint64_t hv) const {        // slot of hv, or ~0
        uint64_t bk = hv >> (64 - hb), lo = hv & ((1ull << (64 - hb)) - 1);
        for (uint64_t a = bstart[bk], b = bstart[bk + 1]; a < b; ++a) {
            uint64_t x = low.get(a);
            if (x == lo) return a;
            if (x > lo) break;
        }
        return ~0ull;
    }
    uint64_t bytes() const { return bstart.size() * 4 + low.bytes() + idf.bytes() + idr.bytes(); }
    void save(FILE *f) const { wr(f, hb); wr(f, size); wvec(f, bstart); low.save(f); idf.save(f); idr.save(f); }
    void load(FILE *f) { rd(f, hb); rd(f, size); rvec(f, bstart); low.load(f); idf.load(f); idr.load(f); }
};

// The block-start rows: Elias-Fano or run-length coded, whichever is smaller.
struct BlockStarts {
    uint8_t rl = 0;
    uint64_t alt_bytes = 0;
    EliasFano ef;
    RLBits r;
    void build(const std::vector<uint64_t> &rows, uint64_t universe) {
        std::vector<std::pair<uint64_t, uint64_t>> runs;
        for (uint64_t x : rows) {
            if (!runs.empty() && runs.back().first + runs.back().second == x) ++runs.back().second;
            else runs.push_back({x, 1});
        }
        r.build(runs, universe);
        ef.build(rows, universe);
        rl = r.bytes() < ef.bytes();
        alt_bytes = rl ? ef.bytes() : r.bytes();
        if (rl) ef = EliasFano(); else r = RLBits();
    }
    inline uint64_t count_less(uint64_t x) const { return rl ? r.rank1(x) : ef.count_less(x); }
    inline uint64_t access(uint64_t k) const { return rl ? r.select1(k) : ef.access(k); }
    uint64_t bytes() const { return rl ? r.bytes() : ef.bytes(); }
    void save(FILE *f) const { wr(f, rl); if (rl) r.save(f); else ef.save(f); }
    void load(FILE *f) { rd(f, rl); if (rl) r.load(f); else ef.load(f); }
};

struct two_level {
    int32_t K = 51, SM = 12;
    uint64_t n = 0;          // rows of the character-level BWT
    uint64_t m = 0;          // phrases; the parse has m + 1 rows
    uint32_t sigma = 0;      // distinct phrases, ids 1..sigma
    uint64_t nblocks = 0;    // blocks; the block BWT has nblocks + 1 rows
    uint64_t nsym = 0;       // distinct block symbols, ids 1..nsym
    uint8_t has_blocks = 1;
    RLBits B;
    EFRank pbwt;
    BlockStarts BS;
    EFRank bbwt;
    PackedInts succp;               // succ(x), x = 0..sigma (0: none)
    EliasFano dstart;               // symbols of d are dstart[d] .. dstart[d+1]-1 (0-based)
    PackedInts jv;                  // j of each symbol, in symbol order within each d
    uint8_t sided = 0;              // 1: symbols have sides (sidev); 0: all +
    std::vector<uint64_t> sidev;
    std::vector<uint64_t> ccont;    // bit x: phrase x is a continuation somewhere in the parse
    PhraseDict dict;

    inline uint32_t succ(uint32_t x) const { return (uint32_t)succp.get(x); }
    inline bool can_continue(uint32_t x) const { return (ccont[x >> 6] >> (x & 63)) & 1; }
    inline bool sym_plus(uint64_t k) const { return !sided || ((sidev[k >> 6] >> (k & 63)) & 1); }
    static inline uint64_t sym_key(uint64_t j, bool plus) { return plus ? (1ull << 32) | (0xffffffffull - j) : j; }
    // the symbol (d, j, side) (1-based), or 0 if it does not occur
    inline uint32_t sym(uint32_t d, uint64_t j, bool plus) const {
        if (!sided && !plus) return 0;
        uint64_t a0, b0;
        dstart.access2(d, a0, b0);
        const uint64_t want = sym_key(j, plus);
        uint64_t a = a0, b = b0;
        while (a < b) {
            uint64_t mid = (a + b) >> 1;
            if (sym_key(jv.get(mid), sym_plus(mid)) < want) a = mid + 1; else b = mid;
        }
        return (a < b0 && sym_key(jv.get(a), sym_plus(a)) == want) ? (uint32_t)a + 1 : 0;
    }
    uint64_t bytes_parse() const { return B.bytes() + pbwt.bytes() + dict.bytes(); }
    uint64_t bytes_blocks() const { return succp.bytes() + ccont.size() * 8 + BS.bytes() + bbwt.bytes() + dstart.bytes() + jv.bytes() + sidev.size() * 8; }

    void save(const std::string &path) const {
        FILE *f = open_or_die(path, "wb");
        wmagic(f, "KTLVL001");
        wr(f, K); wr(f, SM); wr(f, n); wr(f, m); wr(f, sigma); wr(f, nblocks); wr(f, nsym); wr(f, has_blocks); wr(f, sided);
        B.save(f); pbwt.save(f); dict.save(f);
        if (has_blocks) { BS.save(f); bbwt.save(f); succp.save(f); dstart.save(f); jv.save(f); wvec(f, ccont); wvec(f, sidev); }
        fclose(f);
    }
    // parse_only: keep B and the dictionary (for STPD over the parse); the parse BWT and the blocks are read
    // and dropped
    void load(const std::string &path, bool parse_only = false) {
        FILE *f = open_or_die(path, "rb");
        rmagic(f, "KTLVL001", path);
        rd(f, K); rd(f, SM); rd(f, n); rd(f, m); rd(f, sigma); rd(f, nblocks); rd(f, nsym); rd(f, has_blocks); rd(f, sided);
        B.load(f); pbwt.load(f); dict.load(f);
        if (parse_only) { pbwt = EFRank(); has_blocks = 0; fclose(f); return; }
        if (has_blocks) { BS.load(f); bbwt.load(f); succp.load(f); dstart.load(f); jv.load(f); rvec(f, ccont); rvec(f, sidev); }
        fclose(f);
    }
};

// The parse of a pattern against the dictionary.
struct pattern_parse {
    std::vector<int64_t> S;         // syncmer starts
    std::vector<uint32_t> pid;      // pid[i]: id of P[S[i] .. S[i+1]+K), 0 if not a phrase (or i is last)
    std::vector<uint8_t> cont;      // cont[i]: phrase i is a continuation of phrase i-1
    std::vector<int32_t> bstart;    // bstart[i]: the block of phrase i starts at phrase bstart[i]
    void finish(const two_level &L) {
        size_t z = S.size();
        cont.assign(z, 0); bstart.assign(z, 0);
        for (size_t i = 0; i < z; ++i) {
            cont[i] = L.has_blocks && i > 0 && pid[i] && pid[i - 1] && L.succ(pid[i - 1]) == pid[i];
            bstart[i] = cont[i] ? bstart[i - 1] : (int32_t)i;
        }
    }
};

// Parses P (into pi) and rc(P) (into qi) with one scan for syncmers and one dictionary lookup per phrase.
inline void parse_pattern(const two_level &L, const char *Pc, uint64_t len, pattern_parse &pi, pattern_parse &qi) {
    const PhraseDict &D = L.dict;
    const uint8_t *x = reinterpret_cast<const uint8_t *>(Pc);
    const int K = L.K;
    pi.S.clear();
    for_each_syncmer(x, len, K, L.SM, [&](uint64_t q) { pi.S.push_back((int64_t)q); });
    const size_t z = pi.S.size(), np = z ? z - 1 : 0;
    const int64_t m = (int64_t)len;
    qi.S.resize(z);
    for (size_t i = 0; i < z; ++i) qi.S[i] = m - K - pi.S[z - 1 - i];
    static thread_local std::vector<uint64_t> h, sl;
    h.resize(np); sl.resize(np);
    for (size_t i = 0; i < np; ++i) h[i] = phrase_hash(x + pi.S[i], pi.S[i + 1] + K - pi.S[i]);
    constexpr size_t AHEAD = 16;
    for (size_t i = 0; i < np && i < AHEAD; ++i) D.prefetch(h[i]);
    for (size_t i = 0; i < np; ++i) {
        if (i + AHEAD < np) D.prefetch(h[i + AHEAD]);
        sl[i] = D.slot(h[i]);
    }
    pi.pid.assign(z, 0); qi.pid.assign(z, 0);
    for (size_t i = 0; i < np; ++i) {
        if (sl[i] == ~0ull) continue;
        pi.pid[i] = (uint32_t)D.idf.get(sl[i]);
        qi.pid[np - 1 - i] = (uint32_t)D.idr.get(sl[i]);
    }
    pi.finish(L); qi.finish(L);
}

struct search_counters {
    u64 searches = 0, char_steps = 0, parse_steps = 0, block_steps = 0, kmer_tests = 0, kmer_restarts = 0;
    u64 stpd_jumps = 0, stpd_compares = 0;     // STPD: binary searches over the sample, and text comparisons
    // KeBaB with --longest (mems.hpp): pseudo-MEMs in the reads, pseudo-MEMs searched (BML run on them),
    // batches, reads that needed a second batch, reads that ran the whole-read fallback
    u64 kb_pseudo = 0, kb_searched = 0, kb_batches = 0, kb_multibatch = 0, kb_fallbacks = 0;
    // STPD over the parse: character steps at the right end (backward search), phrase comparisons in the
    // middle (one per pattern phrase), parse-sample binary searches and their comparisons, and character
    // steps of the character-level STPD at the left end
    u64 ps_right = 0, ps_mid = 0, ps_jumps = 0, ps_compares = 0, ps_left = 0;
    void add(const search_counters &o) {
        ps_right += o.ps_right; ps_mid += o.ps_mid; ps_jumps += o.ps_jumps; ps_compares += o.ps_compares; ps_left += o.ps_left;
        searches += o.searches; char_steps += o.char_steps; parse_steps += o.parse_steps; block_steps += o.block_steps;
        kmer_tests += o.kmer_tests; kmer_restarts += o.kmer_restarts; stpd_jumps += o.stpd_jumps; stpd_compares += o.stpd_compares;
        kb_pseudo += o.kb_pseudo; kb_searched += o.kb_searched; kb_batches += o.kb_batches; kb_multibatch += o.kb_multibatch;
        kb_fallbacks += o.kb_fallbacks;
    }
};

// The longest suffix X[j..end) of X[begin..end) that occurs in S, and its interval: the same result as
// plain backward search (Theorem "twolevel"; with jumps, Lemmas "blockstart" and "backup").  Returns j.
inline u64 search2(const stepper &st, const two_level &L, const pattern_parse &pi, const char *X, int64_t begin, int64_t end,
                   bool jumps, interval &out, search_counters &cnt) {
    const int K = L.K;
    interval I = st.full();
    int64_t j = end;
    auto char_steps_to = [&](int64_t stop) -> bool {      // false if a step fails
        while (j > stop) {
            ++cnt.char_steps;
            if (!st.step(I, (unsigned char)X[j - 1])) return false;
            --j;
        }
        return true;
    };
    ++cnt.searches;
    const auto &S = pi.S;
    int64_t idx = (int64_t)(std::upper_bound(S.begin(), S.end(), end - K) - S.begin()) - 1;
    if (idx >= 0 && S[idx] >= begin && char_steps_to(S[idx])) {
        // parse rows [plo, phi] (inclusive) for X[S[i] .. end)
        uint64_t plo = L.B.rank1(I.sp) + 1, phi = L.B.rank1(I.ep);
        int64_t i = idx, noblock_above = idx - 1;
        jumps = jumps && L.has_blocks;
        while (i > 0 && S[i - 1] >= begin && pi.pid[i - 1]) {
            if (jumps && i <= noblock_above && !pi.cont[i]) {
                uint64_t blo = L.BS.count_less(plo), bhi = L.BS.count_less(phi + 1);
                int64_t cur = i, taken = 0;
                uint64_t pblo = 0, pbhi = 0; int64_t pcur = i;
                while (cur > 0 && blo < bhi) {
                    int32_t t0 = pi.bstart[cur - 1];
                    if (S[t0] < begin || !pi.pid[t0]) break;
                    const bool plus = pi.pid[cur] > L.succ(pi.pid[cur - 1]);
                    uint32_t y = L.sym(pi.pid[t0], (uint64_t)(cur - 1 - t0), plus);
                    if (!y) break;
                    uint64_t o1, o2;
                    L.bbwt.lf2(y, blo, bhi, o1, o2);
                    ++cnt.block_steps;
                    if (o1 >= o2) break;
                    pblo = blo; pbhi = bhi; pcur = cur;
                    blo = o1; bhi = o2; cur = t0; ++taken;
                }
                if (taken >= 1 && !L.can_continue(pi.pid[cur])) {     // exact (Lemma "backup", first case)
                    plo = L.BS.access(blo); phi = L.BS.access(bhi - 1);
                    i = cur;
                    noblock_above = i;
                    continue;
                }
                if (taken >= 2) {             // back up to the interval before the last block
                    plo = L.BS.access(pblo); phi = L.BS.access(pbhi - 1);
                    i = pcur;
                }
                noblock_above = pi.bstart[i - 1];
            }
            uint64_t o1, o2;
            L.pbwt.lf2(pi.pid[i - 1], plo, phi + 1, o1, o2);
            ++cnt.parse_steps;
            if (o1 >= o2) break;
            plo = o1; phi = o2 - 1; --i;
        }
        I.sp = L.B.select1(plo - 1); I.ep = L.B.select1(phi - 1) + 1;
        j = S[i];
        char_steps_to(begin);
    } else if (idx < 0 || S[idx] < begin) {
        char_steps_to(begin);
    }
    out = I;
    return (u64)j;
}

}  // namespace katka
