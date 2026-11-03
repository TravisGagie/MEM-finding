// A text oracle for S$ that reads it through the two-level parse instead of the plain text: the phrase
// starts (the boundary list, Elias-Fano), the parse (phrase ids, packed) and the phrase dictionary (the
// distinct phrases concatenated, in id order: the two-level index's lexicographic ids).  A position maps to
// its phrase by a predecessor query on the starts; a per-thread cursor keeps the current phrase, so reads
// in or next to it (the STPD comparisons read rightward from a sampled suffix, the LF-like steps leftward)
// cost no query.  Consecutive phrases overlap by k characters, so a position can be read from either.
//
// Space: about (2 + log(n/m)) + ceil(log2(sigma + 1)) bits per phrase plus the dictionary, against 8 bits
// per character for the plain text.
#pragma once
#include "twolevel_build.hpp"
#include <atomic>
#include "libsais.h"

namespace katka {

// The phrase dictionary at 2 bits per base: the phrases concatenated in id order, with A/C/G/T packed and
// any other character stored as A in the packed array; a bit per phrase marks those with other characters,
// and for each marked phrase its runs of X as (offset, length) and the offset of $ (byte 1; the length if
// none).  The text has no other characters (katka prep writes X for every non-A/C/G/T character).
// Decoding is exact, so comparisons with X and $ give the same order as with the bytes.
struct packed_dict {
    u64 total = 0;                    // bases
    std::vector<uint64_t> w;          // 32 bases per word, base i at bits 2(i mod 32)
    BitRank marked;                   // bit id (1..sigma)
    PackedInts roff, rpos, rlen, dpos;
    u64 nmarked = 0, nruns = 0;

    void build(const std::vector<uint8_t> &d, const std::vector<uint64_t> &doffv, uint32_t sigma) {
        static const int8_t code[256] = {
            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
            -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3};
        total = d.size();
        w.assign(total / 32 + 1, 0);
        marked.init(sigma + 1);
        std::vector<uint64_t> ro{0}, rp, rl, dp;
        for (uint32_t id = 1; id <= sigma; ++id) {
            const u64 a = doffv[id - 1], b = doffv[id];
            bool other = false;
            u64 dol = b - a;
            for (u64 i = a; i < b; ++i) {
                const uint8_t c = d[i];
                int x = c < 85 ? code[c] : -1;
                if (x >= 0) { w[i >> 5] |= (uint64_t)x << (2 * (i & 31)); continue; }
                other = true;
                if (c == 1) { dol = i - a; continue; }
                if (c != 'X') die("parse oracle: a character other than A, C, G, T, X and $ in the text");
                if (i > a && d[i - 1] == 'X') rl.back()++;
                else { rp.push_back(i - a); rl.push_back(1); }
            }
            if (other) { marked.set(id); dp.push_back(dol); ro.push_back(rp.size()); }
        }
        marked.finish();
        nmarked = dp.size(); nruns = rp.size();
        roff.build(ro); rpos.build(rp.empty() ? std::vector<uint64_t>{0} : rp); rlen.build(rl.empty() ? std::vector<uint64_t>{0} : rl);
        dpos.build(dp.empty() ? std::vector<uint64_t>{0} : dp);
    }
    // characters [from, to) of phrase id (which starts at base a) into out[0 .. to - from)
    void decode(u64 id, u64 a, u64 from, u64 to, uint8_t *out) const {
        // 4 bases per byte of the packed array, from a table of their 4 characters
        static const struct table { uint32_t t[256]; table() {
            static const uint8_t base[4] = {'A', 'C', 'G', 'T'};
            for (int b = 0; b < 256; ++b) {
                uint8_t c[4];
                for (int q = 0; q < 4; ++q) c[q] = base[(b >> (2 * q)) & 3];
                memcpy(&t[b], c, 4);
            }
        } } T;
        static const uint8_t base[4] = {'A', 'C', 'G', 'T'};
        const uint8_t *bytes = reinterpret_cast<const uint8_t *>(w.data());
        const u64 len = to - from;
        u64 i = a + from, o = 0;
        while (o < len && (i & 3)) { out[o++] = base[(bytes[i >> 2] >> (2 * (i & 3))) & 3]; ++i; }
        while (o + 4 <= len) { memcpy(out + o, &T.t[bytes[i >> 2]], 4); o += 4; i += 4; }
        while (o < len) { out[o++] = base[(bytes[i >> 2] >> (2 * (i & 3))) & 3]; ++i; }
        if (marked.get(id)) {                    // the X runs overlapping [from, to) (sorted, disjoint), and $
            const u64 r = marked.rank1(id), d = dpos.get(r), q1 = roff.get(r + 1);
            u64 lo = roff.get(r), hi = q1;
            while (lo < hi) { u64 mid = (lo + hi) / 2; if (rpos.get(mid) + rlen.get(mid) <= from) lo = mid + 1; else hi = mid; }
            for (u64 q = lo; q < q1; ++q) {
                const u64 p0 = rpos.get(q), p1 = p0 + rlen.get(q);
                if (p0 >= to) break;
                const u64 x0 = std::max(p0, from), x1 = std::min(p1, to);
                memset(out + (x0 - from), 'X', x1 - x0);
            }
            if (d >= from && d < to) out[d - from] = 1;
        }
    }
    u64 bytes() const { return w.size() * 8 + marked.bytes() + roff.bytes() + rpos.bytes() + rlen.bytes() + dpos.bytes(); }
    void save(FILE *f) const { wr(f, total); wr(f, nmarked); wr(f, nruns); wvec(f, w); marked.save(f); roff.save(f); rpos.save(f); rlen.save(f); dpos.save(f); }
    void load(FILE *f) { rd(f, total); rd(f, nmarked); rd(f, nruns); rvec(f, w); marked.load(f); roff.load(f); rpos.load(f); rlen.load(f); dpos.load(f); }
};

// The parse as an RLZAP-style factorization against a reference parse (Cox, Farruggia, Gagie, Puglisi and
// Sirén, "RLZAP: Relative Lempel-Ziv with Adaptive Pointers"): the first R phrases (the reference: in katka's
// texts, the first sequence and its reverse complement) are stored as packed ids with their starts; the
// rest is a sequence of factors, each some literal phrase ids followed by a copy of the reference parse at a
// relative pointer (source = target - rel), with the pointers run-length coded (a value per change).
// Factorization, as tool/results/chr19_rlz/rlz_parse.cpp (whose sizes it reproduces): after each copy, look
// ahead j = 0..J phrases for a match of at least M phrases at rel + delta, |delta| <= D (smallest j, then
// smallest |delta|; the j phrases become literals); otherwise the greedy longest match (among equally long
// ones, the one continuing the pointer), or a one-phrase literal if the phrase is not in the reference.
// Access: phrase id, text start of a phrase, and the phrase at a text position, through the factor that
// holds the position (Elias-Fano over factor starts in parse and in text positions), with a per-thread cache.
struct rlz_parse {
    u64 m = 0, R = 0, n = 0;
    uint32_t K = 0;
    PackedInts ref;                   // Z[0..R)
    EliasFano refstart;               // the starts of phrases 0..R-1
    EliasFano fpos, ftext;            // factor f: its first phrase, and that phrase's text start
    EliasFano litcum;                 // literal phrases before factor f (z + 1 values)
    PackedInts lits;                  // the literal ids
    EliasFano chg; PackedInts relv;   // the factors where the pointer changes, and -rel there
    u64 z = 0, nchg = 0;

    struct fcache { const rlz_parse *owner = nullptr; u64 f = 0, p0 = 0, p1 = 0, t0 = 0, nl = 0, lb = 0, rel = 0; };
    static fcache &fc() { static thread_local fcache c; return c; }
    inline const fcache &factor_at(u64 i) const {          // i >= R
        fcache &c = fc();
        if (c.owner == this && i >= c.p0 && i < c.p1) return c;
        c.owner = this;
        u64 f = fpos.pred_index(i);
        c.f = f;
        if (f + 1 < z) fpos.access2(f, c.p0, c.p1); else { c.p0 = fpos.access(f); c.p1 = m; }
        c.t0 = ftext.access(f);
        u64 a, b; litcum.access2(f, a, b); c.lb = a; c.nl = b - a;
        u64 k = chg.count_less(f + 1);
        c.rel = k ? relv.get(k - 1) : 0;
        return c;
    }
    inline u64 zid(u64 i) const {
        if (i < R) return ref.get(i);
        const fcache &c = factor_at(i);
        u64 o = i - c.p0;
        return o < c.nl ? lits.get(c.lb + o) : ref.get(i - c.rel);
    }
    template <class LEN> inline u64 start(u64 i, LEN len) const {
        if (i < R) return refstart.access(i);
        const fcache &c = factor_at(i);
        u64 o = i - c.p0, s = c.t0, q = 0;
        for (; q < o && q < c.nl; ++q) s += len(lits.get(c.lb + q)) - K;
        if (o <= c.nl) return s;
        const u64 c0 = c.p0 + c.nl;                          // the copy: starts as in the reference
        return s + refstart.access(i - c.rel) - refstart.access(c0 - c.rel);
    }
    template <class LEN> inline u64 pred(u64 x, LEN len) const {    // the last phrase starting at or before x
        if (R == m || x < ftext.access(0)) return refstart.pred_index(x);
        const u64 f = ftext.pred_index(x);
        const fcache &c = factor_at(fpos.access(f));
        u64 s = c.t0;
        for (u64 q = 0; q < c.nl; ++q) {
            u64 nx = s + len(lits.get(c.lb + q)) - K;
            if (x < nx || c.p0 + q + 1 == c.p1) return c.p0 + q;
            s = nx;
        }
        const u64 c0 = c.p0 + c.nl, r0 = refstart.access(c0 - c.rel);
        return refstart.pred_index(x - s + r0) + c.rel;
    }
    u64 bytes() const { return ref.bytes() + refstart.bytes() + fpos.bytes() + ftext.bytes() + litcum.bytes() + lits.bytes() + chg.bytes() + relv.bytes(); }
    u64 bytes_ref() const { return ref.bytes() + refstart.bytes(); }
    void save(FILE *f) const { wr(f, m); wr(f, R); wr(f, n); wr(f, K); wr(f, z); wr(f, nchg); ref.save(f); refstart.save(f); fpos.save(f); ftext.save(f); litcum.save(f); lits.save(f); chg.save(f); relv.save(f); }
    void load(FILE *f) { rd(f, m); rd(f, R); rd(f, n); rd(f, K); rd(f, z); rd(f, nchg); ref.load(f); refstart.load(f); fpos.load(f); ftext.load(f); litcum.load(f); lits.load(f); chg.load(f); relv.load(f); }
};

struct parse_text {
    u64 n = 0, m = 0;                 // |S$|; phrases
    uint32_t K = 0, SM = 0, sigma = 0;
    packed_dict dict;                 // phrase a (1..sigma) is bases doff[a-1] .. doff[a] of dict
    PackedInts doff, ids;             // ids[i]: the id of phrase i
    bool rlz = false;                 // the parse as an RLZAP-style factorization (Zr) instead of ids and starts
    rlz_parse Zr;
    // the parse: the id of phrase i, its start in S$, and the phrase at text position x
    inline u64 zid(u64 i) const { return rlz ? Zr.zid(i) : ids.get(i); }
    inline u64 start(u64 i) const { return rlz ? Zr.start(i, [&](u64 id) { return len_of(id); }) : starts.access(i); }
    inline u64 pred(u64 x) const { return rlz ? Zr.pred(x, [&](u64 id) { return len_of(id); }) : starts.pred_index(x); }
    u64 bytes_parse() const { return rlz ? Zr.bytes() : ids.bytes() + starts.bytes(); }
    // characters [from, to) of phrase id into buf (resized as needed)
    const uint8_t *phrase(u64 id, u64 from, u64 to, std::vector<uint8_t> &buf) const {
        u64 a = doff.get(id - 1);
        if (buf.size() < to - from) buf.resize(to - from);
        dict.decode(id, a, from, to, buf.data());
        return buf.data();
    }
    EliasFano starts;                 // the start of each phrase in S$

    // statistics: seeks served by the current phrase, by a neighbour, by a predecessor query (summed over
    // threads as they exit)
    static std::atomic<u64> *totals() { static std::atomic<u64> t[3] = {{0}, {0}, {0}}; return t; }
    struct cursor {                    // phrase i is S$[s..e)
        const parse_text *owner = nullptr; u64 i = 0, s = 0, e = 0; const uint8_t *p = nullptr;
        u64 id = 0, a = 0;             // its id and first base in the dictionary
        u64 ws = 0, we = 0;            // the decoded window S$[ws..we) of the phrase: p[x - ws]
        std::vector<uint8_t> buf;
        u64 same = 0, step = 0, query = 0;
        ~cursor() { totals()[0] += same; totals()[1] += step; totals()[2] += query; }
    };
    static cursor &cur() { static thread_local cursor c; return c; }

    inline u64 len_of(u64 id) const { return doff.get(id) - doff.get(id - 1); }
    inline void set_phrase(cursor &c, u64 i, u64 s) const {
        u64 id = zid(i), a = doff.get(id - 1);
        c.i = i; c.s = s; c.e = s + (doff.get(id) - a); c.id = id; c.a = a; c.ws = c.we = 0;
    }
    // makes the decoded window contain x (in the cursor's phrase): 64 characters to its left, 192 to its right
    static constexpr u64 WL = 64, WR = 192;
    inline void window(cursor &c, u64 x) const {
        if (x >= c.ws && x < c.we) return;
        const u64 o = x - c.s, from = o > WL ? o - WL : 0, to = std::min(c.e - c.s, o + WR);
        if (c.buf.size() < WL + WR) c.buf.resize(WL + WR);
        dict.decode(c.id, c.a, from, to, c.buf.data());
        c.ws = c.s + from; c.we = c.s + to; c.p = c.buf.data();
    }
    // makes the cursor's phrase contain x
    inline void seek(cursor &c, u64 x) const {
        if (c.owner == this) {
            if (x >= c.s && x < c.e) { ++c.same; return; }
            if (x >= c.e && c.i + 1 < m) {                   // the next phrase starts at e - K
                set_phrase(c, c.i + 1, c.e - K);
                if (x < c.e) { ++c.step; return; }
            } else if (x < c.s && c.i > 0) {                 // the previous phrase ends at s + K
                u64 e = c.s + K, l = len_of(zid(c.i - 1));
                set_phrase(c, c.i - 1, e - l);
                if (x >= c.s) { ++c.step; return; }
            }
        }
        ++c.query;
        c.owner = this;
        u64 i = pred(x);
        set_phrase(c, i, start(i));
    }
    inline uint8_t at(u64 x) const { cursor &c = cur(); seek(c, x); window(c, x); return c.p[x - c.ws]; }
    // the longest common prefix of S$[b..b+len) and x[0..len) (b + len <= n)
    inline u64 lce(u64 b, const uint8_t *x, u64 len) const {
        cursor &c = cur();
        u64 l = 0;
        while (l < len) {
            seek(c, b + l); window(c, b + l);
            u64 avail = std::min(len - l, c.we - (b + l));
            u64 q = lce_bytes(c.p + (b + l - c.ws), x + l, avail);
            l += q;
            if (q < avail) break;
        }
        return l;
    }
    // the longest common prefix of S$[a..a+len) and S$[b..b+len) (with two cursors of its own)
    u64 lce2(u64 a, u64 b, u64 len) const {
        cursor ca, cb;
        u64 l = 0;
        while (l < len) {
            seek(ca, a + l); window(ca, a + l); seek(cb, b + l); window(cb, b + l);
            u64 avail = std::min({len - l, ca.we - (a + l), cb.we - (b + l)});
            u64 q = lce_bytes(ca.p + (a + l - ca.ws), cb.p + (b + l - cb.ws), avail);
            l += q;
            if (q < avail) break;
        }
        return l;
    }
    u64 bytes() const { return dict.bytes() + doff.bytes() + bytes_parse(); }

    // The same parse and ids as build_twolevel (syncmer_parse with (k, s), lexicographic ids), from S only.
    // The parse is spilled to tmp as phrase indices while S is scanned.
    void build(const std::string &Sfile, int k, int s, const std::string &tmp) {
        double t0 = now();
        if (s < 1 || s >= k || s > 32) die("parse oracle: need 1 <= s < k and s <= 32");
        K = k; SM = s;
        std::vector<uint8_t> d;
        std::vector<uint64_t> off, rhash;
        std::unordered_map<uint64_t, uint32_t> seen;
        FILE *ft = open_or_die(tmp, "w+b");
        std::vector<uint32_t> pbuf;
        pbuf.reserve(1 << 20);
        m = 0;
        bool zero;
        n = syncmer_parse(Sfile, k, s, zero, [&](const uint8_t *p, u64 len) {
            u64 hv = phrase_hash(p, len);
            auto it = seen.find(hv);
            uint32_t idx;
            if (it == seen.end()) {
                idx = (uint32_t)rhash.size();
                seen.emplace(hv, idx); rhash.push_back(hv); off.push_back(d.size());
                d.insert(d.end(), p, p + len);
            } else {
                idx = it->second;
                u64 a = off[idx], b = idx + 1 < off.size() ? off[idx + 1] : d.size();
                if (b - a != len || memcmp(d.data() + a, p, len)) throw std::runtime_error("phrase hash collision");
            }
            pbuf.push_back(idx); ++m;
            if (pbuf.size() == pbuf.capacity()) { fwrite(pbuf.data(), 4, pbuf.size(), ft); pbuf.clear(); }
        });
        if (!pbuf.empty()) fwrite(pbuf.data(), 4, pbuf.size(), ft);
        std::vector<uint32_t>().swap(pbuf);
        seen = std::unordered_map<uint64_t, uint32_t>();
        std::vector<uint64_t>().swap(rhash);
        off.push_back(d.size());
        sigma = (uint32_t)(off.size() - 1);
        // lexicographic ids (as in build_twolevel), and the dictionary in id order
        std::vector<uint32_t> order(sigma), idof(sigma);
        for (uint32_t a = 0; a < sigma; ++a) order[a] = a;
        std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t c) {
            u64 la = off[a + 1] - off[a], lc = off[c + 1] - off[c];
            int r = memcmp(d.data() + off[a], d.data() + off[c], std::min(la, lc));
            if (r == 0 && la != lc) throw std::runtime_error("dictionary is not prefix-free");
            return r < 0;
        });
        std::vector<uint64_t> doffv(sigma + 1, 0);
        {
            std::vector<uint8_t> dd;
            dd.reserve(d.size());
            for (uint32_t r = 0; r < sigma; ++r) {
                uint32_t a = order[r];
                idof[a] = r + 1;
                dd.insert(dd.end(), d.begin() + off[a], d.begin() + off[a + 1]);
                doffv[r + 1] = dd.size();
            }
            std::vector<uint8_t>().swap(d);
            dict.build(dd, doffv, sigma);
        }
        std::vector<uint32_t>().swap(order);
        doff.build(doffv);
        // the parse (ids) and the phrase starts, in one pass over tmp
        ids.init(m, PackedInts::width_for(sigma));
        starts.begin(m, n);
        rewind(ft);
        std::vector<uint32_t> buf(1 << 20);
        u64 i = 0, st = 0;
        while (i < m) {
            size_t want = std::min<u64>(buf.size(), m - i);
            if (fread(buf.data(), 4, want, ft) != want) die("cannot read " + tmp);
            for (size_t q = 0; q < want; ++q, ++i) {
                uint32_t id = idof[buf[q]];
                ids.set(i, id);
                starts.push(st);
                u64 len = doffv[id] - doffv[id - 1];
                if (i + 1 < m) st += len - K; else st += len;
            }
        }
        fclose(ft); remove(tmp.c_str());
        starts.finish();
        if (st != n) die("parse oracle: the phrases do not cover S$");
        fprintf(stderr, "parse oracle: n=%lu, %lu phrases, %u distinct; dictionary %.1f MB (%lu bases at 2 bits, %lu phrases with other characters, %lu X runs; offsets %.1f MB), parse %.1f MB, starts %.1f MB (%.2f bits per character in all; %.1f s)\n",
                n, m, sigma, dict.bytes() / 1e6, dict.total, dict.nmarked, dict.nruns, doff.bytes() / 1e6, ids.bytes() / 1e6, starts.bytes() / 1e6, 8.0 * bytes() / n, now() - t0);
    }
    // Replaces the parse (ids, starts) by an RLZAP-style factorization against its first R phrases, those
    // starting before refend (the reference); prints the counts that results/chr19_rlz/rlz_parse.cpp prints.
    void make_rlz(u64 refend, u64 J, int64_t D, u64 M) {
        double t0 = now();
        if (rlz) return;
        std::vector<int32_t> id(m);
        for (u64 i = 0; i < m; ++i) id[i] = (int32_t)ids.get(i);
        u64 R = 0;
        while (R < m && starts.access(R) < refend) ++R;
        if (R == 0) die("RLZ: an empty reference");
        if (R >= (1ull << 31)) die("RLZ: the reference parse is too long");
        std::vector<int32_t> T(id.begin(), id.begin() + R), SA(R), ISA(R);
        if (libsais_int(T.data(), SA.data(), (int32_t)R, (int32_t)sigma + 1, 0) != 0) die("libsais_int failed");
        for (u64 i = 0; i < R; ++i) T[i] = id[i];               // (libsais_int overwrote it)
        for (u64 x = 0; x < R; ++x) ISA[SA[x]] = (int32_t)x;
        auto at = [&](int64_t q) -> int64_t { return q < (int64_t)R ? T[q] : -1; };
        auto match = [&](u64 p, u64 &lo, u64 &hi) {           // the longest match of id[p..) in T, its SA interval
            lo = 0; hi = R; u64 d = 0;
            while (p + d < m) {
                int64_t c = id[p + d];
                u64 a = lo, b = hi;
                while (a < b) { u64 mid = (a + b) / 2; if (at((int64_t)SA[mid] + d) < c) a = mid + 1; else b = mid; }
                u64 l2 = a; b = hi;
                while (a < b) { u64 mid = (a + b) / 2; if (at((int64_t)SA[mid] + d) <= c) a = mid + 1; else b = mid; }
                if (l2 == a) break;
                lo = l2; hi = a; ++d;
            }
            return d;
        };
        auto lce = [&](u64 p, int64_t q) { u64 l = 0; while (p + l < m && q + (int64_t)l < (int64_t)R && id[p + l] == T[q + l]) ++l; return l; };
        std::vector<uint64_t> fp, ft, lc{0}, lv, cf, cv;
        u64 nlit = 0, litonly = 0, changes = 0, explicit_ = 0;
        bool have = false; int64_t rel = 0;                    // source = target + rel
        int64_t stored = INT64_MIN;
        auto factor = [&](u64 p, u64 j, u64 copylen, int64_t r) {
            fp.push_back(p); ft.push_back(starts.access(p));
            for (u64 q = 0; q < j; ++q) lv.push_back((uint64_t)id[p + q]);
            nlit += j; lc.push_back(nlit);
            if (copylen && r != stored) { cf.push_back(fp.size() - 1); cv.push_back((uint64_t)(-r)); stored = r; }
        };
        u64 p = R;
        while (p < m) {
            bool done = false;
            if (have) {
                for (u64 j = 0; j <= J && !done && p + j < m; ++j)
                    for (int64_t a = 0; a <= D && !done; ++a)
                        for (int sg = 0; sg < (a ? 2 : 1) && !done; ++sg) {
                            int64_t dl = sg ? -a : a, src = (int64_t)(p + j) + rel + dl;
                            if (src < 0 || src >= (int64_t)R) continue;
                            u64 l = lce(p + j, src);
                            if (l >= M) {
                                if (dl != 0) { ++changes; rel += dl; }
                                factor(p, j, l, rel);
                                p += j + l; done = true;
                            }
                        }
            }
            if (done) continue;
            u64 lo, hi, d = match(p, lo, hi);
            if (d == 0) { ++litonly; factor(p, 1, 0, rel); p += 1; continue; }
            int64_t q = SA[lo];
            if (have) { int64_t c = (int64_t)p + rel; if (c >= 0 && c < (int64_t)R && (u64)ISA[c] >= lo && (u64)ISA[c] < hi) q = c; }
            int64_t nr = q - (int64_t)p;
            if (!have || nr != rel) { ++explicit_; rel = nr; have = true; }
            factor(p, 0, d, rel);
            p += d;
        }
        Zr = rlz_parse();
        Zr.m = m; Zr.R = R; Zr.n = n; Zr.K = K; Zr.z = fp.size(); Zr.nchg = cf.size();
        {
            std::vector<uint64_t> rv(T.begin(), T.end()), rs(R);
            for (u64 i = 0; i < R; ++i) rs[i] = starts.access(i);
            Zr.ref.build(rv); Zr.refstart.build(rs, n);
        }
        if (fp.empty()) { fp.push_back(m); ft.push_back(n); lc.push_back(0); }   // (no target: a dummy factor)
        Zr.fpos.build(fp, m + 1); Zr.ftext.build(ft, n + 1); Zr.litcum.build(lc, nlit + 1);
        Zr.lits.build(lv.empty() ? std::vector<uint64_t>{0} : lv);
        if (cf.empty()) { cf.push_back(0); cv.push_back(0); }
        Zr.chg.build(cf, fp.size() + 1); Zr.relv.build(cv);
        rlz = true;
        // check every phrase against the plain parse before dropping it
        for (u64 i = 0; i < m; ++i)
            if (zid(i) != ids.get(i) || start(i) != starts.access(i)) die("RLZ parse: phrase " + std::to_string(i) + " differs");
        for (u64 i = 0; i < m; i += 1 + i % 7) { u64 x = starts.access(i) + (i % 5); if (x < n && pred(x) != starts.pred_index(x)) die("RLZ parse: pred differs"); }
        ids = PackedInts(); starts = EliasFano();
        fprintf(stderr, "RLZAP parse (J = %lu, D = %ld, M = %lu): reference %lu phrases; target %lu phrases in %lu factors (%.2f phrases per factor); %lu literal phrases (%lu factors only a literal); %lu pointer changes, %lu explicit pointers; %lu stored pointer values\n",
                J, D, M, R, m - R, (u64)Zr.z, (double)(m - R) / std::max<u64>(1, Zr.z), nlit, litonly, changes, explicit_, (u64)Zr.nchg);
        fprintf(stderr, "  sizes (MB): reference ids %.3f + starts %.3f; factor starts (parse %.3f, text %.3f); literal counts %.3f; literals %.3f; pointers (change positions %.3f, values %.3f): parse part %.3f; dictionary %.3f (+ offsets %.3f) (%.1f s)\n",
                Zr.ref.bytes() / 1e6, Zr.refstart.bytes() / 1e6, Zr.fpos.bytes() / 1e6, Zr.ftext.bytes() / 1e6, Zr.litcum.bytes() / 1e6, Zr.lits.bytes() / 1e6,
                Zr.chg.bytes() / 1e6, Zr.relv.bytes() / 1e6, Zr.bytes() / 1e6, dict.bytes() / 1e6, doff.bytes() / 1e6, now() - t0);
    }
    void save(const std::string &f) const {
        if (rlz) {
            FILE *o = open_or_die(f, "wb"); wmagic(o, "KTPTXR01");
            wr(o, n); wr(o, m); wr(o, K); wr(o, SM); wr(o, sigma);
            dict.save(o); doff.save(o); Zr.save(o);
            fclose(o);
            return;
        }
        FILE *o = open_or_die(f, "wb"); wmagic(o, "KTPTXT02");
        wr(o, n); wr(o, m); wr(o, K); wr(o, SM); wr(o, sigma);
        dict.save(o); doff.save(o); ids.save(o); starts.save(o);
        fclose(o);
    }
    static bool current_format(const std::string &f) {
        FILE *i = fopen(f.c_str(), "rb");
        if (!i) return false;
        char m8[8] = {0};
        bool ok = fread(m8, 1, 8, i) == 8 && !memcmp(m8, "KTPTXT02", 8);
        fclose(i);
        return ok;
    }
    void load(const std::string &f) {
        if (FILE *t = fopen(f.c_str(), "rb")) {                 // the RLZ format
            char m8[8] = {0};
            bool r = fread(m8, 1, 8, t) == 8 && !memcmp(m8, "KTPTXR01", 8);
            if (r) {
                rd(t, n); rd(t, m); rd(t, K); rd(t, SM); rd(t, sigma);
                dict.load(t); doff.load(t); Zr.load(t); rlz = true;
                fclose(t);
                return;
            }
            fclose(t);
        }
        FILE *i = open_or_die(f, "rb"); rmagic(i, "KTPTXT02", f);
        rd(i, n); rd(i, m); rd(i, K); rd(i, SM); rd(i, sigma);
        dict.load(i); doff.load(i); ids.load(i); starts.load(i);
        fclose(i);
    }
};

}  // namespace katka
