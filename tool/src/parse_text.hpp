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

namespace katka {

struct parse_text {
    u64 n = 0, m = 0;                 // |S$|; phrases
    uint32_t K = 0, SM = 0, sigma = 0;
    std::vector<uint8_t> dict;        // phrase a (1..sigma) is dict[doff[a-1] .. doff[a])
    PackedInts doff, ids;             // ids[i]: the id of phrase i
    EliasFano starts;                 // the start of each phrase in S$

    // statistics: seeks served by the current phrase, by a neighbour, by a predecessor query (summed over
    // threads as they exit)
    static std::atomic<u64> *totals() { static std::atomic<u64> t[3] = {{0}, {0}, {0}}; return t; }
    struct cursor {                    // phrase i is S$[s..e)
        const parse_text *owner = nullptr; u64 i = 0, s = 0, e = 0; const uint8_t *p = nullptr;
        u64 same = 0, step = 0, query = 0;
        ~cursor() { totals()[0] += same; totals()[1] += step; totals()[2] += query; }
    };
    static cursor &cur() { static thread_local cursor c; return c; }

    inline u64 len_of(u64 id) const { return doff.get(id) - doff.get(id - 1); }
    inline void set_phrase(cursor &c, u64 i, u64 s) const {
        u64 id = ids.get(i), a = doff.get(id - 1);
        c.i = i; c.s = s; c.e = s + (doff.get(id) - a); c.p = dict.data() + a;
    }
    // makes the cursor's phrase contain x
    inline void seek(cursor &c, u64 x) const {
        if (c.owner == this) {
            if (x >= c.s && x < c.e) { ++c.same; return; }
            if (x >= c.e && c.i + 1 < m) {                   // the next phrase starts at e - K
                set_phrase(c, c.i + 1, c.e - K);
                if (x < c.e) { ++c.step; return; }
            } else if (x < c.s && c.i > 0) {                 // the previous phrase ends at s + K
                u64 e = c.s + K, l = len_of(ids.get(c.i - 1));
                set_phrase(c, c.i - 1, e - l);
                if (x >= c.s) { ++c.step; return; }
            }
        }
        ++c.query;
        c.owner = this;
        u64 i = starts.pred_index(x);
        set_phrase(c, i, starts.access(i));
    }
    inline uint8_t at(u64 x) const { cursor &c = cur(); seek(c, x); return c.p[x - c.s]; }
    // the longest common prefix of S$[b..b+len) and x[0..len) (b + len <= n)
    inline u64 lce(u64 b, const uint8_t *x, u64 len) const {
        cursor &c = cur();
        u64 l = 0;
        while (l < len) {
            seek(c, b + l);
            u64 avail = std::min(len - l, c.e - (b + l));
            u64 q = lce_bytes(c.p + (b + l - c.s), x + l, avail);
            l += q;
            if (q < avail) break;
        }
        return l;
    }
    u64 bytes() const { return dict.size() + doff.bytes() + ids.bytes() + starts.bytes(); }

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
        dict.clear(); dict.reserve(d.size());
        std::vector<uint64_t> doffv(sigma + 1, 0);
        for (uint32_t r = 0; r < sigma; ++r) {
            uint32_t a = order[r];
            idof[a] = r + 1;
            dict.insert(dict.end(), d.begin() + off[a], d.begin() + off[a + 1]);
            doffv[r + 1] = dict.size();
        }
        std::vector<uint8_t>().swap(d);
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
        fprintf(stderr, "parse oracle: n=%lu, %lu phrases, %u distinct; dictionary %.1f MB, parse %.1f MB, starts %.1f MB (%.2f bits per character in all; %.1f s)\n",
                n, m, sigma, dict.size() / 1e6, ids.bytes() / 1e6, starts.bytes() / 1e6, 8.0 * bytes() / n, now() - t0);
    }
    void save(const std::string &f) const {
        FILE *o = open_or_die(f, "wb"); wmagic(o, "KTPTXT01");
        wr(o, n); wr(o, m); wr(o, K); wr(o, SM); wr(o, sigma);
        wvec(o, dict); doff.save(o); ids.save(o); starts.save(o);
        fclose(o);
    }
    void load(const std::string &f) {
        FILE *i = open_or_die(f, "rb"); rmagic(i, "KTPTXT01", f);
        rd(i, n); rd(i, m); rd(i, K); rd(i, SM); rd(i, sigma);
        rvec(i, dict); doff.load(i); ids.load(i); starts.load(i);
        fclose(i);
    }
};

}  // namespace katka
