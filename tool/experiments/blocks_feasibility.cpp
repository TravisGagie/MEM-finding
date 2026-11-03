// katka-blocks: is it worth crossing matches with order-k character-level blocks?  A feasibility study,
// not a search method: it measures what blocks of continuations at the character level would save, on a
// text of both strands (e.g. SILVA), before anyone builds the index.
//
// The model.  For a k-mer kappa over A/C/G/T, succ(kappa) is the character (A/C/G/T) that most often
// follows kappa in S (ties: the smallest).  Position i of a string X is a continuation if X[i-k..i) is a
// k-mer over A/C/G/T and X[i] = succ(X[i-k..i)).  A block is a position that is not a continuation
// followed by its continuations: a block step would prepend a whole block, as the blocks of phrases of
// Section 5 do.  At the ends of a match the search uses T's BWT, as in the back-up lemma: crossing a
// match X[a..e) from its right end costs
//   - character steps from e down to the rightmost block start s >= a + k (a block start whose context
//     X[s-k..s) lies inside the match, so that its occurrences are block starts too),
//   - one block step per complete block between the leftmost and the rightmost such block starts,
//   - the leftmost complete block again in character steps (the back-up: the last block step may miss
//     occurrences whose first character is a continuation in the text), and
//   - character steps from the leftmost such block start s0 down to a.
// With fewer than two complete blocks this saves nothing, and the match costs e - a character steps, as now.
// So with block starts s0 < s1 < ... < sz in [a + k, e), z >= 2, crossing costs (e - a) - (sz - s1)
// character steps and z block steps, and the characters of [s1, sz) are crossed by block steps.
//
// Reported, for each k (default 8, 12, 16):
//   a. the number of distinct k-mers of S (over A/C/G/T);
//   b. the distribution of block lengths in S (and the fraction of characters in blocks of each length);
//   c. for the MEMs of a KATKA listing file (with the reads, to get their bases): steps to cross each MEM
//      under the model, against its length (the steps of plain backward search);
//   d. BML with threshold L (default 30) run on the same reads over the index PREFIX (RLBWT, or RLCSA if
//      PREFIX.csa exists): every backward step classified as in a failed window, at the ends of a
//      successful search (what the model still pays in character steps, including the failing step), or
//      inside a long crossing (what block steps would replace), with the block steps that would replace
//      them.
//
// usage: katka-blocks PREFIX reads listing [-k 8,12,16] [-L 30] [-n reads] [--check]
//   PREFIX.S (the text) and PREFIX.rlbwt (katka index) are needed; reads is FASTQ or FASTA (plain); the
//   listing is katka find's (or katka-classify's) output for those reads; -n uses the first n reads
//   (default 5000); --check also verifies the k-mer counts against a direct array for k <= 12.
// Memory: the text (1 byte per character), the index during the BML pass only (it is freed before the
// k-mer tables are built), and, for each k in turn, the counts: an array of 4^(k+1) 32-bit counts for
// k <= 13, else a hash table of 12-byte slots, at most 3/4 full (with 16-bit saturating counts), which
// briefly needs 1.5 times its size when it doubles.
#include "mems.hpp"
#include <fstream>
#include <map>

using namespace katka;
using std::string;
using std::vector;

// ------------------------------------------------------------------------------------------------
// counts of the character after each k-mer context: an array for small k, else a hash table
struct kmer_counts {
    int k = 0;
    bool direct = false;
    vector<uint32_t> arr;                   // direct: 4 counts per k-mer
    struct slot { uint32_t key; uint16_t c[4]; } __attribute__((packed));
    vector<slot> tab;                       // hash: a slot is used iff some count is non-zero
    u64 used = 0, mask = 0;

    void init(int kk, u64 expected, bool force_hash = false) {
        k = kk; direct = k <= 13 && !force_hash;
        if (direct) { arr.assign((u64)4 << (2 * k), 0); return; }
        u64 sz = 1;
        while (sz < expected * 4 / 3 + 16) sz <<= 1;
        tab.assign(sz, slot{0, {0, 0, 0, 0}}); mask = sz - 1; used = 0;
    }
    static inline u64 hash(u64 x) { return mix64(x); }
    inline slot *find(uint32_t key, bool insert) {
        for (u64 h = hash(key) & mask;; h = (h + 1) & mask) {
            slot &s = tab[h];
            bool empty = !(s.c[0] | s.c[1] | s.c[2] | s.c[3]);
            if (empty) {
                if (!insert) return nullptr;
                s.key = key; ++used;
                return &s;
            }
            if (s.key == key) return &s;
        }
    }
    void grow() {
        vector<slot> old; old.swap(tab);
        tab.assign(old.size() * 2, slot{0, {0, 0, 0, 0}}); mask = tab.size() - 1; used = 0;
        for (auto &s : old) if (s.c[0] | s.c[1] | s.c[2] | s.c[3]) { slot *t = find(s.key, true); memcpy(t->c, s.c, sizeof s.c); }
    }
    inline void add(u64 kmer, int next) {
        if (direct) { arr[kmer * 4 + next]++; return; }
        if (used * 4 >= tab.size() * 3) grow();
        slot *s = find((uint32_t)kmer, true);
        if (s->c[next] < 65535) s->c[next]++;            // saturating: only the argmax matters
    }
    // succ: the most frequent next character (0..3), or -1 if the k-mer never has one
    inline int succ(u64 kmer) {
        const uint32_t *c32 = nullptr; uint32_t c[4];
        if (direct) c32 = &arr[kmer * 4];
        else { slot *s = find((uint32_t)kmer, false); if (!s) return -1; for (int i = 0; i < 4; ++i) c[i] = s->c[i]; c32 = c; }
        int best = -1; uint32_t bc = 0;
        for (int i = 0; i < 4; ++i) if (c32[i] > bc) { bc = c32[i]; best = i; }
        return best;
    }
    u64 distinct() const {
        if (!direct) return used;
        u64 d = 0;
        for (u64 i = 0; i < arr.size(); i += 4) d += (arr[i] | arr[i + 1] | arr[i + 2] | arr[i + 3]) != 0;
        return d;
    }
    u64 bytes() const { return direct ? arr.size() * 4 : tab.size() * sizeof(slot); }
};

// calls f(i, kmer code of s[i-k..i)) for every i such that s[i-k..i) is over A/C/G/T (i in [k, n])
template <class F> void for_each_context(const uint8_t *s, u64 n, int k, F f) {
    u64 code = 0, len = 0, mask = k == 32 ? ~0ull : (1ull << (2 * k)) - 1;
    for (u64 i = 0; i < n; ++i) {
        int c = base_code(s[i]);
        if (len >= (u64)k) f(i, code);
        if (c < 0) { len = 0; code = 0; continue; }
        code = ((code << 2) | (u64)c) & mask; ++len;
    }
    if (len >= (u64)k) f(n, code);
}

// the block starts of X[0..n): positions i that are not continuations (continuation: X[i-k..i) over
// A/C/G/T and X[i] = succ of it)
static void block_starts(kmer_counts &K, const char *X, u64 n, vector<char> &start) {
    start.assign(n, 1);
    for_each_context((const uint8_t *)X, n, K.k, [&](u64 i, u64 code) {
        if (i < n && base_code(X[i]) >= 0 && K.succ(code) == base_code(X[i])) start[i] = 0;
    });
}

// the model's cost of crossing X[a..e) from its right end (see the top of the file), given the block starts
// of X; also the block steps it takes and the characters those steps cover
struct cross_cost { u64 char_steps = 0, block_steps = 0, covered = 0; };   // covered: characters crossed by block steps (not redone)
static cross_cost crossing(const vector<char> &start, u64 a, u64 e, int k) {
    cross_cost c;
    vector<u64> B;
    for (u64 i = a + k; i < e; ++i) if (start[i]) B.push_back(i);
    if (B.size() < 3) { c.char_steps = e - a; return c; }       // fewer than two complete blocks: nothing saved
    c.block_steps = B.size() - 1;                                // one block step per complete block (the last one backed up)
    c.covered = B.back() - B[1];
    c.char_steps = (e - a) - c.covered;                          // right end, the backed-up block, the left end
    return c;
}

// ------------------------------------------------------------------------------------------------
struct reader {
    std::ifstream in; bool fasta = false; string pending;
    explicit reader(const string &f) : in(f) { if (!in) die("cannot open " + f); fasta = in.peek() == '>'; }
    bool next(string &name, string &seq) {
        string h, line;
        if (!fasta) { string p, q; if (!(std::getline(in, h) && std::getline(in, seq) && std::getline(in, p) && std::getline(in, q))) return false; }
        else {
            if (pending.empty()) { while (std::getline(in, line)) if (!line.empty() && line[0] == '>') { pending = line; break; } if (pending.empty()) return false; }
            h.swap(pending); pending.clear(); seq.clear();
            while (std::getline(in, line)) { if (!line.empty() && line[0] == '>') { pending = line; break; } if (!line.empty() && line.back() == '\r') line.pop_back(); seq += line; }
        }
        for (auto &c : seq) c = (char)toupper((unsigned char)c);
        size_t b = h.find_first_of(" \t"); name = h.substr(1, b == string::npos ? string::npos : b - 1);
        return true;
    }
};

// the MEM spans of a listing file, record by record in file order (= the order of the reads): (name,
// [a, b] inclusive).  Records are paired with reads by position, not by name, since names can repeat
// (e.g. reads of several regions or mates in one file); main checks that the names agree.
static vector<std::pair<string, vector<std::pair<u64, u64>>>> read_listing(const string &f) {
    vector<std::pair<string, vector<std::pair<u64, u64>>>> m;
    std::ifstream in(f);
    if (!in) die("cannot open " + f);
    string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') { m.push_back({line.substr(1), {}}); continue; }
        if (m.empty()) continue;
        for (size_t p = line.find('['); p != string::npos; p = line.find('[', p + 1)) {
            u64 a, b;
            if (sscanf(line.c_str() + p, "[%lu,%lu]", &a, &b) == 2) m.back().second.push_back({a, b});
        }
    }
    return m;
}

static void hist_line(const char *what, const vector<u64> &h, const vector<u64> &hc, u64 total, u64 totc) {
    printf("    %s\n", what);
    const char *lab[] = {"1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256", "257-512", "513-1024", ">1024"};
    for (size_t b = 0; b < h.size(); ++b)
        if (h[b]) printf("      %-9s %12lu blocks (%6.2f%%)  %14lu characters (%6.2f%%)\n", lab[b], h[b], 100.0 * h[b] / total, hc[b], 100.0 * hc[b] / totc);
}
static size_t bucket(u64 len) { size_t b = 0; u64 x = 1; while (b < 11 && len > x) { x <<= 1; ++b; } return b; }

int main(int argc, char **argv) {
    if (argc < 4) die("usage: katka-blocks PREFIX reads listing [-k 8,12,16] [-L 30] [-n reads] [--check]");
    const string prefix = argv[1], readfile = argv[2], listfile = argv[3];
    vector<int> ks{8, 12, 16};
    u64 L = 30, maxreads = 5000;
    bool check = false;
    for (int i = 4; i < argc; ++i) {
        string a = argv[i];
        if (a == "-k" && i + 1 < argc) { ks.clear(); string v = argv[++i], w; for (char c : v + ",") { if (c == ',') { if (!w.empty()) ks.push_back(std::stoi(w)); w.clear(); } else w += c; } }
        else if (a == "-L" && i + 1 < argc) L = std::stoull(argv[++i]);
        else if (a == "-n" && i + 1 < argc) maxreads = std::stoull(argv[++i]);
        else if (a == "--check") check = true;
        else die("unknown option " + a);
    }
    for (int k : ks) if (k < 1 || k > 16) die("k must be between 1 and 16");
    double t0 = now();
    // the text
    const u64 n = file_size(prefix + ".S");
    vector<uint8_t> S(n);
    { FILE *f = open_or_die(prefix + ".S", "rb"); if (fread(S.data(), 1, n, f) != n) die("cannot read the text"); fclose(f); }
    // the reads and the listing
    vector<string> names, seqs;
    { reader R(readfile); string nm, sq; while (names.size() < maxreads && R.next(nm, sq)) { names.push_back(nm); seqs.push_back(sq); } }
    auto listing = read_listing(listfile);
    if (listing.size() < names.size()) die("the listing has " + std::to_string(listing.size()) + " records for " + std::to_string(names.size()) + " reads");
    for (u64 r = 0; r < names.size(); ++r)
        if (listing[r].first != names[r]) die("listing record " + std::to_string(r) + " is " + listing[r].first + ", read " + std::to_string(r) + " is " + names[r] + ": the listing must be of these reads, in order");
    // BML with threshold L on the reads, recording each search's stretch (needed for d)
    rlbwt B;
    { FILE *f = open_or_die(prefix + ".rlbwt", "rb"); rmagic(f, "KTRLBW01", prefix + ".rlbwt"); B.load(f); fclose(f); }
    rlcsa X;
    stepper st; st.bwt = &B;
    if (file_exists(prefix + ".csa")) { X.load(prefix + ".csa"); if (X.n == B.n) st.csa = &X; }
    struct stretch { u64 read, lo, hi; bool rc, ok, stopped; };   // a backward search matching X[lo..hi) (X = read or rc); ok: not a failed
                                                                   // window; stopped: it ended with a step that failed
    vector<stretch> stretches;
    vector<string> rcs(seqs.size());
    u64 total_steps = 0, nmems_bml = 0, nmems_list = 0, mem_mismatch = 0;
    for (u64 r = 0; r < seqs.size(); ++r) {
        rcs[r] = reverse_complement(seqs[r]);
        vector<std::pair<u64, u64>> found;
        // BML exactly as mems.hpp's bml_segment, with each search's steps attributed
        const string &P = seqs[r];
        for (u64 i = 0; i < P.size();) {
            while (i < P.size() && base_code(P[i]) < 0) ++i;
            u64 r0 = i;
            while (i < P.size() && base_code(P[i]) >= 0) ++i;
            u64 m = i - r0;
            if (m < L) continue;
            const char *R = P.data() + r0;
            u64 x = 0;
            while (x + L <= m) {
                interval I = st.full();
                u64 j = x + L;
                while (j > x) { ++total_steps; if (!st.step(I, R[j - 1])) break; --j; }
                if (j > x) { stretches.push_back({r, r0 + j, r0 + x + L, false, false, true}); x = j; continue; }   // failed window
                u64 a = j;                                                // extend left (the same search)
                while (a > 0) { ++total_steps; if (!st.step(I, R[a - 1])) break; --a; }
                stretches.push_back({r, r0 + a, r0 + x + L, false, true, a > 0});
                interval J = st.full();                                   // forward: backward on the rc
                u64 e = a;
                while (e < m) { ++total_steps; if (!st.step(J, complement(R[e]))) break; ++e; }
                        stretches.push_back({r, P.size() - (r0 + e), P.size() - (r0 + a), true, true, e < m});
                found.push_back({r0 + a, r0 + e - 1});
                x = e + 1 >= L ? e + 1 - L : 0;
            }
        }
        nmems_bml += found.size();
        nmems_list += listing[r].second.size();
        if (listing[r].second != found) ++mem_mismatch;
    }
    B = rlbwt(); X = rlcsa();                       // the index is not needed any more: free it for the k-mer tables
    printf("text: %lu characters; reads: %zu; BML with L = %lu: %lu MEMs, %lu backward steps (%.1f s)\n", n, seqs.size(), L, nmems_bml, total_steps, now() - t0);
    printf("listing %s: %lu MEMs for these reads; %lu reads whose MEMs differ from BML's (L = %lu)%s\n", listfile.c_str(), nmems_list, mem_mismatch, L,
           mem_mismatch ? " -- the listing was made with other options; (c) uses the listing, (d) BML's own MEMs" : "");

    for (int k : ks) {
        double t1 = now();
        printf("\n== k = %d\n", k);
        kmer_counts K;
        K.init(k, std::min<u64>(n, (u64)1 << 28) / 4);
        for_each_context(S.data(), n, k, [&](u64 i, u64 code) { if (i < n) { int c = base_code(S[i]); if (c >= 0) K.add(code, c); } });
        // a. distinct k-mers (those followed by a base)
        printf("  a. distinct k-mers followed by a base: %lu  (counts table %.2f GB)\n", K.distinct(), K.bytes() / 1e9);
        if (check && k <= 12) {          // the hash table against direct counts
            kmer_counts D; D.init(k, 1 << 16, true);
            for_each_context(S.data(), n, k, [&](u64 i, u64 code) { if (i < n) { int c = base_code(S[i]); if (c >= 0) D.add(code, c); } });
            u64 bad = 0;
            for (u64 x = 0; x < ((u64)1 << (2 * k)); ++x) if (K.succ(x) != D.succ(x)) ++bad;
            printf("     check: hash table and array agree on succ for all k-mers: %s (%lu differences)\n", bad ? "NO" : "yes", bad);
        }
        // b. block lengths in S
        {
            vector<u64> h(12, 0), hc(12, 0);
            u64 nb = 0, run = 0, chars = 0;
            auto close = [&]() { if (run) { size_t b = bucket(run); h[b]++; hc[b] += run; ++nb; chars += run; } };
            // one pass with a rolling code: S[i] continues the block if S[i-k..i) is over A/C/G/T and S[i] = succ of it
            u64 code = 0, len = 0, mask = (1ull << (2 * k)) - 1;
            for (u64 i = 0; i < n; ++i) {
                int c = base_code(S[i]);
                if (c < 0) { close(); run = 0; len = 0; code = 0; continue; }   // separators are not in blocks
                bool cont = len >= (u64)k && K.succ(code) == c;
                if (!cont) { close(); run = 1; } else ++run;
                code = ((code << 2) | (u64)c) & mask; ++len;
            }
            close();
            printf("  b. blocks in S: %lu over %lu bases, %.2f bases per block on average\n", nb, chars, (double)chars / std::max<u64>(1, nb));
            hist_line("block length       blocks                       bases in them", h, hc, nb, chars);
        }
        // c. crossing the listing's MEMs
        {
            u64 nm = 0, len = 0, steps = 0, bsteps = 0;
            vector<u64> ratio_h(5, 0);      // cost / length: <0.1, <0.25, <0.5, <0.75, >=0.75
            vector<u64> lenb(12, 0), stepb(12, 0), cntb(12, 0);
            vector<char> start;
            for (u64 r = 0; r < seqs.size(); ++r) {
                block_starts(K, seqs[r].data(), seqs[r].size(), start);
                for (auto [a, b] : listing[r].second) {
                    u64 e = b + 1;
                    if (e > seqs[r].size()) continue;
                    cross_cost c = crossing(start, a, e, k);
                    u64 cost = c.char_steps + c.block_steps;
                    ++nm; len += e - a; steps += cost; bsteps += c.block_steps;
                    double q = (double)cost / (e - a);
                    ratio_h[q < 0.1 ? 0 : q < 0.25 ? 1 : q < 0.5 ? 2 : q < 0.75 ? 3 : 4]++;
                    size_t lb = bucket(e - a); lenb[lb] += e - a; stepb[lb] += cost; cntb[lb]++;
                }
            }
            printf("  c. the listing's MEMs: %lu, %lu characters; crossing them takes %lu steps under the model (%lu of them block steps): %.3f steps per character (plain backward search: 1)\n",
                   nm, len, steps, bsteps, (double)steps / std::max<u64>(1, len));
            printf("     steps / length per MEM:  <0.1: %lu   0.1-0.25: %lu   0.25-0.5: %lu   0.5-0.75: %lu   >=0.75: %lu\n", ratio_h[0], ratio_h[1], ratio_h[2], ratio_h[3], ratio_h[4]);
            const char *lab[] = {"1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256", "257-512", "513-1024", ">1024"};
            for (size_t b = 0; b < 12; ++b) if (cntb[b]) printf("     MEM length %-9s %8lu MEMs, %.3f steps per character\n", lab[b], cntb[b], (double)stepb[b] / lenb[b]);
        }
        // d. BML's backward steps
        {
            u64 failed = 0, ends = 0, inside = 0, blocks = 0;
            vector<char> sf, sr;
            u64 cur = UINT64_MAX;
            for (auto &s : stretches) {
                if (s.read != cur) { cur = s.read; block_starts(K, seqs[cur].data(), seqs[cur].size(), sf); block_starts(K, rcs[cur].data(), rcs[cur].size(), sr); }
                u64 len = s.hi - s.lo;
                if (!s.ok) { failed += len + 1; continue; }               // the window's steps, and the step that failed
                cross_cost c = crossing(s.rc ? sr : sf, s.lo, s.hi, k);
                ends += c.char_steps + s.stopped; inside += c.covered; blocks += c.block_steps;
                if (c.char_steps + c.covered != len) die("internal error: crossing accounting");
            }
            u64 tot = failed + ends + inside;
            if (tot != total_steps) die("internal error: " + std::to_string(tot) + " steps classified, " + std::to_string(total_steps) + " taken");
            printf("  d. BML (L = %lu): its %lu backward steps, by where they fall:\n", L, tot);
            printf("     in failed windows:          %12lu (%5.1f%%)\n", failed, 100.0 * failed / tot);
            printf("     at the ends of searches:    %12lu (%5.1f%%)  (still character steps under the model)\n", ends, 100.0 * ends / tot);
            printf("     inside long crossings:      %12lu (%5.1f%%)  -> %lu block steps instead\n", inside, 100.0 * inside / tot, blocks);
            printf("     total under the model: %lu steps instead of %lu (%.2f times fewer)\n", failed + ends + blocks, tot, (double)tot / std::max<u64>(1, failed + ends + blocks));
        }
        printf("  (k = %d: %.1f s)\n", k, now() - t1);
    }
    return 0;
}
