// A range-minimum structure that does not keep the array (for Muthukrishnan-Sadakane listing, Section 4
// of the paper, "Counting"): about 4 bits per element plus o(n).
//
// For A[0..n), let node i's parent be the last j < i with A[j] < A[i] (a virtual root, -infinity, if
// none), with children in index order, and store the balanced parentheses (BP) of that tree, the root
// first; the open parenthesis of node i is the (i+2)-th '(' and the excess after it is depth(i) + 1.
// For l < r, let m be the rightmost minimum of A[l..r] and p its parent.  Every node in [l, r] is in p's
// subtree, below p; the rightmost minimum of the excess over [open(l), open(r)] is depth(m) at the closing
// parenthesis just before open(m), unless m = l, when it is depth(l) + 1, its value at open(l).  So
//   rmq(l, r) = l if that minimum equals the excess at open(l), else the node opened just after it.
// Ties among minima are resolved to the rightmost; listing does not care which minimum it gets, as long as
// it recurses into the left part before the right part (tag_array::list).
//
// The excess is stored per 64-bit word (its value before the word, 32 bits) and a segment tree over the
// minima of the words finds the rightmost minimum among whole words.
#pragma once
#include "bits.hpp"

namespace katka {

struct succinct_rmq {
    uint64_t n = 0;
    BitRank bp;                         // 2n + 2 parentheses, 1 = '('
    std::vector<int32_t> before;        // excess before each word
    std::vector<int32_t> tree;          // segment tree over word minima (absolute), size 2P
    uint64_t P = 1, nw = 0;

    // calls next() n times to get A[0], A[1], ...
    template <class G> void build(uint64_t size, G next) {
        n = size;
        if (2 * n + 2 >= (1ull << 31)) throw std::runtime_error("rmq: too many elements");
        bp.init(2 * n + 2);
        uint64_t pos = 0;
        std::vector<uint64_t> stk;      // values of the open ancestors
        bp.set(pos++);                  // the root
        for (uint64_t i = 0; i < n; ++i) {
            uint64_t a = next();
            while (!stk.empty() && stk.back() >= a) { stk.pop_back(); ++pos; }
            bp.set(pos++); stk.push_back(a);
        }
        pos += stk.size() + 1;
        if (pos != 2 * n + 2) throw std::runtime_error("rmq: internal error");
        bp.finish();
        nw = (2 * n + 2 + 63) / 64;
        before.assign(nw + 1, 0);
        while (P < nw) P <<= 1;
        tree.assign(2 * P, INT32_MAX);
        int32_t e = 0;
        for (uint64_t w = 0; w < nw; ++w) {
            before[w] = e;
            int32_t mn = INT32_MAX;
            for (uint64_t b = 0; b < 64 && w * 64 + b < 2 * n + 2; ++b) {
                e += bp.get(w * 64 + b) ? 1 : -1;
                mn = std::min(mn, e);
            }
            tree[P + w] = mn;
        }
        before[nw] = e;
        for (uint64_t k = P - 1; k >= 1; --k) tree[k] = std::min(tree[2 * k], tree[2 * k + 1]);
    }

    // rightmost minimum of E over positions [x, y] of one word: (value, position)
    std::pair<int32_t, uint64_t> word_min(uint64_t x, uint64_t y) const {
        uint64_t w = x >> 6, word = bp.bits[w];
        int32_t e = before[w];
        for (uint64_t b = 0; b < (x & 63); ++b) e += ((word >> b) & 1) ? 1 : -1;
        int32_t mn = INT32_MAX; uint64_t at = x;
        for (uint64_t p = x; p <= y; ++p) {
            e += ((word >> (p & 63)) & 1) ? 1 : -1;
            if (e <= mn) { mn = e; at = p; }
        }
        return {mn, at};
    }
    int32_t tree_min(uint64_t a, uint64_t b) const {     // min over words [a, b]
        int32_t mn = INT32_MAX;
        for (a += P, b += P + 1; a < b; a >>= 1, b >>= 1) {
            if (a & 1) mn = std::min(mn, tree[a++]);
            if (b & 1) mn = std::min(mn, tree[--b]);
        }
        return mn;
    }
    // rightmost word in [a, b] whose minimum is v (one exists)
    uint64_t rightmost_word(uint64_t node, uint64_t nl, uint64_t nr, uint64_t a, uint64_t b, int32_t v) const {
        if (nr < a || nl > b || tree[node] > v) return UINT64_MAX;
        if (nl == nr) return nl;
        uint64_t mid = (nl + nr) >> 1;
        uint64_t r = rightmost_word(2 * node + 1, mid + 1, nr, a, b, v);
        return r != UINT64_MAX ? r : rightmost_word(2 * node, nl, mid, a, b, v);
    }
    // rightmost minimum of E over [x, y]: (value, position)
    std::pair<int32_t, uint64_t> range_min(uint64_t x, uint64_t y) const {
        uint64_t wx = x >> 6, wy = y >> 6;
        if (wx == wy) return word_min(x, y);
        auto L = word_min(x, wx * 64 + 63), R = word_min(wy * 64, y);
        int32_t mid = wy > wx + 1 ? tree_min(wx + 1, wy - 1) : INT32_MAX;
        int32_t mn = std::min(L.first, std::min(mid, R.first));
        if (R.first == mn) return R;
        if (mid == mn) {
            uint64_t w = rightmost_word(1, 0, P - 1, wx + 1, wy - 1, mn);
            return word_min(w * 64, w * 64 + 63);
        }
        return L;
    }
    // position of a minimum of A[l..r]
    uint64_t operator()(uint64_t l, uint64_t r) const {
        if (l >= r) return l;
        uint64_t ol = bp.select1(l + 1), orr = bp.select1(r + 1);
        auto mn = range_min(ol, orr);
        int32_t el = before[ol >> 6];
        for (uint64_t b = 0; b <= (ol & 63); ++b) el += ((bp.bits[ol >> 6] >> b) & 1) ? 1 : -1;
        if (mn.first == el) return l;
        return bp.rank1(mn.second + 2) - 2;
    }
    uint64_t bytes() const { return bp.bytes() + before.size() * 4 + tree.size() * 4 + 24; }
    void save(FILE *f) const { wr(f, n); wr(f, P); wr(f, nw); bp.save(f); wvec(f, before); wvec(f, tree); }
    void load(FILE *f) { rd(f, n); rd(f, P); rd(f, nw); bp.load(f); rvec(f, before); rvec(f, tree); }
};

}  // namespace katka
