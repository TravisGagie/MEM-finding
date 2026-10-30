// Construction of the shortcuts (paper, Section 6).
#include "build.hpp"

#include <algorithm>
#include <stdexcept>

#include "text.hpp"

namespace rs {

namespace {

// An order-statistics structure over the rows of an entrance (a Fenwick
// tree over 0/1 values): delete a row, and find the k-th remaining row.
class RemainingRows {
 public:
  explicit RemainingRows(int64_t n) : n_(n), tree_(static_cast<size_t>(n) + 1, 0) {
    for (int64_t i = 1; i <= n; ++i) {
      tree_[static_cast<size_t>(i)] += 1;
      int64_t parent = i + (i & -i);
      if (parent <= n) tree_[static_cast<size_t>(parent)] += tree_[static_cast<size_t>(i)];
    }
    while (top_ * 2 <= n) top_ *= 2;
  }

  void remove(int64_t row) {
    for (int64_t i = row + 1; i <= n_; i += i & -i) --tree_[static_cast<size_t>(i)];
  }

  // The k-th remaining row, counting from 0.
  int64_t kth(int64_t k) const {
    int64_t pos = 0;
    ++k;
    for (int64_t step = top_; step > 0; step >>= 1) {
      if (pos + step <= n_ && tree_[static_cast<size_t>(pos + step)] < k) {
        pos += step;
        k -= tree_[static_cast<size_t>(pos)];
      }
    }
    return pos;
  }

 private:
  int64_t n_;
  int64_t top_ = 1;
  std::vector<int32_t> tree_;
};

struct Window {
  int64_t anchor, first_row, last_row;
};

bool all_acgt(const std::vector<uint8_t> &r, int64_t from, int64_t to) {
  for (int64_t k = from; k < to; ++k) {
    if (!is_acgt(r[static_cast<size_t>(k)])) return false;
  }
  return true;
}

// The coordinate in which anchors are multiples of ell and tiers are
// chosen: p itself, or, for a mirrored reference of length N (the reverse
// of the reference of another index), N - p, the boundary of the other
// index that p corresponds to.
int64_t boundary(int64_t p, int64_t mirror_length) { return mirror_length < 0 ? p : mirror_length - p; }

// Every anchor whose window is all A/C/G/T, with the window's interval.
// Anchors whose windows contain another character are dropped (the second
// pruning rule of Section 5), before collisions are discarded.
std::vector<Window> search_windows(const RLBWT &rlbwt, const std::vector<uint8_t> &R, int64_t ell,
                                   int64_t mirror_length, ConstructionStats &stats) {
  std::vector<Window> windows;
  // The first anchor: the smallest p >= ell whose boundary is a multiple of ell.
  const int64_t first = mirror_length < 0 ? ell : ell + ((mirror_length - ell) % ell + ell) % ell;
  for (int64_t p = first; p + ell <= static_cast<int64_t>(R.size()); p += ell) {
    ++stats.anchors;
    if (!all_acgt(R, p, p + ell)) {
      ++stats.dropped_unusable;  // the window is not all A/C/G/T
      continue;
    }
    int64_t lo = 0, hi = rlbwt.size() - 1;
    for (int64_t k = p + ell - 1; k >= p; --k) {
      rlbwt.backward_step(R[static_cast<size_t>(k)], lo, hi);
      ++stats.backward_steps;
    }
    if (hi < lo) throw std::runtime_error("a window of the reference does not occur in the text");
    windows.push_back({p, lo, hi});
  }
  return windows;
}

// Sorts the windows by interval and discards those whose window equals
// another's: equal windows have equal intervals, and different windows of
// the same length have disjoint intervals.
std::vector<Window> discard_collisions(std::vector<Window> windows, ConstructionStats &stats) {
  std::sort(windows.begin(), windows.end(), [](const Window &a, const Window &b) {
    return a.first_row != b.first_row ? a.first_row < b.first_row : a.anchor < b.anchor;
  });
  std::vector<Window> kept;
  for (size_t i = 0; i < windows.size();) {
    size_t j = i + 1;
    while (j < windows.size() && windows[j].first_row == windows[i].first_row) ++j;
    if (j == i + 1) {
      kept.push_back(windows[i]);
    } else {
      stats.discarded_colliding += static_cast<int64_t>(j - i);
    }
    i = j;
  }
  return kept;
}

// tau_p: the largest t <= max_tier with ell * 2^t dividing p's boundary
// and ell * 2^t <= p, lowered so that the tier-tau_p string is all
// A/C/G/T. -1 if even tier 0 is not.
int64_t top_tier_of(const std::vector<uint8_t> &R, int64_t p, int64_t ell, int64_t max_tier, int64_t mirror_length) {
  const int64_t b = boundary(p, mirror_length);
  int64_t tau = 0;
  while (tau < max_tier && b % (ell << (tau + 1)) == 0 && (ell << (tau + 1)) <= p) ++tau;
  int64_t clean = 0;  // number of A/C/G/T characters directly before p, up to ell * 2^tau
  while (clean < (ell << tau) && is_acgt(R[static_cast<size_t>(p - 1 - clean)])) ++clean;
  while (tau >= 0 && (ell << tau) > clean) --tau;
  return tau;
}

// Continues the backward search for the window of anchor p through
// R[p - ell * 2^tau .. p), tracking which rows of the entrance die at
// each step. Appends the rows' symbols to `symbols` and f_{p,t} for
// t = 0, ..., tau to `targets`.
void extend_anchor(const RLBWT &rlbwt, const std::vector<uint8_t> &R, const Window &w, int64_t ell, int64_t tau,
                   std::vector<uint8_t> &symbols, std::vector<int64_t> &targets, ConstructionStats &stats) {
  const int64_t p = w.anchor, rows = w.last_row - w.first_row + 1;
  const size_t sym0 = symbols.size();
  symbols.resize(sym0 + static_cast<size_t>(rows), static_cast<uint8_t>(tau + 1));  // alive to the end
  RemainingRows alive(rows);  // the set Q_k of the paper
  std::vector<std::pair<int64_t, int64_t>> dying;  // ranges of ranks in J_k

  int64_t lo = w.first_row, hi = w.last_row;  // J_k
  int64_t tiers_reached = 0;                  // |{t <= tau : ell * 2^t <= k}|
  for (int64_t k = 0; k < (ell << tau); ++k) {
    while (tiers_reached <= tau && (ell << tiers_reached) <= k) ++tiers_reached;
    const uint8_t c = R[static_cast<size_t>(p - 1 - k)];

    // The rows of J_k in runs of characters other than c die.
    dying.clear();
    for (int64_t run = rlbwt.run_containing(lo), row = lo; row <= hi; ++run) {
      int64_t end = std::min(hi, rlbwt.run_start(run + 1) - 1);
      if (rlbwt.run_char(run) != c) dying.emplace_back(row - lo, end - lo);
      row = end + 1;
    }
    // Delete from right to left, so the ranks still to come do not shift.
    for (auto it = dying.rbegin(); it != dying.rend(); ++it) {
      for (int64_t rank = it->second; rank >= it->first; --rank) {
        int64_t row = alive.kth(rank);
        symbols[sym0 + static_cast<size_t>(row)] = static_cast<uint8_t>(tiers_reached);
        alive.remove(row);
        ++stats.deaths;
      }
    }

    rlbwt.backward_step(c, lo, hi);
    ++stats.backward_steps;
    if (hi < lo) throw std::runtime_error("internal error: the reference does not occur in the text");

    // After k + 1 = ell * 2^t steps, J is I(alpha_{p,t} beta_p).
    int64_t done = k + 1;
    if (done % ell == 0 && ((done / ell) & (done / ell - 1)) == 0) targets.push_back(lo);
  }
}

}  // namespace

UnprunedShortcuts compute_shortcuts(const RLBWT &rlbwt, const std::vector<uint8_t> &reference, int64_t ell,
                                    int64_t max_tier, ConstructionStats &stats, bool mirrored) {
  if (ell <= 0 || max_tier < 0) throw std::invalid_argument("ell must be positive and max_tier non-negative");
  UnprunedShortcuts u;
  u.ell = ell;
  u.max_tier = max_tier;
  u.text_size = rlbwt.size();
  u.reference = reference;
  u.target_start.push_back(0);
  u.symbol_start.push_back(0);

  const int64_t mirror_length = mirrored ? static_cast<int64_t>(reference.size()) : -1;
  std::vector<Window> windows = discard_collisions(search_windows(rlbwt, reference, ell, mirror_length, stats), stats);
  for (const Window &w : windows) {
    int64_t tau = top_tier_of(reference, w.anchor, ell, max_tier, mirror_length);
    if (tau < 0) {
      ++stats.dropped_unusable;  // the tier-0 string is not all A/C/G/T
      continue;
    }
    extend_anchor(rlbwt, reference, w, ell, tau, u.symbols, u.target_first_row, stats);
    u.entrances.push_back({w.first_row, w.last_row, w.anchor});
    u.top_tier.push_back(static_cast<int8_t>(tau));
    u.target_start.push_back(static_cast<int64_t>(u.target_first_row.size()));
    u.symbol_start.push_back(static_cast<int64_t>(u.symbols.size()));
  }
  return u;
}

Shortcuts prune(const UnprunedShortcuts &u, int64_t theta, int64_t &dropped) {
  std::vector<Entrance> entrances;
  std::vector<int8_t> top_tier;
  std::vector<int64_t> targets;
  std::vector<uint16_t> symbols;
  dropped = 0;

  for (size_t e = 0; e < u.entrances.size(); ++e) {
    const uint8_t *sym = u.symbols.data() + u.symbol_start[e];
    const int64_t rows = u.symbol_start[e + 1] - u.symbol_start[e];
    const int64_t tau = u.top_tier[e];

    // survivors[t] = number of rows with symbol > t.
    std::vector<int64_t> survivors(static_cast<size_t>(tau) + 2, 0);
    for (int64_t i = 0; i < rows; ++i) {
      for (int64_t t = 0; t < sym[i]; ++t) ++survivors[static_cast<size_t>(t)];
    }
    int64_t top = -1;  // the largest tier kept; -1 drops the anchor
    while (top < tau && survivors[static_cast<size_t>(top + 1)] * (u.ell << (top + 1)) >= theta) ++top;
    if (top < 0) {
      ++dropped;
      continue;
    }

    entrances.push_back(u.entrances[e]);
    top_tier.push_back(static_cast<int8_t>(top));
    for (int64_t t = 0; t <= top; ++t) targets.push_back(u.target_first_row[u.target_start[e] + t]);
    for (int64_t i = 0; i < rows; ++i) symbols.push_back(std::min<uint16_t>(sym[i], static_cast<uint16_t>(top + 1)));
  }
  return Shortcuts::build(u.ell, u.max_tier, u.reference, entrances, top_tier, targets, symbols, u.text_size);
}

}  // namespace rs
