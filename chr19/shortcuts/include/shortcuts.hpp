// The shortcuts for a linear reference R (paper, Sections 4 and 5).
//
// Notation, as in the paper:
//   ell          the window length;
//   max_tier     tau, so shortcuts cross ell * 2^t characters, 0 <= t <= tau;
//   anchor p     a multiple of ell; its window is beta_p = R[p .. p + ell)
//                and its tier-t string is alpha_{p,t} = R[p - ell*2^t .. p);
//   entrance     the BWT interval I(beta_p) of an anchor's window;
//   top tier     tau_p, the largest tier anchor p offers, after it is
//                lowered for characters other than A/C/G/T and by pruning;
//   sym(i)       for a row i of an entrance, the number of tiers whose
//                strings precede row i's suffix in T; row i survives
//                tier t if sym(i) > t;
//   f_{p,t}      the first row of I(alpha_{p,t} beta_p).
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "bitvector.hpp"
#include "int_vector.hpp"
#include "wavelet_tree.hpp"

namespace rs {

struct Entrance {
  int64_t first_row, last_row;  // I(beta_p)
  int64_t anchor;               // p
};

// The data structure of Section 5. Entrance e is the e-th kept anchor in
// the order of the entrances' first rows; entrances are disjoint.
class Shortcuts {
 public:
  int64_t ell = 0;
  int64_t max_tier = 0;

  // `targets` holds f_{p,t} for t = 0, ..., top_tier[e] for each entrance
  // in turn, and `symbols` sym(i) for every row of every entrance, in row
  // order. text_size is n, the number of rows of the BWT.
  static Shortcuts build(int64_t ell, int64_t max_tier, const std::vector<uint8_t> &reference,
                         const std::vector<Entrance> &entrances, const std::vector<int8_t> &top_tier,
                         const std::vector<int64_t> &targets, const std::vector<uint16_t> &symbols,
                         int64_t text_size);

  int64_t num_entrances() const { return anchors_.size(); }
  int64_t num_rows() const { return symbols_.size(); }  // W, the rows of all entrances
  int64_t first_row(int64_t e) const { return entrance_starts_.select1(e); }
  int64_t rows(int64_t e) const { return symbol_start_[e + 1] - symbol_start_[e]; }
  int64_t last_row(int64_t e) const { return first_row(e) + rows(e) - 1; }
  int64_t anchor(int64_t e) const { return anchors_[e]; }
  int64_t top_tier(int64_t e) const { return top_tier_[e]; }  // tau_p, after capping and pruning
  int64_t target(int64_t e, int64_t t) const { return targets_[target_start_[e] + t]; }  // f_{p,t}
  uint8_t reference_char(int64_t i) const { return reference_[i]; }  // R[i]

  // The entrance containing [lo, hi], or -1: one rank query on
  // entrance_starts_ gives the last entrance starting at or before lo, the
  // only one that can contain [lo, hi].
  // If `first` is not null, the entrance's first row is stored there.
  int64_t find_entrance(int64_t lo, int64_t hi, int64_t *first = nullptr) const {
    int64_t f = 0;
    const int64_t e = entrance_starts_.rank1_and_predecessor(lo + 1, f) - 1;
    if (e < 0 || hi > f + rows(e) - 1) return -1;
    if (first) *first = f;
    return e;
  }

  // Lemma 3: given the interval [lo, hi] of X, contained in entrance e
  // whose first row is `first`, computes the interval of alpha_{p,t} X.
  // The marked rows are those that survive tier t. Returns false if the
  // new interval is empty.
  bool jump(int64_t e, int64_t first, int64_t t, int64_t lo, int64_t hi, int64_t &new_lo, int64_t &new_hi) const;
  bool jump(int64_t e, int64_t t, int64_t lo, int64_t hi, int64_t &new_lo, int64_t &new_hi) const {
    return jump(e, first_row(e), t, lo, hi, new_lo, new_hi);
  }

  void save(const std::string &path) const;
  static Shortcuts load(const std::string &path);

 private:
  PackedDna reference_;        // R, 2 bits per character
  Bitvector entrance_starts_;  // sparse, over the n rows: the first row of each entrance
  IntVector anchors_;          // p
  IntVector top_tier_;
  IntVector symbol_start_;     // E + 1 entries: where each entrance's rows start in symbols_
  IntVector target_start_;     // E + 1 entries: where each entrance's f_{p,t} start in targets_
  IntVector targets_;
  WaveletTree symbols_;        // sym(i) for every row of every entrance
};

}  // namespace rs
