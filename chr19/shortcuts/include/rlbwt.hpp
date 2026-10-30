// The run-length compressed BWT (paper, Section 2), supporting backward
// steps and listing the runs that meet an interval.
//
// BWT[i] = T[SA[i] - 1], with BWT[i] = $ for the row with SA[i] = 0. The
// r runs are stored as
//   - run_starts_: a bitvector over the n rows marking the first row of
//     each run;
//   - run_char_: the character of each run, and run_char_wt_, a wavelet
//     tree over the same sequence;
//   - sorted_run_starts_: a bitvector marking where each run starts once
//     the runs are stably sorted by character, i.e., where its
//     characters are in the first column F.
// rank_c(BWT, i) is then: find the run containing row i-1, count the runs
// of c before it with the wavelet tree, and look up where that many runs
// of c end in F.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "bitvector.hpp"
#include "wavelet_tree.hpp"

namespace rs {

class RLBWT {
 public:
  RLBWT() = default;

  // From T, computing its suffix array with libsais (8n bytes of memory).
  static RLBWT from_text(const std::vector<uint8_t> &text);

  // From a BWT file with one byte per row, such as Big-BWT writes for a
  // collection file with one sequence per line: `separator` (normally
  // '\n') plays the role of #, and the end-of-string byte, the only byte
  // smaller than 'A' other than `separator` that occurs exactly once,
  // plays the role of $. Reads the file twice and never holds the BWT.
  static RLBWT from_bwt_file(const std::string &path, uint8_t separator = '\n');

  int64_t size() const { return n_; }
  int64_t num_runs() const { return static_cast<int64_t>(run_char_.size()); }

  // Number of occurrences of character c in BWT[0..i).
  int64_t rank(uint8_t c, int64_t i) const;

  // One backward step: from the interval [lo, hi] of X to that of cX.
  // Empty intervals have hi < lo.
  void backward_step(uint8_t c, int64_t &lo, int64_t &hi) const;

  // The runs, for the construction of the shortcuts.
  int64_t run_containing(int64_t row) const { return run_starts_.rank1(row + 1) - 1; }
  int64_t run_start(int64_t run) const { return run < num_runs() ? run_starts_.select1(run) : n_; }
  uint8_t run_char(int64_t run) const { return alphabet_[run_char_[static_cast<size_t>(run)]]; }

  void save(const std::string &path) const;
  static RLBWT load(const std::string &path);

 private:
  int64_t n_ = 0;
  std::vector<uint8_t> alphabet_;          // the distinct characters, sorted
  int16_t code_[256];                      // character -> index in alphabet_, or -1
  std::vector<int64_t> C_;                 // C_[a] = number of characters smaller than alphabet_[a]
  std::vector<int64_t> runs_before_;       // runs_before_[a] = number of runs of characters smaller than alphabet_[a]
  std::vector<uint16_t> run_char_;         // index in alphabet_ of each run's character
  Bitvector run_starts_;
  WaveletTree run_char_wt_;
  Bitvector sorted_run_starts_;

  class Builder;
  void set_codes();
};

}  // namespace rs
