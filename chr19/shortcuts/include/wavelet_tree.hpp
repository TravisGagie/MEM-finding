// A balanced wavelet tree over a small alphabet {0, ..., sigma-1}. The
// RLBWT uses one over its run characters; the shortcuts use one over the
// symbols sym(i) of the rows of the entrances (paper, Section 5).
#pragma once
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <vector>

#include "bitvector.hpp"

namespace rs {

class WaveletTree {
 public:
  WaveletTree() = default;
  WaveletTree(WaveletTree &&) = default;
  WaveletTree &operator=(WaveletTree &&) = default;

  // Every symbol must be in [0, sigma).
  static WaveletTree build(const std::vector<uint16_t> &symbols, int sigma, BitvectorKind kind);

  int64_t size() const { return size_; }

  // Number of occurrences of `symbol` among the first i symbols.
  int64_t rank(uint16_t symbol, int64_t i) const;

  // Number of symbols greater than `threshold` among the first i symbols,
  // in one root-to-leaves pass.
  int64_t count_greater(uint16_t threshold, int64_t i) const;

  void save(std::ostream &out) const;
  static WaveletTree load(std::istream &in);

 private:
  struct Node {
    Bitvector bits;  // bit 1 sends a symbol to the right child
    std::unique_ptr<Node> left, right;
  };

  int sigma_ = 1;
  int levels_ = 0;  // ceil(log2 sigma)
  int64_t size_ = 0;
  std::unique_ptr<Node> root_;

  static std::unique_ptr<Node> build_node(const std::vector<uint16_t> &symbols, int depth, int levels,
                                          BitvectorKind kind);
  static int64_t count_at_most(const Node *node, int depth, int levels, int64_t i, uint16_t threshold,
                               uint16_t smallest);
  static void save_node(const Node *node, std::ostream &out);
  static std::unique_ptr<Node> load_node(std::istream &in, int depth, int levels);
};

}  // namespace rs
