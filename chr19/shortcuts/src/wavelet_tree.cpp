#include "wavelet_tree.hpp"

#include "io.hpp"

namespace rs {

std::unique_ptr<WaveletTree::Node> WaveletTree::build_node(const std::vector<uint16_t> &symbols, int depth,
                                                           int levels, BitvectorKind kind) {
  if (depth >= levels) return nullptr;
  auto node = std::make_unique<Node>();
  int shift = levels - 1 - depth;
  Bitvector::Builder bits;
  std::vector<uint16_t> left, right;
  for (uint16_t s : symbols) {
    bool bit = (s >> shift) & 1;
    bits.push_back(bit);
    (bit ? right : left).push_back(s);
  }
  node->bits = std::move(bits).finish(kind);
  node->left = build_node(left, depth + 1, levels, kind);
  node->right = build_node(right, depth + 1, levels, kind);
  return node;
}

WaveletTree WaveletTree::build(const std::vector<uint16_t> &symbols, int sigma, BitvectorKind kind) {
  WaveletTree wt;
  wt.sigma_ = sigma < 1 ? 1 : sigma;
  wt.size_ = static_cast<int64_t>(symbols.size());
  while ((1 << wt.levels_) < wt.sigma_) ++wt.levels_;
  wt.root_ = build_node(symbols, 0, wt.levels_, kind);
  return wt;
}

int64_t WaveletTree::rank(uint16_t symbol, int64_t i) const {
  const Node *node = root_.get();
  for (int depth = 0; depth < levels_; ++depth) {
    bool bit = (symbol >> (levels_ - 1 - depth)) & 1;
    int64_t ones = node->bits.rank1(i);
    i = bit ? ones : i - ones;
    node = bit ? node->right.get() : node->left.get();
  }
  return i;
}

int64_t WaveletTree::count_greater(uint16_t threshold, int64_t i) const {
  if (i <= 0) return 0;
  if (threshold >= sigma_ - 1) return 0;
  return i - count_at_most(root_.get(), 0, levels_, i, threshold, 0);
}

// Number of symbols <= threshold among the first i symbols of this node,
// whose subtree holds the symbols [smallest, smallest + 2^(levels-depth)).
// A child whose whole range is <= threshold is counted without descending.
int64_t WaveletTree::count_at_most(const Node *node, int depth, int levels, int64_t i, uint16_t threshold,
                                   uint16_t smallest) {
  if (i <= 0) return 0;
  if (depth >= levels) return smallest <= threshold ? i : 0;
  uint16_t half = static_cast<uint16_t>(1u << (levels - 1 - depth));
  int64_t right_i = node->bits.rank1(i), left_i = i - right_i;
  int64_t total = 0;
  for (int side = 0; side < 2; ++side) {
    int64_t child_i = side ? right_i : left_i;
    uint16_t lo = static_cast<uint16_t>(smallest + side * half), hi = static_cast<uint16_t>(lo + half - 1);
    if (child_i == 0 || lo > threshold) continue;
    if (hi <= threshold) {
      total += child_i;
    } else {
      total += count_at_most(side ? node->right.get() : node->left.get(), depth + 1, levels, child_i, threshold, lo);
    }
  }
  return total;
}

void WaveletTree::save_node(const Node *node, std::ostream &out) {
  node->bits.save(out);
  if (node->left) save_node(node->left.get(), out);
  if (node->right) save_node(node->right.get(), out);
}

std::unique_ptr<WaveletTree::Node> WaveletTree::load_node(std::istream &in, int depth, int levels) {
  if (depth >= levels) return nullptr;
  auto node = std::make_unique<Node>();
  node->bits = Bitvector::load(in);
  node->left = load_node(in, depth + 1, levels);
  node->right = load_node(in, depth + 1, levels);
  return node;
}

void WaveletTree::save(std::ostream &out) const {
  io::write(out, static_cast<int64_t>(sigma_));
  io::write(out, static_cast<int64_t>(levels_));
  io::write(out, size_);
  if (root_) save_node(root_.get(), out);
}

WaveletTree WaveletTree::load(std::istream &in) {
  WaveletTree wt;
  wt.sigma_ = static_cast<int>(io::read<int64_t>(in));
  wt.levels_ = static_cast<int>(io::read<int64_t>(in));
  wt.size_ = io::read<int64_t>(in);
  wt.root_ = load_node(in, 0, wt.levels_);
  return wt;
}

}  // namespace rs
