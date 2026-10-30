#include "shortcuts.hpp"

#include <fstream>
#include <stdexcept>

#include "io.hpp"

namespace rs {

Shortcuts Shortcuts::build(int64_t ell, int64_t max_tier, const std::vector<uint8_t> &reference,
                           const std::vector<Entrance> &entrances, const std::vector<int8_t> &top_tier,
                           const std::vector<int64_t> &targets, const std::vector<uint16_t> &symbols,
                           int64_t text_size) {
  Shortcuts s;
  s.ell = ell;
  s.max_tier = max_tier;
  s.reference_ = PackedDna(reference);
  std::vector<int64_t> first_rows, anchors, tiers, symbol_start = {0}, target_start = {0};
  for (size_t e = 0; e < entrances.size(); ++e) {
    first_rows.push_back(entrances[e].first_row);
    anchors.push_back(entrances[e].anchor);
    tiers.push_back(top_tier[e]);
    symbol_start.push_back(symbol_start.back() + entrances[e].last_row - entrances[e].first_row + 1);
    target_start.push_back(target_start.back() + top_tier[e] + 1);
  }
  if (symbol_start.back() != static_cast<int64_t>(symbols.size()) ||
      target_start.back() != static_cast<int64_t>(targets.size()))
    throw std::invalid_argument("Shortcuts::build: inconsistent sizes");
  s.entrance_starts_ = Bitvector::from_ones(first_rows, text_size, BitvectorKind::Sparse);
  s.anchors_ = IntVector(anchors);
  s.top_tier_ = IntVector(tiers);
  s.symbol_start_ = IntVector(symbol_start);
  s.target_start_ = IntVector(target_start);
  s.targets_ = IntVector(targets);
  s.symbols_ = WaveletTree::build(symbols, static_cast<int>(max_tier) + 2, BitvectorKind::Dense);
  return s;
}

bool Shortcuts::jump(int64_t e, int64_t first, int64_t t, int64_t lo, int64_t hi, int64_t &new_lo,
                     int64_t &new_hi) const {
  const int64_t base = symbol_start_[e];
  const uint16_t tier = static_cast<uint16_t>(t);
  const int64_t earlier = symbols_.count_greater(tier, base);  // survivors in earlier entrances
  auto survivors_before = [&](int64_t row) {  // rows of the entrance before `row` that survive tier t
    return symbols_.count_greater(tier, base + (row - first)) - earlier;
  };
  const int64_t before_lo = survivors_before(lo);
  const int64_t through_hi = survivors_before(hi + 1);
  if (through_hi == before_lo) return false;  // no row of [lo, hi] survives
  const int64_t f = target(e, t);
  new_lo = f + before_lo;
  new_hi = f + through_hi - 1;
  return true;
}

void Shortcuts::save(const std::string &path) const {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  io::write(out, ell);
  io::write(out, max_tier);
  reference_.save(out);
  entrance_starts_.save(out);
  anchors_.save(out);
  top_tier_.save(out);
  symbol_start_.save(out);
  target_start_.save(out);
  targets_.save(out);
  symbols_.save(out);
  if (!out) throw std::runtime_error("cannot write " + path);
}

Shortcuts Shortcuts::load(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  Shortcuts s;
  s.ell = io::read<int64_t>(in);
  s.max_tier = io::read<int64_t>(in);
  s.reference_ = PackedDna::load(in);
  s.entrance_starts_ = Bitvector::load(in);
  s.anchors_ = IntVector::load(in);
  s.top_tier_ = IntVector::load(in);
  s.symbol_start_ = IntVector::load(in);
  s.target_start_ = IntVector::load(in);
  s.targets_ = IntVector::load(in);
  s.symbols_ = WaveletTree::load(in);
  return s;
}

}  // namespace rs
