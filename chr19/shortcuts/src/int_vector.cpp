#include "int_vector.hpp"

#include <algorithm>
#include <stdexcept>

#include "io.hpp"

namespace rs {

IntVector::IntVector(const std::vector<int64_t> &values, int width) : size_(static_cast<int64_t>(values.size())) {
  if (width == 0) {
    int64_t largest = 0;
    for (int64_t v : values) {
      if (v < 0) throw std::invalid_argument("IntVector values must be non-negative");
      largest = std::max(largest, v);
    }
    while (width < 63 && (largest >> width) != 0) ++width;
  }
  width_ = width;
  mask_ = width_ >= 64 ? ~uint64_t{0} : (uint64_t{1} << width_) - 1;
  words_.assign(static_cast<size_t>((size_ * width_ + 63) / 64) + 1, 0);
  for (int64_t i = 0; i < size_; ++i) {
    const uint64_t v = static_cast<uint64_t>(values[static_cast<size_t>(i)]) & mask_;
    const uint64_t bit = static_cast<uint64_t>(i) * static_cast<uint64_t>(width_), w = bit >> 6, off = bit & 63;
    if (width_ == 0) continue;
    words_[w] |= v << off;
    if (off + static_cast<uint64_t>(width_) > 64) words_[w + 1] |= v >> (64 - off);
  }
}

void IntVector::save(std::ostream &out) const {
  io::write(out, size_);
  io::write(out, static_cast<int64_t>(width_));
  io::write_vector(out, words_);
}

IntVector IntVector::load(std::istream &in) {
  IntVector v;
  v.size_ = io::read<int64_t>(in);
  v.width_ = static_cast<int>(io::read<int64_t>(in));
  v.mask_ = v.width_ >= 64 ? ~uint64_t{0} : (uint64_t{1} << v.width_) - 1;
  v.words_ = io::read_vector<uint64_t>(in);
  return v;
}

PackedDna::PackedDna(const std::vector<uint8_t> &s) : size_(static_cast<int64_t>(s.size())) {
  words_.assign(static_cast<size_t>((size_ + 31) / 32), 0);
  for (int64_t i = 0; i < size_; ++i) {
    uint64_t code = 0;
    switch (s[static_cast<size_t>(i)]) {
      case 'C': code = 1; break;
      case 'G': code = 2; break;
      case 'T': code = 3; break;
      default: code = 0;
    }
    words_[static_cast<size_t>(i >> 5)] |= code << (2 * (i & 31));
  }
}

void PackedDna::save(std::ostream &out) const {
  io::write(out, size_);
  io::write_vector(out, words_);
}

PackedDna PackedDna::load(std::istream &in) {
  PackedDna d;
  d.size_ = io::read<int64_t>(in);
  d.words_ = io::read_vector<uint64_t>(in);
  return d;
}

}  // namespace rs
