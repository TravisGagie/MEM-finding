// Fixed-width packed integer arrays, for the per-entrance data of the
// shortcuts, and 2-bit packed DNA, for the reference.
#pragma once
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace rs {

// Non-negative integers of `width` bits each, packed into 64-bit words.
class IntVector {
 public:
  IntVector() = default;
  // width = 0 chooses the smallest width that holds the largest value.
  explicit IntVector(const std::vector<int64_t> &values, int width = 0);

  int64_t size() const { return size_; }
  int width() const { return width_; }
  int64_t operator[](int64_t i) const {
    if (width_ == 0) return 0;
    const uint64_t bit = static_cast<uint64_t>(i) * static_cast<uint64_t>(width_);
    const uint64_t w = bit >> 6, off = bit & 63;
    uint64_t v = words_[w] >> off;
    if (off + static_cast<uint64_t>(width_) > 64) v |= words_[w + 1] << (64 - off);
    return static_cast<int64_t>(v & mask_);
  }

  void save(std::ostream &out) const;
  static IntVector load(std::istream &in);

 private:
  int64_t size_ = 0;
  int width_ = 0;
  uint64_t mask_ = 0;
  std::vector<uint64_t> words_;
};

// A string over A, C, G, T in 2 bits per character. Other characters are
// stored as A; the shortcuts only compare characters of tier strings,
// which are all A/C/G/T.
class PackedDna {
 public:
  PackedDna() = default;
  explicit PackedDna(const std::vector<uint8_t> &s);

  int64_t size() const { return size_; }
  uint8_t operator[](int64_t i) const {
    static constexpr uint8_t kBases[4] = {'A', 'C', 'G', 'T'};
    return kBases[(words_[static_cast<uint64_t>(i) >> 5] >> (2 * (static_cast<uint64_t>(i) & 31))) & 3];
  }

  void save(std::ostream &out) const;
  static PackedDna load(std::istream &in);

 private:
  int64_t size_ = 0;
  std::vector<uint64_t> words_;
};

}  // namespace rs
