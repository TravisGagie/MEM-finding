// Bitvectors with rank and select, used inside the RLBWT and the wavelet
// trees. Two encodings with the same interface:
//   Dense:  one bit per position, in blocks of 512 bits (8 words). For
//           each block we store the number of 1s before it and, packed in
//           one word, the number of 1s before each of its words (as in
//           Vigna's rank9); and we store the block holding every 512th 1
//           and every 512th 0. rank is two lookups and one popcount;
//           select binary-searches the blocks between two samples, then
//           finds the word from the packed counts and the bit within it.
//   Sparse: Elias-Fano. The 1 at position p is split into its low
//           low_width bits, stored in a packed array, and its high part
//           p >> low_width, stored in unary in a dense bitvector: the k-th
//           1 sets bit (p >> low_width) + k there. select is one dense
//           select; rank is one dense select0, to find the start of the
//           bucket p >> low_width, and a scan of that bucket.
#pragma once
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace rs {

enum class BitvectorKind : uint8_t { Dense = 0, Sparse = 1 };

class DenseBits {
 public:
  DenseBits() = default;
  DenseBits(const std::vector<int64_t> &ones, int64_t size);

  int64_t size() const { return size_; }
  int64_t num_ones() const { return block_rank_.back(); }
  bool get(int64_t i) const { return (words_[static_cast<size_t>(i >> 6)] >> (i & 63)) & 1; }
  int64_t rank1(int64_t pos) const;  // number of 1s in [0, pos)
  int64_t select1(int64_t k) const;  // position of the k-th 1, counting from 0
  int64_t select0(int64_t k) const;  // position of the k-th 0, counting from 0

  void save(std::ostream &out) const;
  static DenseBits load(std::istream &in);

 private:
  static constexpr int64_t kBlockBits = 512, kWordsPerBlock = 8, kSampleRate = 512;
  int64_t size_ = 0;
  std::vector<uint64_t> words_;
  std::vector<int64_t> block_rank_;       // 1s before each block; one extra entry at the end
  std::vector<uint64_t> word_rank_;       // per block: 1s before word w of the block, in bits 9(w-1) .. 9w-1, for w = 1..7
  std::vector<int64_t> select1_samples_;  // block holding the (j * kSampleRate)-th 1
  std::vector<int64_t> select0_samples_;  // block holding the (j * kSampleRate)-th 0

  int64_t zeros_before_block(int64_t b) const { return b * kBlockBits - block_rank_[static_cast<size_t>(b)]; }
  // 1s in block b before its word w, for 0 <= w < 8.
  int64_t ones_before_word(int64_t b, int64_t w) const {
    return w == 0 ? 0 : static_cast<int64_t>((word_rank_[static_cast<size_t>(b)] >> (9 * (w - 1))) & 511);
  }
};

class Bitvector {
 public:
  Bitvector() = default;

  // A bitvector of length `size` whose 1s are at the given strictly
  // increasing positions.
  static Bitvector from_ones(const std::vector<int64_t> &ones, int64_t size, BitvectorKind kind);

  // Collects bits left to right, then builds.
  class Builder {
   public:
    void push_back(bool bit) {
      if (bit) ones_.push_back(size_);
      ++size_;
    }
    Bitvector finish(BitvectorKind kind) && { return from_ones(ones_, size_, kind); }

   private:
    std::vector<int64_t> ones_;
    int64_t size_ = 0;
  };

  int64_t size() const { return size_; }
  int64_t num_ones() const { return num_ones_; }

  int64_t rank1(int64_t pos) const;  // number of 1s in [0, pos)
  int64_t select1(int64_t k) const;  // position of the k-th 1, counting from 0
  // rank1(pos), and in `last` the position of the last 1 before pos (if
  // there is one), with less work than rank1 and select1 separately.
  int64_t rank1_and_predecessor(int64_t pos, int64_t &last) const;

  void save(std::ostream &out) const;
  static Bitvector load(std::istream &in);

 private:
  BitvectorKind kind_ = BitvectorKind::Dense;
  int64_t size_ = 0;
  int64_t num_ones_ = 0;
  DenseBits bits_;  // Dense: the bitvector itself. Sparse: the high parts, in unary.
  int low_width_ = 0;          // Sparse only
  std::vector<uint64_t> low_;  // Sparse only: the low parts, packed

  uint64_t low(int64_t k) const;
};

}  // namespace rs
