#include "bitvector.hpp"

#if defined(__BMI2__)
#include <immintrin.h>
#endif
#include <istream>
#include <ostream>

#include "io.hpp"

namespace rs {

namespace {

// Position of the k-th 1 (counting from 0) in a word with more than k 1s.
inline int select_in_word(uint64_t word, int64_t k) {
#if defined(__BMI2__)
  return __builtin_ctzll(_pdep_u64(uint64_t{1} << k, word));
#else
  for (int byte = 0;; ++byte) {
    const uint64_t w = (word >> (8 * byte)) & 0xff;
    const int c = __builtin_popcountll(w);
    if (k < c) {
      uint64_t x = w;
      for (; k > 0; --k) x &= x - 1;  // drop the lowest 1s
      return 8 * byte + __builtin_ctzll(x);
    }
    k -= c;
  }
#endif
}

}  // namespace

// ---------------------------------------------------------------- DenseBits

DenseBits::DenseBits(const std::vector<int64_t> &ones, int64_t size) : size_(size) {
  const int64_t blocks = (size + kBlockBits - 1) / kBlockBits;
  words_.assign(static_cast<size_t>(blocks * kWordsPerBlock), 0);
  for (int64_t p : ones) words_[static_cast<size_t>(p >> 6)] |= uint64_t{1} << (p & 63);

  block_rank_.assign(static_cast<size_t>(blocks) + 1, 0);
  word_rank_.assign(static_cast<size_t>(blocks), 0);
  for (int64_t b = 0; b < blocks; ++b) {
    int64_t c = 0;
    uint64_t packed = 0;
    for (int64_t w = 0; w < kWordsPerBlock; ++w) {
      if (w > 0) packed |= static_cast<uint64_t>(c) << (9 * (w - 1));
      c += __builtin_popcountll(words_[static_cast<size_t>(b * kWordsPerBlock + w)]);
    }
    word_rank_[static_cast<size_t>(b)] = packed;
    block_rank_[static_cast<size_t>(b) + 1] = block_rank_[static_cast<size_t>(b)] + c;
  }
  // Every kSampleRate-th 1 and 0: the block holding it. (Zeros in the
  // padding after `size` are counted but never asked for.)
  for (int64_t b = 0; b < blocks; ++b) {
    for (int64_t k = (block_rank_[b] + kSampleRate - 1) / kSampleRate * kSampleRate; k < block_rank_[b + 1];
         k += kSampleRate) {
      select1_samples_.push_back(b);
    }
    for (int64_t k = (zeros_before_block(b) + kSampleRate - 1) / kSampleRate * kSampleRate;
         k < zeros_before_block(b + 1); k += kSampleRate) {
      select0_samples_.push_back(b);
    }
  }
}

int64_t DenseBits::rank1(int64_t pos) const {
  if (pos <= 0) return 0;
  if (pos >= size_) return num_ones();
  const int64_t b = pos / kBlockBits, w = (pos >> 6) & 7;
  const uint64_t below = (uint64_t{1} << (pos & 63)) - 1;
  return block_rank_[static_cast<size_t>(b)] + ones_before_word(b, w) +
         __builtin_popcountll(words_[static_cast<size_t>(pos >> 6)] & below);
}

int64_t DenseBits::select1(int64_t k) const {
  // The block holding the k-th 1 lies between the sampled blocks of the
  // samples around k: the last b there with block_rank_[b] <= k.
  const size_t j = static_cast<size_t>(k / kSampleRate);
  int64_t lo = select1_samples_[j];
  int64_t hi = j + 1 < select1_samples_.size() ? select1_samples_[j + 1] : static_cast<int64_t>(block_rank_.size()) - 2;
  while (lo < hi) {
    const int64_t mid = (lo + hi + 1) / 2;
    if (block_rank_[static_cast<size_t>(mid)] <= k) lo = mid; else hi = mid - 1;
  }
  k -= block_rank_[static_cast<size_t>(lo)];
  int64_t w = 0;  // the number of words of the block with at most k 1s before them, minus 1
  for (int64_t i = 1; i < kWordsPerBlock; ++i) w += ones_before_word(lo, i) <= k;
  return (lo * kWordsPerBlock + w) * 64 + select_in_word(words_[static_cast<size_t>(lo * kWordsPerBlock + w)],
                                                         k - ones_before_word(lo, w));
}

int64_t DenseBits::select0(int64_t k) const {
  const size_t j = static_cast<size_t>(k / kSampleRate);
  int64_t lo = select0_samples_[j];
  int64_t hi = j + 1 < select0_samples_.size() ? select0_samples_[j + 1] : static_cast<int64_t>(block_rank_.size()) - 2;
  while (lo < hi) {
    const int64_t mid = (lo + hi + 1) / 2;
    if (zeros_before_block(mid) <= k) lo = mid; else hi = mid - 1;
  }
  k -= zeros_before_block(lo);
  int64_t w = 0;
  for (int64_t i = 1; i < kWordsPerBlock; ++i) w += 64 * i - ones_before_word(lo, i) <= k;
  return (lo * kWordsPerBlock + w) * 64 + select_in_word(~words_[static_cast<size_t>(lo * kWordsPerBlock + w)],
                                                         k - (64 * w - ones_before_word(lo, w)));
}

void DenseBits::save(std::ostream &out) const {
  io::write(out, size_);
  io::write_vector(out, words_);
  io::write_vector(out, block_rank_);
  io::write_vector(out, word_rank_);
  io::write_vector(out, select1_samples_);
  io::write_vector(out, select0_samples_);
}

DenseBits DenseBits::load(std::istream &in) {
  DenseBits d;
  d.size_ = io::read<int64_t>(in);
  d.words_ = io::read_vector<uint64_t>(in);
  d.block_rank_ = io::read_vector<int64_t>(in);
  d.word_rank_ = io::read_vector<uint64_t>(in);
  d.select1_samples_ = io::read_vector<int64_t>(in);
  d.select0_samples_ = io::read_vector<int64_t>(in);
  return d;
}

// ---------------------------------------------------------------- Bitvector

Bitvector Bitvector::from_ones(const std::vector<int64_t> &ones, int64_t size, BitvectorKind kind) {
  Bitvector bv;
  bv.kind_ = kind;
  bv.size_ = size;
  bv.num_ones_ = static_cast<int64_t>(ones.size());
  if (kind == BitvectorKind::Dense) {
    bv.bits_ = DenseBits(ones, size);
    return bv;
  }
  // Elias-Fano with low_width = floor(log2(size / m)).
  const int64_t m = bv.num_ones_;
  int width = 0;
  if (m > 0 && size > m) {
    while ((int64_t{1} << (width + 1)) <= size / m) ++width;
  } else if (m == 0) {  // no ones: one bucket, so that the high part takes O(1) bits
    while (width < 62 && ((size - 1) >> width) > 0) ++width;
  }
  bv.low_width_ = width;
  bv.low_.assign(static_cast<size_t>((m * width + 63) / 64) + 1, 0);
  std::vector<int64_t> high_ones;
  high_ones.reserve(static_cast<size_t>(m));
  for (int64_t k = 0; k < m; ++k) {
    const int64_t p = ones[static_cast<size_t>(k)];
    if (width > 0) {  // pack the low bits of p into field k
      const uint64_t v = static_cast<uint64_t>(p) & ((uint64_t{1} << width) - 1);
      const int64_t bit = k * width;
      bv.low_[static_cast<size_t>(bit >> 6)] |= v << (bit & 63);
      if ((bit & 63) + width > 64) bv.low_[static_cast<size_t>(bit >> 6) + 1] |= v >> (64 - (bit & 63));
    }
    high_ones.push_back((p >> width) + k);
  }
  const int64_t buckets = size <= 0 ? 1 : ((size - 1) >> width) + 1;
  bv.bits_ = DenseBits(high_ones, buckets + m);
  return bv;
}

uint64_t Bitvector::low(int64_t k) const {
  if (low_width_ == 0) return 0;
  const int64_t bit = k * low_width_;
  uint64_t v = low_[static_cast<size_t>(bit >> 6)] >> (bit & 63);
  if ((bit & 63) + low_width_ > 64) v |= low_[static_cast<size_t>(bit >> 6) + 1] << (64 - (bit & 63));
  return v & ((uint64_t{1} << low_width_) - 1);
}

int64_t Bitvector::select1(int64_t k) const {
  if (kind_ == BitvectorKind::Dense) return bits_.select1(k);
  const int64_t high = bits_.select1(k) - k;
  return (high << low_width_) | static_cast<int64_t>(low(k));
}

int64_t Bitvector::rank1(int64_t pos) const {
  if (kind_ == BitvectorKind::Dense) return bits_.rank1(pos);
  if (pos <= 0) return 0;
  if (pos >= size_) return num_ones_;
  // The 1s in buckets before pos's bucket come before the bucket's 0 in
  // unary; then scan the 1s of pos's bucket.
  const int64_t bucket = pos >> low_width_;
  const uint64_t pos_low = static_cast<uint64_t>(pos) & ((uint64_t{1} << low_width_) - 1);
  int64_t h = bucket == 0 ? 0 : bits_.select0(bucket - 1) + 1;  // where the bucket starts in unary
  int64_t k = h - bucket;                                       // 1s in earlier buckets
  while (h < bits_.size() && bits_.get(h) && low(k) < pos_low) {
    ++h;
    ++k;
  }
  return k;
}

int64_t Bitvector::rank1_and_predecessor(int64_t pos, int64_t &last) const {
  if (kind_ == BitvectorKind::Dense || pos >= size_) {
    const int64_t k = rank1(pos);
    if (k > 0) last = select1(k - 1);
    return k;
  }
  if (pos <= 0) return 0;
  const int64_t bucket = pos >> low_width_;
  const uint64_t pos_low = static_cast<uint64_t>(pos) & ((uint64_t{1} << low_width_) - 1);
  int64_t h = bucket == 0 ? 0 : bits_.select0(bucket - 1) + 1;
  int64_t k = h - bucket;
  while (h < bits_.size() && bits_.get(h) && low(k) < pos_low) {
    ++h;
    ++k;
  }
  if (k > 0) {
    // The (k-1)-th 1 is the last 1 before position h in unary; between
    // them are only the 0s that close empty buckets.
    int64_t x = h - 1;
    while (!bits_.get(x)) --x;
    last = ((x - (k - 1)) << low_width_) | static_cast<int64_t>(low(k - 1));
  }
  return k;
}

void Bitvector::save(std::ostream &out) const {
  io::write(out, static_cast<uint8_t>(kind_));
  io::write(out, size_);
  io::write(out, num_ones_);
  bits_.save(out);
  io::write(out, static_cast<int64_t>(low_width_));
  io::write_vector(out, low_);
}

Bitvector Bitvector::load(std::istream &in) {
  Bitvector bv;
  bv.kind_ = static_cast<BitvectorKind>(io::read<uint8_t>(in));
  bv.size_ = io::read<int64_t>(in);
  bv.num_ones_ = io::read<int64_t>(in);
  bv.bits_ = DenseBits::load(in);
  bv.low_width_ = static_cast<int>(io::read<int64_t>(in));
  bv.low_ = io::read_vector<uint64_t>(in);
  return bv;
}

}  // namespace rs
