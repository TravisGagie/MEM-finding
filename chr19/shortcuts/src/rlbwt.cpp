#include "rlbwt.hpp"

#include <fstream>
#include <stdexcept>

#include "io.hpp"
#include "libsais64.h"
#include "text.hpp"

namespace rs {

// The bitvectors over rows are sparse, since runs are long; the wavelet
// tree over the few run characters uses dense bitvectors.
constexpr BitvectorKind kRowBitvectors = BitvectorKind::Sparse;
constexpr BitvectorKind kRunCharBitvectors = BitvectorKind::Dense;

void RLBWT::set_codes() {
  for (int c = 0; c < 256; ++c) code_[c] = -1;
  for (size_t a = 0; a < alphabet_.size(); ++a) code_[alphabet_[a]] = static_cast<int16_t>(a);
}

// Builds an RLBWT from its characters, given row by row. `count[c]` must
// be the number of occurrences of c.
class RLBWT::Builder {
 public:
  Builder(const int64_t count[256]) {
    for (int c = 0; c < 256; ++c) {
      if (count[c] > 0) r_.alphabet_.push_back(static_cast<uint8_t>(c));
    }
    r_.set_codes();
    r_.C_.assign(r_.alphabet_.size() + 1, 0);
    for (size_t a = 0; a < r_.alphabet_.size(); ++a) r_.C_[a + 1] = r_.C_[a] + count[r_.alphabet_[a]];
    r_.n_ = r_.C_.back();
  }

  void push(uint8_t c) {
    int a = r_.code_[c];
    bool new_run = a != previous_;
    starts_.push_back(new_run);
    if (new_run) r_.run_char_.push_back(static_cast<uint16_t>(a));
    previous_ = a;
  }

  RLBWT finish() && {
    RLBWT &r = r_;
    r.run_starts_ = std::move(starts_).finish(kRowBitvectors);
    if (r.run_starts_.size() != r.n_) throw std::runtime_error("RLBWT: wrong number of rows");
    int sigma = static_cast<int>(r.alphabet_.size());
    int64_t R = r.num_runs();
    r.run_char_wt_ = WaveletTree::build(r.run_char_, sigma, kRunCharBitvectors);

    r.runs_before_.assign(static_cast<size_t>(sigma) + 1, 0);
    for (uint16_t a : r.run_char_) ++r.runs_before_[a + 1];
    for (int a = 0; a < sigma; ++a) r.runs_before_[a + 1] += r.runs_before_[a];

    // The k-th run of character a, of length len, starts in F at
    // C[a] + (total length of the runs of a before it).
    std::vector<int64_t> sorted_starts(static_cast<size_t>(R));
    std::vector<int64_t> runs_so_far(static_cast<size_t>(sigma), 0), length_so_far(static_cast<size_t>(sigma), 0);
    for (int64_t run = 0; run < R; ++run) {
      int a = r.run_char_[static_cast<size_t>(run)];
      int64_t len = r.run_start(run + 1) - r.run_start(run);
      sorted_starts[static_cast<size_t>(r.runs_before_[a] + runs_so_far[a]++)] = r.C_[a] + length_so_far[a];
      length_so_far[a] += len;
    }
    r.sorted_run_starts_ = Bitvector::from_ones(sorted_starts, r.n_, kRowBitvectors);
    return std::move(r_);
  }

 private:
  RLBWT r_;
  Bitvector::Builder starts_;
  int previous_ = -1;
};

RLBWT RLBWT::from_text(const std::vector<uint8_t> &text) {
  int64_t n = static_cast<int64_t>(text.size());
  if (n == 0 || text.back() != kTerminator) throw std::runtime_error("the text must end with the terminator");
  std::vector<int64_t> sa(static_cast<size_t>(n));
  if (libsais64(text.data(), sa.data(), n, 0, nullptr) != 0) throw std::runtime_error("libsais64 failed");
  int64_t count[256] = {0};
  for (uint8_t c : text) ++count[c];
  Builder b(count);
  for (int64_t i = 0; i < n; ++i) {
    int64_t pos = sa[static_cast<size_t>(i)];
    b.push(text[static_cast<size_t>(pos == 0 ? n - 1 : pos - 1)]);
  }
  return std::move(b).finish();
}

RLBWT RLBWT::from_bwt_file(const std::string &path, uint8_t separator) {
  const size_t kBuffer = size_t{1} << 24;
  std::vector<char> buf(kBuffer);
  auto for_each_byte = [&](auto f) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    while (in) {
      in.read(buf.data(), static_cast<std::streamsize>(kBuffer));
      std::streamsize got = in.gcount();
      for (std::streamsize i = 0; i < got; ++i) f(static_cast<uint8_t>(buf[static_cast<size_t>(i)]));
    }
  };

  // First pass: count the characters and find the end-of-string byte.
  int64_t count[256] = {0};
  for_each_byte([&](uint8_t c) { ++count[c]; });
  int eos = -1;
  for (int c = 0; c < 'A'; ++c) {
    if (c == separator || count[c] == 0) continue;
    if (count[c] != 1 || eos != -1) throw std::runtime_error("cannot identify the end-of-string byte in " + path);
    eos = c;
  }
  if (eos == -1) throw std::runtime_error("no end-of-string byte in " + path);

  // Map separator -> #, end-of-string -> $. Both are smaller than every
  // other character present, so the order of the characters is kept.
  uint8_t map[256];
  for (int c = 0; c < 256; ++c) map[c] = static_cast<uint8_t>(c);
  map[separator] = kSeparator;
  map[eos] = kTerminator;
  int64_t mapped_count[256] = {0};
  for (int c = 0; c < 256; ++c) mapped_count[map[c]] += count[c];

  // Second pass: build.
  Builder b(mapped_count);
  for_each_byte([&](uint8_t c) { b.push(map[c]); });
  return std::move(b).finish();
}

int64_t RLBWT::rank(uint8_t c, int64_t i) const {
  int a = code_[c];
  if (a < 0 || i <= 0) return 0;
  if (i >= n_) return C_[a + 1] - C_[a];
  int64_t run = run_containing(i - 1);
  int64_t k = run_char_wt_.rank(static_cast<uint16_t>(a), run);  // runs of c before this run
  // In F, the runs sorted by character, the occurrences of c in those k runs
  // end where run number runs_before_[a] + k starts: the (k+1)-st run of c,
  // or, if there is none, the first run of the next character, which starts
  // at C[a+1]. Only for the last character is there no such run.
  int64_t before = runs_before_[a] + k < num_runs() ? sorted_run_starts_.select1(runs_before_[a] + k) - C_[a]
                                                     : C_[a + 1] - C_[a];
  if (run_char_[static_cast<size_t>(run)] == a) before += i - run_start(run);
  return before;
}

void RLBWT::backward_step(uint8_t c, int64_t &lo, int64_t &hi) const {
  int a = code_[c];
  if (a < 0) {
    lo = 1;
    hi = 0;
    return;
  }
  int64_t new_lo = C_[a] + rank(c, lo);
  hi = C_[a] + rank(c, hi + 1) - 1;
  lo = new_lo;
}

void RLBWT::save(const std::string &path) const {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  io::write(out, n_);
  io::write_vector(out, alphabet_);
  io::write_vector(out, C_);
  io::write_vector(out, runs_before_);
  io::write_vector(out, run_char_);
  run_starts_.save(out);
  run_char_wt_.save(out);
  sorted_run_starts_.save(out);
  if (!out) throw std::runtime_error("cannot write " + path);
}

RLBWT RLBWT::load(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  RLBWT r;
  r.n_ = io::read<int64_t>(in);
  r.alphabet_ = io::read_vector<uint8_t>(in);
  r.C_ = io::read_vector<int64_t>(in);
  r.runs_before_ = io::read_vector<int64_t>(in);
  r.run_char_ = io::read_vector<uint16_t>(in);
  r.run_starts_ = Bitvector::load(in);
  r.run_char_wt_ = WaveletTree::load(in);
  r.sorted_run_starts_ = Bitvector::load(in);
  r.set_codes();
  return r;
}

}  // namespace rs
