// Checks rank and select of both bitvector encodings, and select0 of the
// plain one, against the bits themselves, for random densities and sizes.
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "bitvector.hpp"

using namespace rs;

int main() {
  std::mt19937_64 rng(1);
  int64_t checks = 0, failures = 0;
  const double densities[] = {0.0005, 0.01, 0.1, 0.5, 0.99};
  for (int trial = 0; trial < 520; ++trial) {
    // The last 20 trials are larger, with the density changing every
    // 50,000 bits, so that the blocks between select samples vary widely.
    const bool large = trial >= 500;
    const int64_t n = large ? 1000000 + static_cast<int64_t>(rng() % 1000000) : 1 + static_cast<int64_t>(rng() % 20000);
    double density = densities[trial % 5];
    std::vector<int64_t> ones;
    for (int64_t i = 0; i < n; ++i) {
      if (large && i % 50000 == 0) density = densities[rng() % 5];
      if (std::uniform_real_distribution<double>(0, 1)(rng) < density) ones.push_back(i);
    }

    for (BitvectorKind kind : {BitvectorKind::Dense, BitvectorKind::Sparse}) {
      Bitvector bv = Bitvector::from_ones(ones, n, kind);
      size_t next = 0;  // index of the next 1
      for (int64_t p = 0; p <= n; ++p) {
        ++checks;
        if (bv.rank1(p) != static_cast<int64_t>(next)) ++failures;
        int64_t last = -1;
        ++checks;
        if (bv.rank1_and_predecessor(p, last) != static_cast<int64_t>(next) || (next > 0 && last != ones[next - 1]))
          ++failures;
        if (next < ones.size() && ones[next] == p) {
          ++checks;
          if (bv.select1(static_cast<int64_t>(next)) != p) ++failures;
          ++next;
        }
      }
    }

    DenseBits bits(ones, n);
    size_t next = 0;
    int64_t zeros = 0;
    for (int64_t p = 0; p < n; ++p) {
      if (next < ones.size() && ones[next] == p) {
        ++next;
        continue;
      }
      ++checks;
      if (bits.select0(zeros++) != p) ++failures;
    }
  }
  // Bitvectors with no ones; a sparse one must also be small when saved.
  for (int64_t n : {int64_t{1}, int64_t{2}, int64_t{63}, int64_t{64}, int64_t{65}, int64_t{1000}, int64_t{1} << 20,
                    int64_t{3000000}}) {
    for (BitvectorKind kind : {BitvectorKind::Dense, BitvectorKind::Sparse}) {
      Bitvector bv = Bitvector::from_ones({}, n, kind);
      for (int64_t p : {int64_t{0}, n / 3, n - 1, n}) {
        int64_t last = -1;
        ++checks;
        if (bv.rank1(p) != 0 || bv.rank1_and_predecessor(p, last) != 0) ++failures;
      }
      if (kind == BitvectorKind::Sparse) {
        std::ostringstream out;
        bv.save(out);
        ++checks;
        if (out.str().size() > 1000) ++failures;
      }
    }
  }
  std::cout << checks << " checks, " << failures << " failed\n";
  return failures == 0 ? 0 : 1;
}
