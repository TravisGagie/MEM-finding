// Brute-force tests of MEM finding with two indexes. For random
// collections and patterns, checks that forward-backward and
// Boyer-Moore-Li, each with and without shortcuts on either index, find
// exactly the MEMs of the definition (mems.hpp), with their intervals;
// and that the reversed collection is the text of the collection file
// with its lines reversed and in reverse order.
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "build.hpp"
#include "libsais64.h"
#include "mems.hpp"
#include "text.hpp"

using namespace rs;

namespace {

int64_t checks = 0, failures = 0;

#define CHECK(cond, msg)                                                         \
  do {                                                                           \
    ++checks;                                                                    \
    if (!(cond)) {                                                               \
      if (++failures <= 30) std::cerr << "FAIL: " << msg << " (" #cond ")\n";    \
    }                                                                            \
  } while (0)

std::string random_dna(std::mt19937_64 &rng, int64_t len) {
  std::string s(static_cast<size_t>(len), 'A');
  for (char &c : s) c = "ACGT"[rng() % 4];
  return s;
}

std::string mutate(std::mt19937_64 &rng, const std::string &s, double sub, double indel) {
  std::uniform_real_distribution<double> U(0, 1);
  std::string out;
  for (char c : s) {
    double u = U(rng);
    if (u < indel / 2) continue;
    if (u < indel) out.push_back("ACGT"[rng() % 4]);
    out.push_back(U(rng) < sub ? "ACGT"[rng() % 4] : c);
  }
  return out;
}

// Intervals in the BWT of T, from its suffix array.
struct Brute {
  std::vector<uint8_t> text;
  std::vector<int64_t> sa;

  std::pair<int64_t, int64_t> interval(const std::string &s) const {
    const int64_t n = static_cast<int64_t>(sa.size());
    auto cmp = [&](int64_t row) {
      int64_t pos = sa[static_cast<size_t>(row)];
      for (size_t k = 0; k < s.size(); ++k) {
        if (pos + static_cast<int64_t>(k) >= n) return -1;
        uint8_t a = text[static_cast<size_t>(pos) + k], b = static_cast<uint8_t>(s[k]);
        if (a != b) return a < b ? -1 : 1;
      }
      return 0;
    };
    int64_t lo = 0, hi = n;
    while (lo < hi) { int64_t mid = (lo + hi) / 2; if (cmp(mid) < 0) lo = mid + 1; else hi = mid; }
    int64_t first = lo;
    hi = n;
    while (lo < hi) { int64_t mid = (lo + hi) / 2; if (cmp(mid) <= 0) lo = mid + 1; else hi = mid; }
    return {first, lo - 1};
  }
  bool occurs(const std::string &s) const {
    auto iv = interval(s);
    return iv.first <= iv.second;
  }
};

// All MEMs of P, straight from the definition.
std::vector<Mem> brute_mems(const Brute &bf, const std::string &P) {
  const int64_t m = static_cast<int64_t>(P.size());
  std::vector<Mem> out;
  for (int64_t s = 0; s < m; ++s) {
    for (int64_t e = s + 1; e <= m; ++e) {
      const std::string X = P.substr(static_cast<size_t>(s), static_cast<size_t>(e - s));
      if (!bf.occurs(X)) break;  // nor does any longer string starting at s
      if (s > 0 && bf.occurs(P.substr(static_cast<size_t>(s - 1), static_cast<size_t>(e - s + 1)))) continue;
      if (e < m && bf.occurs(P.substr(static_cast<size_t>(s), static_cast<size_t>(e - s + 1)))) continue;
      auto iv = bf.interval(X);
      out.push_back({s, e, iv.first, iv.second});
    }
  }
  return out;
}

void run_trial(std::mt19937_64 &rng, int trial) {
  const int64_t ref_len = 40 + static_cast<int64_t>(rng() % 300);
  std::string R0 = random_dna(rng, ref_len);
  if (rng() % 4 == 0) {
    for (int k = 1 + rng() % 4; k > 0; --k) R0[rng() % R0.size()] = 'N';
  }
  std::vector<std::string> seqs = {R0};
  double sub = (1 + rng() % 4) * 0.01, indel = (rng() % 3) * 0.004;
  for (int c = 1 + rng() % 8; c > 0; --c) seqs.push_back(mutate(rng, R0, sub, indel));

  const std::string path = "/tmp/rs_mems_" + std::to_string(trial);
  {
    std::ofstream f(path + ".txt");
    for (auto &s : seqs) f << s << '\n';
  }
  Brute bf;
  bf.text = read_text(path + ".txt");
  const int64_t n = static_cast<int64_t>(bf.text.size());
  bf.sa.resize(static_cast<size_t>(n));
  libsais64(bf.text.data(), bf.sa.data(), n, 0, nullptr);

  // The reversed collection is the file with its lines reversed and in
  // reverse order (tac | rev).
  const std::vector<uint8_t> rev_text = reverse_collection(bf.text);
  {
    std::ofstream f(path + ".rev.txt");
    for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) f << std::string(it->rbegin(), it->rend()) << '\n';
  }
  CHECK(read_text(path + ".rev.txt") == rev_text, "trial " << trial << ": reversed collection");

  // R spans the first sequence and, sometimes, part of the second.
  int64_t R_len = ref_len + 1 + (rng() % 2 ? static_cast<int64_t>(rng() % 30) : 0);
  R_len = std::min(R_len, n - 1);
  const std::vector<uint8_t> R = read_reference(path + ".txt", R_len);
  const std::vector<uint8_t> R_rev(R.rbegin(), R.rend());
  for (const char *ext : {".txt", ".rev.txt"}) std::remove((path + ext).c_str());

  const RLBWT text_bwt = RLBWT::from_text(bf.text), rev_bwt = RLBWT::from_text(rev_text);
  const int64_t ell = 1 + static_cast<int64_t>(rng() % 6), max_tier = static_cast<int64_t>(rng() % 4);
  const int64_t theta = rng() % 3 == 0 ? 0 : static_cast<int64_t>(rng() % 60);
  ConstructionStats stats;
  int64_t dropped = 0;
  const Shortcuts sc = prune(compute_shortcuts(text_bwt, R, ell, max_tier, stats), theta, dropped);
  const Shortcuts sc_rev = prune(compute_shortcuts(rev_bwt, R_rev, ell, max_tier, stats, true), theta, dropped);

  // The reversed index's anchors are at boundaries of R: |R| - p is a
  // multiple of ell, and of ell * 2^tau_p.
  for (int64_t e = 0; e < sc_rev.num_entrances(); ++e) {
    const int64_t b = R_len - sc_rev.anchor(e);
    const int64_t tau = sc_rev.top_tier(e);
    CHECK(b % (ell << tau) == 0 && sc_rev.anchor(e) >= (ell << tau),
          "trial " << trial << ": aligned anchor " << sc_rev.anchor(e));
  }

  const TwoIndexes plain{&text_bwt, &rev_bwt, nullptr, nullptr};
  const TwoIndexes both{&text_bwt, &rev_bwt, &sc, &sc_rev};
  const TwoIndexes only_text{&text_bwt, &rev_bwt, &sc, nullptr};
  const TwoIndexes only_rev{&text_bwt, &rev_bwt, nullptr, &sc_rev};

  for (int q = 0; q < 25; ++q) {
    const std::string &src = seqs[rng() % seqs.size()];
    size_t a = rng() % src.size();
    std::string P = src.substr(a, 1 + rng() % std::min<size_t>(120, src.size() - a));
    int kind = static_cast<int>(rng() % 4);
    if (kind == 1) P = mutate(rng, P, 0.03, 0.01);
    if (kind == 2) P = random_dna(rng, 1 + rng() % 6) + mutate(rng, P, 0.05, 0) + random_dna(rng, rng() % 6);
    if (kind == 3 && P.size() > 2) P[rng() % P.size()] = 'N';
    if (P.empty()) P = "A";

    const std::vector<Mem> all = brute_mems(bf, P);
    for (int64_t L : {int64_t{1}, 1 + static_cast<int64_t>(rng() % 20)}) {
      std::vector<Mem> expected;
      for (const Mem &x : all) if (x.end - x.start >= L) expected.push_back(x);
      for (const TwoIndexes *idx : {&plain, &both, &only_text, &only_rev}) {
        CHECK(mems_forward_backward(*idx, P, L) == expected,
              "trial " << trial << ": forward-backward, L = " << L << ", |P| = " << P.size());
        CHECK(mems_bml(*idx, P, L) == expected, "trial " << trial << ": BML, L = " << L << ", |P| = " << P.size());
      }
    }
  }
}

}  // namespace

int main(int argc, char **argv) {
  int trials = argc > 1 ? std::atoi(argv[1]) : 300;
  std::mt19937_64 rng(4242);
  for (int t = 0; t < trials; ++t) run_trial(rng, t);
  std::cout << checks << " checks, " << failures << " failed\n";
  return failures == 0 ? 0 : 1;
}
