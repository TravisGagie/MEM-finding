// Brute-force tests. For random collections (a reference R and mutated
// copies, sometimes with internal repeats in R and with N characters),
// checks against definitions computed directly from the suffix array:
//   - plain backward search counts correctly;
//   - an RLBWT read from a BWT file is identical to one built from the text;
//   - the anchors kept, their entrances and top tiers (Section 4);
//   - each row's symbol and each f_{p,t} (Sections 5 and 6);
// and, after pruning with random thresholds:
//   - the entrance found for an interval is the one containing it;
//   - search with shortcuts gives the same result as plain search;
//   - every report entry's certificate holds;
//   - saving and loading changes nothing.
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "build.hpp"
#include "libsais64.h"
#include "search.hpp"
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
    if (u < indel / 2) continue;                    // deletion
    if (u < indel) out.push_back("ACGT"[rng() % 4]);  // insertion
    out.push_back(U(rng) < sub ? "ACGT"[rng() % 4] : c);
  }
  return out;
}

struct Brute {
  std::vector<uint8_t> text;
  std::vector<int64_t> sa;

  // The interval of rows whose suffixes start with s.
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
};

void run_trial(std::mt19937_64 &rng, int trial) {
  // ---- the collection; R is its first sequence ----
  const int64_t ref_len = 40 + static_cast<int64_t>(rng() % 400);
  std::string R = random_dna(rng, ref_len);
  if (rng() % 3 == 0 && ref_len > 60) {  // an internal repeat, so that windows collide
    size_t len = 10 + rng() % 30, from = rng() % (ref_len - len), to = rng() % (ref_len - len);
    R.replace(to, len, R.substr(from, len));
  }
  if (rng() % 4 == 0) {
    for (int k = 1 + rng() % 5; k > 0; --k) R[rng() % R.size()] = 'N';
  }
  std::vector<std::string> seqs = {R};
  double sub = (rng() % 4) * 0.01, indel = (rng() % 3) * 0.004;
  for (int c = 1 + rng() % 8; c > 0; --c) seqs.push_back(mutate(rng, R, sub, indel));

  const std::string path = "/tmp/rs_test_" + std::to_string(trial);
  {
    std::ofstream f(path + ".txt");
    for (auto &s : seqs) f << s << '\n';
  }
  Brute bf;
  bf.text = read_text(path + ".txt");
  const int64_t n = static_cast<int64_t>(bf.text.size());
  bf.sa.resize(static_cast<size_t>(n));
  libsais64(bf.text.data(), bf.sa.data(), n, 0, nullptr);
  RLBWT rlbwt = RLBWT::from_text(bf.text);

  // ---- an RLBWT read from a BWT file, with '\n' for # and byte 2 for $ ----
  {
    std::ofstream f(path + ".bwt", std::ios::binary);
    for (int64_t i = 0; i < n; ++i) {
      int64_t pos = bf.sa[static_cast<size_t>(i)];
      uint8_t c = bf.text[static_cast<size_t>(pos == 0 ? n - 1 : pos - 1)];
      f.put(static_cast<char>(c == kSeparator ? '\n' : c == kTerminator ? 2 : c));
    }
  }
  RLBWT from_file = RLBWT::from_bwt_file(path + ".bwt");
  rlbwt.save(path + ".a");
  from_file.save(path + ".b");
  {
    std::ifstream a(path + ".a", std::ios::binary), b(path + ".b", std::ios::binary);
    std::string bytes_a((std::istreambuf_iterator<char>(a)), {}), bytes_b((std::istreambuf_iterator<char>(b)), {});
    CHECK(bytes_a == bytes_b, "trial " << trial << ": RLBWT from BWT file");
  }
  for (const char *ext : {".bwt", ".a", ".b"}) std::remove((path + ext).c_str());

  // Plain search returns the longest suffix of P that occurs, and its interval.
  for (int q = 0; q < 20; ++q) {
    std::string P = random_dna(rng, 1 + rng() % 4) + seqs[rng() % seqs.size()].substr(rng() % 20, 1 + rng() % 15);
    SearchResult r = search_plain(rlbwt, P);
    auto iv = bf.interval(P.substr(P.size() - static_cast<size_t>(r.matched)));
    bool longest = true;
    if (!r.found()) {  // one character more must not occur
      auto longer = bf.interval(P.substr(P.size() - static_cast<size_t>(r.matched) - 1));
      longest = longer.second < longer.first;
    }
    CHECK(iv.first == r.lo && iv.second == r.hi && longest, "trial " << trial << ": plain search");
  }

  // ---- shortcuts against their definitions ----
  const int64_t ell = 1 + static_cast<int64_t>(rng() % 6);
  const int64_t max_tier = static_cast<int64_t>(rng() % 4);
  std::vector<uint8_t> refv = read_reference(path + ".txt", ref_len);
  std::remove((path + ".txt").c_str());
  ConstructionStats stats;
  UnprunedShortcuts u = compute_shortcuts(rlbwt, refv, ell, max_tier, stats);

  auto usable = [&](int64_t from, int64_t to) {
    for (int64_t k = from; k < to; ++k)
      if (!is_acgt(static_cast<uint8_t>(R[static_cast<size_t>(k)]))) return false;
    return true;
  };
  std::map<std::string, int> copies;
  for (int64_t p = ell; p + ell <= ref_len; p += ell)
    if (usable(p, p + ell)) ++copies[R.substr(static_cast<size_t>(p), static_cast<size_t>(ell))];
  struct Expected { int64_t p, lo, hi, tau; };
  std::vector<Expected> expected;
  for (int64_t p = ell; p + ell <= ref_len; p += ell) {
    std::string w = R.substr(static_cast<size_t>(p), static_cast<size_t>(ell));
    if (!usable(p, p + ell) || copies[w] != 1) continue;
    int64_t tau = 0;
    while (tau < max_tier && p % (ell << (tau + 1)) == 0) ++tau;
    while (tau >= 0 && !usable(p - (ell << tau), p)) --tau;
    if (tau < 0) continue;
    auto iv = bf.interval(w);
    expected.push_back({p, iv.first, iv.second, tau});
  }
  std::sort(expected.begin(), expected.end(), [](const Expected &a, const Expected &b) { return a.lo < b.lo; });
  CHECK(u.entrances.size() == expected.size(), "trial " << trial << ": number of entrances");

  for (size_t e = 0; e < std::min(u.entrances.size(), expected.size()); ++e) {
    const Entrance &en = u.entrances[e];
    const Expected &x = expected[e];
    bool ok = en.anchor == x.p && en.first_row == x.lo && en.last_row == x.hi && u.top_tier[e] == x.tau;
    CHECK(ok, "trial " << trial << ": entrance and top tier of anchor " << x.p);
    if (!ok) continue;
    for (int64_t row = en.first_row; row <= en.last_row; ++row) {
      // Number of characters before the row's suffix that agree with R before p.
      int64_t pos = bf.sa[static_cast<size_t>(row)], agree = 0;
      while (agree < (ell << x.tau) && pos - 1 - agree >= 0 &&
             bf.text[static_cast<size_t>(pos - 1 - agree)] == static_cast<uint8_t>(R[static_cast<size_t>(x.p - 1 - agree)]))
        ++agree;
      int64_t sym = 0;
      for (int64_t t = 0; t <= x.tau; ++t)
        if ((ell << t) <= agree) sym = t + 1;
      CHECK(u.symbols[static_cast<size_t>(u.symbol_start[e] + row - en.first_row)] == sym,
            "trial " << trial << ": symbol of row " << row);
    }
    for (int64_t t = 0; t <= x.tau; ++t) {
      auto iv = bf.interval(R.substr(static_cast<size_t>(x.p - (ell << t)), static_cast<size_t>((ell << t) + ell)));
      CHECK(u.target_first_row[static_cast<size_t>(u.target_start[e] + t)] == iv.first,
            "trial " << trial << ": f_{p,t} for anchor " << x.p << ", tier " << t);
    }
  }

  // ---- search ----
  int64_t dropped = 0;
  Shortcuts sc = prune(u, rng() % 3 == 0 ? 0 : static_cast<int64_t>(rng() % 60), dropped);
  sc.save(path + ".shortcuts");
  Shortcuts loaded = Shortcuts::load(path + ".shortcuts");
  std::remove((path + ".shortcuts").c_str());

  // The packed structure holds the unpruned data of the anchors kept.
  {
    std::map<int64_t, size_t> by_anchor;
    for (size_t e = 0; e < u.entrances.size(); ++e) by_anchor[u.entrances[e].anchor] = e;
    for (int64_t e = 0; e < sc.num_entrances(); ++e) {
      auto it = by_anchor.find(sc.anchor(e));
      bool ok = it != by_anchor.end();
      if (ok) {
        const size_t k = it->second;
        ok = sc.first_row(e) == u.entrances[k].first_row && sc.last_row(e) == u.entrances[k].last_row &&
             sc.top_tier(e) <= u.top_tier[k] && sc.rows(e) == u.symbol_start[k + 1] - u.symbol_start[k];
        for (int64_t t = 0; ok && t <= sc.top_tier(e); ++t) ok = sc.target(e, t) == u.target_first_row[u.target_start[k] + t];
        for (int64_t i = 0; ok && i < static_cast<int64_t>(R.size()); ++i)
          if (is_acgt(static_cast<uint8_t>(R[static_cast<size_t>(i)]))) ok = loaded.reference_char(i) == static_cast<uint8_t>(R[static_cast<size_t>(i)]);
      }
      CHECK(ok, "trial " << trial << ": packed entrance " << e);
    }
  }

  // find_entrance against a scan of the entrances, on random intervals and
  // on the entrances themselves and their neighbours.
  const int64_t rows = rlbwt.size();
  auto expected_entrance = [&](int64_t lo, int64_t hi) -> int64_t {
    for (int64_t e = 0; e < sc.num_entrances(); ++e) {
      if (sc.first_row(e) <= lo && hi <= sc.last_row(e)) return e;
    }
    return -1;
  };
  std::vector<std::pair<int64_t, int64_t>> intervals;
  for (int q = 0; q < 100; ++q) {
    int64_t lo = static_cast<int64_t>(rng() % static_cast<uint64_t>(rows));
    intervals.emplace_back(lo, std::min(rows - 1, lo + static_cast<int64_t>(rng() % 20)));
  }
  for (int64_t e = 0; e < sc.num_entrances(); ++e) {
    const int64_t first = sc.first_row(e), last = sc.last_row(e);
    intervals.emplace_back(first, last);
    if (first > 0) intervals.emplace_back(first - 1, last);
    if (last + 1 < rows) intervals.emplace_back(first, last + 1);
  }
  for (auto [lo, hi] : intervals) {
    CHECK(sc.find_entrance(lo, hi) == expected_entrance(lo, hi) && loaded.find_entrance(lo, hi) == expected_entrance(lo, hi),
          "trial " << trial << ": find_entrance(" << lo << ", " << hi << ")");
  }

  for (int q = 0; q < 60; ++q) {
    const std::string &src = seqs[rng() % seqs.size()];
    std::string P;
    int kind = static_cast<int>(rng() % 4);
    if (kind == 3) {
      P = random_dna(rng, 1 + static_cast<int64_t>(rng() % 30));
    } else {
      size_t a = rng() % src.size();
      P = src.substr(a, 1 + rng() % (src.size() - a));
      if (kind == 1) P = mutate(rng, P, 0.02, 0);
      if (kind == 2) P = random_dna(rng, 1 + static_cast<int64_t>(rng() % 8)) + P;
    }
    SearchResult plain = search_plain(rlbwt, P);
    std::vector<ReportEntry> report;
    SearchResult fast = search_with_shortcuts(rlbwt, sc, P, &report);
    CHECK(plain.lo == fast.lo && plain.hi == fast.hi && plain.matched == fast.matched,
          "trial " << trial << ": search with shortcuts on a pattern of length " << P.size());
    for (const ReportEntry &r : report) {
      int64_t span = ell << r.tier;
      bool ok = r.pattern_pos - span >= 0 && r.pattern_pos + ell <= static_cast<int64_t>(P.size()) &&
                P.substr(static_cast<size_t>(r.pattern_pos - span), static_cast<size_t>(span + ell)) ==
                    R.substr(static_cast<size_t>(r.anchor - span), static_cast<size_t>(span + ell));
      CHECK(ok, "trial " << trial << ": certificate of a report entry");
    }
    SearchResult again = search_with_shortcuts(rlbwt, loaded, P);
    CHECK(again.lo == fast.lo && again.hi == fast.hi && again.matched == fast.matched,
          "trial " << trial << ": search after saving and loading");
  }
}

}  // namespace

int main(int argc, char **argv) {
  int trials = argc > 1 ? std::atoi(argv[1]) : 400;
  std::mt19937_64 rng(12345);
  for (int t = 0; t < trials; ++t) run_trial(rng, t);
  std::cout << checks << " checks, " << failures << " failed\n";
  return failures == 0 ? 0 : 1;
}
