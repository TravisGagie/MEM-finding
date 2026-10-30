// Counts patterns with and without shortcuts, checks that the two
// searches agree, and times them.
//
//   query_index --rlbwt FILE --shortcuts FILE --pattern STRING
//       prints both results and the report;
//   query_index --rlbwt FILE --shortcuts FILE --patterns FILE [--repeat K] [--csv FILE]
//       one pattern per line; times each search K times (default 5),
//       alternating, keeps the median, and collects the report in a separate untimed
//       run. Exits with status 2 if any search with shortcuts disagrees
//       with plain backward search.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "search.hpp"

using namespace rs;
using Clock = std::chrono::steady_clock;

namespace {

int usage() {
  std::cerr << "usage: query_index --rlbwt FILE --shortcuts FILE (--pattern P | --patterns FILE [--repeat K] [--csv FILE])\n";
  return 1;
}

bool same(const SearchResult &a, const SearchResult &b) {
  return a.lo == b.lo && a.hi == b.hi && a.matched == b.matched;
}

// Runs search() once and returns its time in nanoseconds.
template <typename F>
double time_ns(F search, SearchResult &result) {
  auto t0 = Clock::now();
  result = search();
  return std::chrono::duration<double, std::nano>(Clock::now() - t0).count();
}

}  // namespace

int main(int argc, char **argv) {
  std::string rlbwt_path, shortcuts_path, pattern, patterns_path, csv_path;
  bool single = false;
  int64_t repeat = 5;
  for (int i = 1; i + 1 < argc; i += 2) {
    std::string flag = argv[i], value = argv[i + 1];
    if (flag == "--rlbwt") rlbwt_path = value;
    else if (flag == "--shortcuts") shortcuts_path = value;
    else if (flag == "--pattern") { pattern = value; single = true; }
    else if (flag == "--patterns") patterns_path = value;
    else if (flag == "--repeat") repeat = std::max<int64_t>(1, std::stoll(value));
    else if (flag == "--csv") csv_path = value;
    else return usage();
  }
  if (argc % 2 == 0 || rlbwt_path.empty() || shortcuts_path.empty() || (!single && patterns_path.empty()))
    return usage();

  try {
    const RLBWT rlbwt = RLBWT::load(rlbwt_path);
    const Shortcuts sc = Shortcuts::load(shortcuts_path);

    if (single) {
      std::vector<ReportEntry> report;
      SearchResult plain = search_plain(rlbwt, pattern);
      SearchResult fast = search_with_shortcuts(rlbwt, sc, pattern, &report);
      std::cout << "plain:          interval [" << plain.lo << ", " << plain.hi << "], matched " << plain.matched
                << ", " << plain.backward_steps << " backward steps\n"
                << "with shortcuts: interval [" << fast.lo << ", " << fast.hi << "], matched " << fast.matched
                << ", " << fast.backward_steps << " backward steps, " << fast.shortcuts << " shortcuts\n"
                << "report (pattern position, anchor, tier):\n";
      for (const ReportEntry &r : report) std::cout << "  " << r.pattern_pos << ' ' << r.anchor << ' ' << r.tier << '\n';
      return same(plain, fast) ? 0 : 2;
    }

    std::ifstream in(patterns_path);
    if (!in) throw std::runtime_error("cannot open " + patterns_path);
    std::ofstream csv;
    if (!csv_path.empty()) {
      csv.open(csv_path);
      csv << "pattern,length,found,matched,plain_ns,shortcut_ns,plain_steps,shortcut_steps,shortcuts,"
             "chars_crossed,distinct_anchors\n";
    }
    int64_t count = 0, mismatches = 0;
    double plain_total = 0, shortcut_total = 0;
    for (std::string P; std::getline(in, P);) {
      if (!P.empty() && P.back() == '\r') P.pop_back();
      if (P.empty()) continue;
      SearchResult plain, fast;
      std::vector<double> plain_times, fast_times;
      // Interleave the two searches, as in the experiments.
      for (int64_t k = 0; k < repeat; ++k) {
        double a = time_ns([&] { return search_plain(rlbwt, P); }, plain);
        double b = time_ns([&] { return search_with_shortcuts(rlbwt, sc, P); }, fast);
        plain_times.push_back(a);
        fast_times.push_back(b);
      }
      auto median = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        const size_t k = v.size() / 2;
        return v.size() % 2 ? v[k] : (v[k - 1] + v[k]) / 2;
      };
      const double plain_ns = median(plain_times), fast_ns = median(fast_times);
      std::vector<ReportEntry> report;
      search_with_shortcuts(rlbwt, sc, P, &report);
      std::vector<int64_t> anchors;
      for (const ReportEntry &r : report) anchors.push_back(r.anchor);
      std::sort(anchors.begin(), anchors.end());
      int64_t distinct = std::unique(anchors.begin(), anchors.end()) - anchors.begin();

      if (!same(plain, fast)) {
        ++mismatches;
        std::cerr << "MISMATCH on pattern " << count << '\n';
      }
      plain_total += plain_ns;
      shortcut_total += fast_ns;
      if (csv.is_open()) {
        csv << count << ',' << P.size() << ',' << plain.found() << ',' << plain.matched << ',' << plain_ns << ','
            << fast_ns << ',' << plain.backward_steps << ',' << fast.backward_steps << ',' << fast.shortcuts << ','
            << fast.chars_crossed << ',' << distinct << '\n';
      }
      ++count;
    }
    std::cout << count << " patterns, " << mismatches << " mismatches; total time " << plain_total / 1e6
              << " ms plain, " << shortcut_total / 1e6 << " ms with shortcuts (speedup "
              << (shortcut_total > 0 ? plain_total / shortcut_total : 0) << ")\n";
    return mismatches ? 2 : 0;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}
