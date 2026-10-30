// Finds the MEMs of patterns with two indexes (mems.hpp), without
// shortcuts and with each of several sets of shortcuts, checks that all
// find the same MEMs, and times them.
//
//   query_mems --index DIR --reversed-index DIR --shortcuts NAME[,NAME...] --patterns FILE
//              [--algorithm fb,bml] [--min-length 1,50] [--repeat K] [--csv FILE]
//
// DIR/rlbwt.bin and DIR/NAME.shortcuts are the index of T, as build_index
// writes it; the reversed index, with shortcuts of the same names, is
// built with build_index --reversed. For each algorithm (fb:
// forward-backward; bml: Boyer-Moore-Li), minimum MEM length and pattern,
// finds the MEMs with plain backward search and then with each set of
// shortcuts in turn, K times (default 5), and reports the median time of
// each, with the fastest and slowest. Plain search thus runs once per
// round, however many sets of shortcuts are compared. Exits with status 2
// if any results disagree.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "mems.hpp"

using namespace rs;
using Clock = std::chrono::steady_clock;

namespace {

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const size_t k = v.size() / 2;
  return v.size() % 2 ? v[k] : (v[k - 1] + v[k]) / 2;
}

int usage() {
  std::cerr << "usage: query_mems --index DIR --reversed-index DIR --shortcuts NAME[,NAME...] --patterns FILE\n"
               "                  [--algorithm fb,bml] [--min-length 1,50] [--repeat K] [--csv FILE]\n";
  return 1;
}

std::vector<std::string> split(const std::string &s) {
  std::vector<std::string> v;
  std::stringstream ss(s);
  for (std::string item; std::getline(ss, item, ',');) v.push_back(item);
  return v;
}

}  // namespace

int main(int argc, char **argv) {
  std::string index_dir, reversed_dir, patterns_path, csv_path;
  std::vector<std::string> names, algorithms = {"fb", "bml"};
  std::vector<int64_t> min_lengths = {1};
  int64_t repeat = 5;
  for (int i = 1; i + 1 < argc; i += 2) {
    std::string flag = argv[i], value = argv[i + 1];
    if (flag == "--index") index_dir = value;
    else if (flag == "--reversed-index") reversed_dir = value;
    else if (flag == "--shortcuts") names = split(value);
    else if (flag == "--patterns") patterns_path = value;
    else if (flag == "--algorithm") algorithms = split(value);
    else if (flag == "--min-length") {
      min_lengths.clear();
      for (const std::string &x : split(value)) min_lengths.push_back(std::stoll(x));
    } else if (flag == "--repeat") repeat = std::max<int64_t>(1, std::stoll(value));
    else if (flag == "--csv") csv_path = value;
    else return usage();
  }
  if (argc % 2 == 0 || index_dir.empty() || reversed_dir.empty() || names.empty() || patterns_path.empty())
    return usage();
  for (const std::string &a : algorithms)
    if (a != "fb" && a != "bml") return usage();

  try {
    const RLBWT text = RLBWT::load(index_dir + "/rlbwt.bin");
    const RLBWT reversed = RLBWT::load(reversed_dir + "/rlbwt.bin");
    std::vector<Shortcuts> sc, sc_rev;
    for (const std::string &name : names) {
      sc.push_back(Shortcuts::load(index_dir + "/" + name + ".shortcuts"));
      sc_rev.push_back(Shortcuts::load(reversed_dir + "/" + name + ".shortcuts"));
    }
    const TwoIndexes plain{&text, &reversed, nullptr, nullptr};
    std::vector<TwoIndexes> fast;
    for (size_t c = 0; c < names.size(); ++c) fast.push_back({&text, &reversed, &sc[c], &sc_rev[c]});

    std::vector<std::string> patterns;
    {
      std::ifstream in(patterns_path);
      if (!in) throw std::runtime_error("cannot open " + patterns_path);
      for (std::string P; std::getline(in, P);) {
        if (!P.empty() && P.back() == '\r') P.pop_back();
        if (!P.empty()) patterns.push_back(P);
      }
    }
    std::ofstream csv;
    if (!csv_path.empty()) {
      csv.open(csv_path);
      csv << "config,algorithm,min_length,pattern,length,mems,mem_chars,plain_ns,shortcut_ns,plain_steps,"
             "shortcut_steps,searches,shortcuts,chars_crossed,repeats,plain_min_ns,plain_max_ns,shortcut_min_ns,"
             "shortcut_max_ns\n";
    }

    int64_t mismatches = 0;
    for (const std::string &algorithm : algorithms) {
      auto find = [&](const TwoIndexes &idx, const std::string &P, int64_t L, MemCounters *c) {
        return algorithm == "fb" ? mems_forward_backward(idx, P, L, c) : mems_bml(idx, P, L, c);
      };
      for (int64_t L : min_lengths) {
        double plain_total = 0;
        std::vector<double> fast_total(names.size(), 0);
        for (size_t q = 0; q < patterns.size(); ++q) {
          const std::string &P = patterns[q];
          std::vector<Mem> a;
          MemCounters ca;
          std::vector<std::vector<Mem>> b(names.size());
          std::vector<MemCounters> cb(names.size());
          // Rounds of plain search followed by each set of shortcuts.
          std::vector<double> plain_times;
          std::vector<std::vector<double>> fast_times(names.size());
          for (int64_t k = 0; k < repeat; ++k) {
            ca = MemCounters();
            auto t0 = Clock::now();
            a = find(plain, P, L, &ca);
            plain_times.push_back(std::chrono::duration<double, std::nano>(Clock::now() - t0).count());
            for (size_t c = 0; c < names.size(); ++c) {
              cb[c] = MemCounters();
              auto t1 = Clock::now();
              b[c] = find(fast[c], P, L, &cb[c]);
              fast_times[c].push_back(std::chrono::duration<double, std::nano>(Clock::now() - t1).count());
              if (b[c] != a) {
                ++mismatches;
                std::cerr << "MISMATCH: " << names[c] << ", " << algorithm << ", L = " << L << ", pattern " << q
                          << '\n';
              }
            }
          }
          const double plain_ns = median(plain_times);
          std::vector<double> fast_ns(names.size());
          for (size_t c = 0; c < names.size(); ++c) fast_ns[c] = median(fast_times[c]);
          int64_t mem_chars = 0;
          for (const Mem &x : a) mem_chars += x.end - x.start;
          plain_total += plain_ns;
          for (size_t c = 0; c < names.size(); ++c) {
            fast_total[c] += fast_ns[c];
            if (csv.is_open()) {
              csv << names[c] << ',' << algorithm << ',' << L << ',' << q << ',' << P.size() << ',' << a.size() << ','
                  << mem_chars << ',' << plain_ns << ',' << fast_ns[c] << ',' << ca.backward_steps << ','
                  << cb[c].backward_steps << ',' << cb[c].searches << ',' << cb[c].shortcuts << ','
                  << cb[c].chars_crossed << ',' << repeat << ','
                  << *std::min_element(plain_times.begin(), plain_times.end()) << ','
                  << *std::max_element(plain_times.begin(), plain_times.end()) << ','
                  << *std::min_element(fast_times[c].begin(), fast_times[c].end()) << ','
                  << *std::max_element(fast_times[c].begin(), fast_times[c].end()) << '\n';
            }
          }
        }
        for (size_t c = 0; c < names.size(); ++c) {
          std::cout << names[c] << ", " << algorithm << ", L = " << L << ": " << patterns.size()
                    << " patterns; total of median times " << plain_total / 1e6 << " ms plain, " << fast_total[c] / 1e6
                    << " ms with shortcuts (speedup " << plain_total / fast_total[c] << ")\n";
        }
      }
    }
    if (mismatches > 0) {
      std::cerr << mismatches << " mismatches\n";
      return 2;
    }
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
