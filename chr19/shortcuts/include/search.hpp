// Counting with and without shortcuts (paper, Section 4, "Searching").
// Both find the longest suffix of P that occurs in T and its interval;
// P occurs in T if that suffix is all of P.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "rlbwt.hpp"
#include "shortcuts.hpp"

namespace rs {

struct SearchResult {
  int64_t lo = 0, hi = -1;     // interval of the longest suffix of P that occurs in T
  int64_t matched = 0;         // its length
  int64_t pattern_length = 0;
  int64_t backward_steps = 0;
  int64_t shortcuts = 0;       // shortcuts taken
  int64_t chars_crossed = 0;   // characters of P crossed by shortcuts
  bool found() const { return matched == pattern_length; }
};

// One entry of the report: a tier-`tier` shortcut taken when the current
// match was P[pattern_pos .. |P|), certifying
//   P[pattern_pos - ell*2^tier .. pattern_pos + ell) = R[anchor - ell*2^tier .. anchor + ell).
struct ReportEntry {
  int64_t pattern_pos;
  int64_t anchor;
  int64_t tier;
};

SearchResult search_plain(const RLBWT &rlbwt, const std::string &P);

// If `report` is not null, the shortcuts taken are appended to it.
SearchResult search_with_shortcuts(const RLBWT &rlbwt, const Shortcuts &sc, const std::string &P,
                                   std::vector<ReportEntry> *report = nullptr);

// Backward search for P[begin .. end): finds its longest suffix that occurs
// in T, and the interval of that suffix, with shortcuts if `sc` is not
// null and with plain backward steps otherwise. pattern_length is
// end - begin, and report entries give positions in P.
SearchResult search_range(const RLBWT &rlbwt, const Shortcuts *sc, const std::string &P, int64_t begin, int64_t end,
                          std::vector<ReportEntry> *report = nullptr);

}  // namespace rs
