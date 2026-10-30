// Finding maximal exact matches (MEMs) with two indexes: the RLBWT of T,
// searched backward, and the RLBWT of the reversed collection
// (text.hpp), in which backward search for the reverse of a string reads
// the string left to right, so that matches can be extended forward.
// Either index may have shortcuts; the MEMs found are the same, since a
// search with shortcuts finds the same longest matching suffix and
// interval as plain backward search (Theorem 2).
//
// A MEM of a pattern P is a substring P[s .. e) that occurs in T while
// P[s-1 .. e) (if s > 0) and P[s .. e+1) (if e < |P|) do not. MEMs are
// never nested, so ordering them by start orders them by end.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "rlbwt.hpp"
#include "search.hpp"
#include "shortcuts.hpp"

namespace rs {

struct Mem {
  int64_t start, end;  // P[start .. end)
  int64_t lo, hi;      // its interval in the BWT of T, e.g. for a tag array
  bool operator==(const Mem &o) const { return start == o.start && end == o.end && lo == o.lo && hi == o.hi; }
};

struct TwoIndexes {
  const RLBWT *text = nullptr;               // the RLBWT of T
  const RLBWT *reversed = nullptr;           // the RLBWT of the reversed collection
  const Shortcuts *text_sc = nullptr;        // shortcuts for each, or null for plain search
  const Shortcuts *reversed_sc = nullptr;
};

// Totals over the searches a MEM finder makes.
struct MemCounters {
  int64_t searches = 0, backward_steps = 0, shortcuts = 0, chars_crossed = 0;
  void add(const SearchResult &r) {
    ++searches;
    backward_steps += r.backward_steps;
    shortcuts += r.shortcuts;
    chars_crossed += r.chars_crossed;
  }
};

// The MEMs of P of length at least min_length, in order, by forward-
// backward search: from a position x, a forward search finds the longest
// P[x .. e) that occurs; a backward search from e finds the MEM P[s .. e)
// and its interval; and a backward search from e + 1 finds the leftmost
// position after s where a match reaching past e can start.
std::vector<Mem> mems_forward_backward(const TwoIndexes &idx, const std::string &P, int64_t min_length,
                                       MemCounters *counters = nullptr);

// The same MEMs, by Boyer-Moore-Li: from a position x, a backward search
// for P[x .. x + min_length) either fails at some position k, so that no
// MEM of length at least min_length starts in [x, k] and the next
// position tried is k + 1, or succeeds and continues as forward-backward.
std::vector<Mem> mems_bml(const TwoIndexes &idx, const std::string &P, int64_t min_length,
                          MemCounters *counters = nullptr);

}  // namespace rs
