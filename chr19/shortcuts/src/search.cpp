#include "search.hpp"

#include <algorithm>

namespace rs {

SearchResult search_plain(const RLBWT &rlbwt, const std::string &P) {
  SearchResult r;
  r.pattern_length = static_cast<int64_t>(P.size());
  r.lo = 0;
  r.hi = rlbwt.size() - 1;  // the interval of the empty suffix
  for (int64_t j = r.pattern_length - 1; j >= 0; --j) {
    int64_t lo = r.lo, hi = r.hi;
    rlbwt.backward_step(static_cast<uint8_t>(P[static_cast<size_t>(j)]), lo, hi);
    ++r.backward_steps;
    if (hi < lo) break;  // P[j .. |P|) does not occur
    r.lo = lo;
    r.hi = hi;
    ++r.matched;
  }
  return r;
}

SearchResult search_with_shortcuts(const RLBWT &rlbwt, const Shortcuts &sc, const std::string &P,
                                   std::vector<ReportEntry> *report) {
  return search_range(rlbwt, &sc, P, 0, static_cast<int64_t>(P.size()), report);
}

SearchResult search_range(const RLBWT &rlbwt, const Shortcuts *sc, const std::string &P, int64_t begin, int64_t end,
                          std::vector<ReportEntry> *report) {
  SearchResult r;
  r.pattern_length = end - begin;
  const int64_t ell = sc ? sc->ell : 0;
  int64_t lo = 0, hi = rlbwt.size() - 1;  // interval of the current match X = P[j .. end)
  int64_t j = end;

  while (j > begin) {
    // If |X| >= ell and I(X) lies in the entrance of anchor p, then beta_p
    // is a prefix of X (Lemma 1).
    int64_t first = 0;
    int64_t e = sc && r.matched >= ell ? sc->find_entrance(lo, hi, &first) : -1;
    if (e >= 0) {
      // Compare P[j-1], P[j-2], ... with R[p-1], R[p-2], ..., without
      // going past P[begin].
      const int64_t top = sc->top_tier(e);
      const int64_t longest = ell << top;
      const int64_t p = sc->anchor(e);  // alpha_{p,top} = R[p - longest .. p)
      int64_t agree = 0;
      while (agree < std::min(longest, j - begin) &&
             static_cast<uint8_t>(P[static_cast<size_t>(j - 1 - agree)]) == sc->reference_char(p - 1 - agree)) {
        ++agree;
      }
      // Take the largest tier whose string P matches and whose interval is
      // not empty.
      bool jumped = false;
      for (int64_t t = top; t >= 0 && !jumped; --t) {
        const int64_t len = ell << t;
        int64_t new_lo, new_hi;
        if (len > agree || !sc->jump(e, first, t, lo, hi, new_lo, new_hi)) continue;
        if (report) report->push_back({j, p, t});
        lo = new_lo;
        hi = new_hi;
        j -= len;
        r.matched += len;
        ++r.shortcuts;
        r.chars_crossed += len;
        jumped = true;
      }
      if (jumped) continue;
    }

    int64_t new_lo = lo, new_hi = hi;
    rlbwt.backward_step(static_cast<uint8_t>(P[static_cast<size_t>(j - 1)]), new_lo, new_hi);
    ++r.backward_steps;
    if (new_hi < new_lo) break;  // P[j-1 .. end) does not occur
    lo = new_lo;
    hi = new_hi;
    ++r.matched;
    --j;
  }
  r.lo = lo;
  r.hi = hi;
  return r;
}

}  // namespace rs
