#include "mems.hpp"

#include <algorithm>
#include <stdexcept>

namespace rs {

namespace {

class MemFinder {
 public:
  MemFinder(const TwoIndexes &idx, const std::string &P, MemCounters *counters)
      : idx_(idx), P_(P), Q_(P.rbegin(), P.rend()), m_(static_cast<int64_t>(P.size())), counters_(counters) {
    if (!idx.text || !idx.reversed) throw std::invalid_argument("both indexes are needed");
  }

  int64_t length() const { return m_; }

  // The largest e such that P[x .. e) occurs: backward search in the
  // reversed collection for reverse(P[x .. |P|)) = Q[0 .. |P| - x), which
  // reads P[x], P[x+1], ...
  int64_t forward(int64_t x) {
    SearchResult r = search_range(*idx_.reversed, idx_.reversed_sc, Q_, 0, m_ - x);
    if (counters_) counters_->add(r);
    return x + r.matched;
  }

  // The longest suffix of P[begin .. end) that occurs, and its interval.
  SearchResult backward(int64_t begin, int64_t end) {
    SearchResult r = search_range(*idx_.text, idx_.text_sc, P_, begin, end);
    if (counters_) counters_->add(r);
    return r;
  }

  // P[x .. e) occurs and e is as large as possible, with e > x. Appends
  // the MEM ending at e if it is long enough, and returns the next
  // position to try: every MEM not yet found starts there or later.
  int64_t report_and_advance(int64_t x, int64_t e, int64_t min_length, std::vector<Mem> &out) {
    // P[s .. e) is the longest suffix of P[0 .. e) that occurs; s <= x.
    // P[s .. e+1) does not occur, since P[x .. e+1) does not, so P[s .. e)
    // is a MEM.
    SearchResult r = backward(0, e);
    const int64_t s = e - r.matched;
    if (e - s >= min_length) out.push_back({s, e, r.lo, r.hi});
    if (e == m_) return m_;
    // A later MEM P[a .. b) has b > e, so P[a .. e+1) occurs and a >= s1,
    // where P[s1 .. e+1) is the longest suffix of P[0 .. e+1) that occurs.
    // It also has a > x, since P[x .. b) does not occur.
    const int64_t s1 = e + 1 - backward(0, e + 1).matched;
    return std::max(s1, x + 1);
  }

 private:
  const TwoIndexes &idx_;
  const std::string &P_;
  const std::string Q_;  // P reversed
  const int64_t m_;
  MemCounters *counters_;
};

}  // namespace

std::vector<Mem> mems_forward_backward(const TwoIndexes &idx, const std::string &P, int64_t min_length,
                                       MemCounters *counters) {
  MemFinder f(idx, P, counters);
  std::vector<Mem> out;
  for (int64_t x = 0; x < f.length();) {
    const int64_t e = f.forward(x);
    if (e == x) {  // P[x] does not occur in T, so no MEM contains it
      ++x;
      continue;
    }
    x = f.report_and_advance(x, e, min_length, out);
  }
  return out;
}

std::vector<Mem> mems_bml(const TwoIndexes &idx, const std::string &P, int64_t min_length,
                          MemCounters *counters) {
  if (min_length < 1) throw std::invalid_argument("min_length must be positive");
  MemFinder f(idx, P, counters);
  std::vector<Mem> out;
  for (int64_t x = 0; x + min_length <= f.length();) {
    const int64_t window_end = x + min_length;
    const SearchResult r = f.backward(x, window_end);
    if (r.matched < min_length) {
      // P[k .. window_end) does not occur, for k = window_end - matched - 1,
      // and every MEM of length at least min_length starting in [x, k]
      // would contain it.
      x = window_end - r.matched;
      continue;
    }
    x = f.report_and_advance(x, f.forward(x), min_length, out);
  }
  return out;
}

}  // namespace rs
