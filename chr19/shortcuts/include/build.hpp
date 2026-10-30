// Building the shortcuts from the RLBWT and R (paper, Section 6), and
// pruning them (paper, Section 5).
#pragma once
#include <cstdint>
#include <vector>

#include "rlbwt.hpp"
#include "shortcuts.hpp"

namespace rs {

struct ConstructionStats {
  int64_t anchors = 0;             // multiples p of ell with p + ell <= |R|
  int64_t dropped_unusable = 0;    // window or tier-0 string not all A/C/G/T
  int64_t discarded_colliding = 0; // window equal to another anchor's window
  int64_t backward_steps = 0;
  int64_t deaths = 0;              // rows deleted while tracking deaths
};

// The shortcuts before pruning: one symbol per row, stored plainly.
struct UnprunedShortcuts {
  int64_t ell = 0, max_tier = 0;
  int64_t text_size = 0;  // n, the number of rows of the BWT
  std::vector<uint8_t> reference;
  std::vector<Entrance> entrances;
  std::vector<int8_t> top_tier;
  std::vector<int64_t> target_start, target_first_row;
  std::vector<int64_t> symbol_start;
  std::vector<uint8_t> symbols;
};

// Section 6. Independent of the pruning threshold theta, so it can be
// computed once and pruned with several thresholds.
//
// If `mirrored`, the reference is the reverse of another index's reference
// R (the index of the reversed collection, mems.hpp), and the anchors and
// tiers are aligned with that index's: anchor p here is the boundary
// |R| - p of R, so it is an anchor when |R| - p is a multiple of ell, and
// its tier strings are the reverses of R[|R|-p .. |R|-p + ell*2^t). Then
// both indexes have shortcuts at the same boundaries of R, crossing to the
// left in the index of T and to the right in the other.
UnprunedShortcuts compute_shortcuts(const RLBWT &rlbwt, const std::vector<uint8_t> &reference, int64_t ell,
                                    int64_t max_tier, ConstructionStats &stats, bool mirrored = false);

// Pruning by area (Section 5): for each entrance, keeps tiers 0, 1, ...
// while (rows surviving tier t) * ell * 2^t >= theta, caps the symbols at
// the first tier that falls short, and drops the anchor if tier 0 does.
// (The second rule, on characters other than A/C/G/T, is applied by
// compute_shortcuts.)
Shortcuts prune(const UnprunedShortcuts &unpruned, int64_t theta, int64_t &dropped);

}  // namespace rs
