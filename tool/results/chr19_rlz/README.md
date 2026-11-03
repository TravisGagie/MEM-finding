# RLZ of the chr19 two-level parse (16 individuals, both strands) against the CHM13 parse

Desktop, 2026-11-03. The tool was at 7500ed0, unchanged.
- **The parse:** `katka index --parse-oracle c33` (syncmer parse, k = 51, s = 12, the chr19 preset's), built in
  a separate directory with links to `c33.S` and `c33.tbl`: 113 s, peak 1.8 GB.
- **The analysis:** `rlz_parse.cpp`, an analysis program in this directory and not part of the tool. Build line
  at its top; it took 85 s, peak 4.1 GB (`run.log`, `run.sh`).
- **Check:** on the self-test's chromosome-like collection, every count agrees exactly with the brute-force
  `brute_force.py`: H₀, H₁, RLZ, GDC with its pointer runs, and both RLZAP settings.
- **Reference parse:** the phrases starting in S[0 .. 2(|CHM13| + 1)), that is CHM13 and its reverse complement.
  The target is the rest.

## The parse

| | |
|---|---|
| characters n (S$) | 4,065,485,531 |
| phrases m | 220,466,619 (18.4 characters per phrase step; phrases overlap by 51) |
| distinct phrases σ | 8,512,838 (log₂ σ = 23.02) |
| reference parse (CHM13, both strands) | 6,691,777 phrases |
| target | 213,774,842 phrases |
| **H₀** | **22.28 bits per phrase** |
| **H₁** | **0.150 bits per phrase** (8,512,837 contexts, 9,342,836 distinct pairs) |

## Factorizations

| | factors | phrases per factor | literal phrases | pointers |
|---|---|---|---|---|
| greedy RLZ (Kuruppu et al.) | 22,811,810 | **9.37** (37.0 per non-literal factor) | 17,510,796 | 5,301,014 explicit |
| GDC-style (match + one mismatch phrase) | 18,356,612 | 11.65 | one per factor | 3,517,221 runs of equal relative pointers |
| RLZAP-style, J = 4, D = 4, M = 2 | 12,352,342 | 17.31 | 18,644,761 (7,090,969 factors are only a literal) | 1,837,448 changes, 1,611,960 explicit |
| RLZAP-style, J = 8, D = 8, M = 4 | 10,329,469 | 20.70 | 19,370,412 (5,457,307) | 1,793,440 changes, 1,474,034 explicit |

**Rules:**
- **RLZ:** a phrase that does not occur in the reference is a one-phrase literal.
- **GDC and RLZAP:** among equally long matches, the one continuing the previous relative pointer is chosen when
  it is in the suffix-array interval.
- **RLZAP-style:** after each match, look ahead j = 0..J phrases for a match of at least M phrases at relative
  offset rel + δ with |δ| ≤ D (smallest j, then smallest |δ|). The j skipped phrases become literals and δ a
  pointer change. Otherwise, a greedy longest match with an explicit pointer.

## Estimated size of a parse oracle (MB)

| | parse part | dictionary | total |
|---|---|---|---|
| current oracle: ids 661.4 + phrase starts (EF) 191.8 | 853.2 | 630.5 | 1,483.7 |
| greedy RLZ: factor starts 28.5 + offsets 65.6 + literals 52.5 + reference parse 26.0 | **172.6** | 630.5 | **803.1** |
| GDC + relative pointers: factor starts 22.9 + mismatch phrases 55.1 + pointer runs 12.7 + reference 25.9 | 116.6 | 630.5 | 747.2 |
| RLZAP (4,4,2): factor starts 17.0 + literal counts 4.9 + literals 55.9 + pointer differences 7.1 + reference 25.9 | 110.9 | 630.5 | 741.4 |
| RLZAP (8,8,4) | 109.1 | 630.5 | 739.6 |

How each part is sized:
- **Factor starts:** Elias–Fano over positions in S, z(2 + ⌈log₂(n/z)⌉) bits. Text positions are needed for random
  access: a position maps to its factor, then to the reference's phrase starts.
- **Offsets:** the offsets into the reference parse are not monotone, so Elias–Fano does not apply to them
  directly. They are stored at ⌈log₂ R⌉ = 23 bits each.
- **Literals and the reference parse's ids:** ⌈log₂(σ + 1)⌉ = 24 bits each.
- **Reference phrase starts:** Elias–Fano.
- **Dictionary:** the same as the current oracle's: the distinct phrases, plus their offsets.
- **RLZAP's literal counts:** γ codes per factor.
- **RLZAP's pointer differences:** run-length coded. For each non-zero difference, a γ code of the gap since the
  previous one, then sign + γ(|δ|) for a change, or a flag + 23 bits for an explicit pointer.
- **GDC's pointer runs:** each run stores its value plus an Elias–Fano entry for where it starts.

## Reading

- **RLZ shrinks the parse part:** about 4.9 times for greedy RLZ (853 → 173 MB) and 7.7 times for RLZAP
  (853 → 111 MB). The average factor covers 9.4 phrases with greedy RLZ and 17.3 with RLZAP.
- **Literals:** 77% of greedy RLZ's factors are literals, phrases that never occur in CHM13's parse. At 24 bits
  each, they are the largest part after the offsets. RLZAP removes most explicit pointers, not literals.
- **The dictionary (630.5 MB) is now most of the oracle.** It is the same for every variant, so further savings
  must come from it.
- **H₁ = 0.150 bits per phrase** is about 4.1 MB for the whole parse, against about 110 MB for RLZAP's parse part.
  Given the previous phrase, the next one is almost always determined: the property that two-level blocks use.
