# Reportable shortcuts with RLBWTs

Code for the paper *Reportable Shortcuts with RLBWTs*. It counts the
occurrences of a pattern with backward search over a run-length
compressed BWT (RLBWT), taking shortcuts across stretches of the pattern
that match a reference, and reports the shortcuts it took.

## Where the paper is in the code

| Paper | Code |
|---|---|
| Section 2: the RLBWT and backward steps | `rlbwt.hpp`, `rlbwt.cpp` |
| The text T, with a separator after each sequence, and R as a prefix of T (Sections 4 and 7) | `text.hpp` |
| Lemma 1 (detection by interval containment) | `Shortcuts::find_entrance` in `shortcuts.hpp`, called from `search.cpp` |
| Lemma 3 (crossing alpha with rank queries) | `Shortcuts::jump` in `shortcuts.cpp` |
| Section 4: anchors, windows, tier strings, collisions, top tiers | `search_windows`, `discard_collisions`, `top_tier_of` in `build.cpp` |
| Section 4, "Searching" | `search_with_shortcuts` in `search.cpp` |
| Section 5: the data structure | `Shortcuts` in `shortcuts.hpp` |
| Section 5: pruning by area | `prune` in `build.cpp` |
| Section 5: dropping and capping for characters other than A/C/G/T | `search_windows`, `top_tier_of` in `build.cpp` |
| Section 6: construction by tracking deaths | `extend_anchor` in `build.cpp` |
| MEM-finding with two indexes (forward-backward, Boyer-Moore-Li) | `mems.hpp`, `mems.cpp`; the reversed collection in `text.hpp` |
| Section 7: experiments | `tools/`, `scripts/` |

`shortcuts.hpp` lists the paper's notation and the names it has in the
code. The other files, in the order they build on each other:

- `io.hpp`: reading and writing binary files;
- `bitvector.hpp`: bitvectors with rank and select, plain or Elias-Fano;
- `wavelet_tree.hpp`: the wavelet trees over run characters and over
  symbols.

Some details differ from the paper's description:

- R and the per-entrance arrays are stored plainly, not packed (as
  Section 7 says).
- The RLBWT is built either in memory with libsais (as for cere in the
  paper) or from a BWT computed by Big-BWT, which uses prefix-free
  parsing (as for the bacterial collection).

The RLBWT is the usual one: Elias-Fano bitvectors marking where the runs
start, in BWT order and with the runs sorted by character, and a wavelet
tree with plain bitvectors over the run characters (`rlbwt.hpp`). The
plain bitvectors, which also hold the high parts of the Elias-Fano ones,
store rank counts per 512-bit block and per 64-bit word, as in Vigna's
rank9, and sample every 512th 1 and 0 for select (`bitvector.hpp`).

## Building

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j

The build uses `-march=native`, for hardware popcount and BMI2; add
`-DNATIVE=OFF` to the first command for a portable build.
    ./build/test_bitvector
    ./build/test_shortcuts
    ./build/test_mems

The tests build random collections and check every entrance, symbol and
f_{p,t} against definitions computed from the suffix array. They also
check that search with shortcuts gives the same result as plain backward
search, that every report entry's certificate holds, and that an RLBWT
read from a BWT file equals one built from the text. `test_mems` checks
that both MEM finders, with and without shortcuts on either index, find
exactly the MEMs of the definition, with their intervals.

libsais (Apache 2.0) is included under `third_party/`.

## Using it

A collection file has one sequence per line, or is in FASTA format. The
reference R is its first N characters, counting one separator after
each sequence.

    build_index --text COLLECTION --reference-length N --out DIR \
                [--bwt BWT] [--ell 25,50] [--theta 0,400] [--max-tier 3]

This writes `DIR/rlbwt.bin` and `DIR/ell<L>_theta<T>.shortcuts` for each
pair, and appends one line per pair to `DIR/build.csv`. Without
`--bwt`, the BWT is computed in memory with libsais, which takes about 9
bytes per character. For larger collections, compute the BWT with
[Big-BWT](https://gitlab.com/manzai/Big-BWT) (`bigbwt COLLECTION`) and
pass `--bwt COLLECTION.bwt`; COLLECTION must then have one sequence per
line, with no blank lines and a final newline, as `prepare_collection.py`
writes. `--rlbwt DIR/rlbwt.bin` reuses an RLBWT built earlier.

    query_index --rlbwt DIR/rlbwt.bin --shortcuts DIR/ell25_theta0.shortcuts --pattern ACGT...

prints the result of plain search, the result with shortcuts, and the
report: one line per shortcut, with the pattern position j, the anchor p
and the tier t, meaning that P[j - ell*2^t .. j + ell) = R[p - ell*2^t .. p + ell).

    query_index --rlbwt ... --shortcuts ... --patterns FILE [--repeat 3] [--csv OUT]

times both searches on each pattern of FILE, as in the experiments.

For MEM-finding, build a second index of the reversed collection (the
lines of COLLECTION in reverse order, each reversed), with the reverse
of R as reference:

    build_index --text COLLECTION --reference-length N --out DIR_REV --reversed \
                [--bwt BWT_OF_REVERSED] [--ell 25,50] [--theta 0]

where `scripts/reverse_collection.py COLLECTION REVERSED` writes the
reversed collection for Big-BWT. Then

    query_mems --index DIR --reversed-index DIR_REV --shortcuts ell25_theta0 \
               --patterns FILE [--algorithm fb,bml] [--min-length 1,50] [--repeat 3] [--csv OUT]

finds the MEMs of each pattern with forward-backward (`fb`) or
Boyer-Moore-Li (`bml`), with and without shortcuts on both indexes,
checks that they agree, and times them. `scripts/run_mem_experiment.sh`
runs this for a grid of configurations and summarizes the results.

## Reproducing the experiments

    scripts/run_experiment.sh COLLECTION N OUTDIR [BWT]

builds the index for a grid of (ell, theta), generates the patterns,
runs the queries, and writes `OUTDIR/results.csv` and
`OUTDIR/summary.txt`. The grid is set with the environment variables
`ELLS`, `THETAS`, `PERCENTS` and `REPEAT`; see the top of the script.
For the paper, the queries were pinned to one performance core with
`PIN=<cpu>` (rerunning only them with `QUERY_ONLY=1`).

- **cere** (Pizza&Chili): `run_experiment.sh cere 12000000 out` with the
  defaults, given about 4 GB of memory.
- **Bacterial collection**: from the seven species files of SPUMONI 2's
  reference,

      scripts/prepare_collection.py zymo.txt --ref-min-len 2000000 \
          --max-bases-per-file 2000000000 --rest <the seven .fa files>

  prints N. Then `bigbwt zymo.txt`, and

      ELLS=16,25,50,100 THETAS=0,400,800 PERCENTS=0,0.1,1 scripts/run_experiment.sh zymo.txt N out zymo.txt.bwt

- **MEMs of chromosome 19** (HPRC release 2): `scripts/prepare_hprc_chr19.py
  DIR` downloads chromosome 19 of CHM13, as the reference, and of the
  haplotypes of the first 16 suitable individuals, and two query
  haplotypes, fetching only those bases from the bgzipped assemblies on
  S3. Then build both indexes, from Big-BWT runs on `DIR/collection.txt`
  and on its reversal by `scripts/reverse_collection.py`:

      build_index --text collection.txt --bwt collection.txt.bwt --reference-length $(cat ref_length.txt) \
                  --out index --ell 16,25,50 --theta 0 --max-tier 6
      build_index ... --bwt collection_rev.txt.bwt --out index_rev ... --reversed

  and run `query_mems --index index --reversed-index index_rev --shortcuts
  ell25_theta0 --patterns queries.txt --algorithm bml --min-length 1000
  --csv out.csv`; `scripts/summarize_mems.py` summarizes the CSVs.
