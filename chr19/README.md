# More RLBWT Tricks

Code for the paper *More RLBWT Tricks* (Travis Gagie). It finds maximal exact
matches (MEMs) with forward–backward search or Boyer–Moore–Li over the
run-length compressed BWTs (RLBWTs) of a collection and of the reversed
collection, in three ways, all returning the same MEMs and BWT intervals as plain
backward search:

- **shortcuts** across windows of a reference sequence (`shortcuts/`);
- **two-level indexing** over a parse of the collection at closed syncmers
  (k = 51, s = 12), bouncing between the RLBWT and an index of the parse;
- **block jumps**, a third index over blocks of phrases, each a phrase followed
  by its continuations, that lets a search cross a whole block in one step
  (`twolevel/`).

## Layout

- `shortcuts/`: the RLBWT, shortcuts and their tools (`build_index`, `query_index`,
  `query_mems`), with libsais in `shortcuts/third_party/`.
- `twolevel/`: the two-level index (`build_twolevel`, `query_twolevel`) and the
  scripts for the paper's experiments.

## Building and testing

Requires CMake 3.10 or later, a C++17 compiler and Python 3. Everything is
compiled with `-march=native`.

    bash twolevel/selftest.sh

This builds everything and checks, on a small synthetic collection, that every
method finds exactly the MEMs and intervals of plain backward search.

## Usage

The collection is a file with one sequence per line. Build the RLBWTs of the
collection and of the reversed collection, then the parse indexes:

    B=twolevel/build
    $B/rs_build/build_index --text coll.txt --reference-length LEN --out idx  --ell 25 --theta 0
    $B/rs_build/build_index --text coll.txt --reference-length LEN --out idxr --ell 25 --theta 0 --reversed
    $B/build_twolevel --text coll.txt --rlbwt idx/rlbwt.bin --reversed-rlbwt idxr/rlbwt.bin --out lv

The reference is the first sequence, of length `LEN`. Then find the MEMs of the
patterns in a file, one per line, over {A, C, G, T}:

    $B/query_twolevel --index idx --reversed-index idxr --level lv.text.lvl --reversed-level lv.rev.lvl \
      --patterns patterns.txt --configs plain,2L,2LJ --algorithm bml --min-length 1000

`2L` is two-level indexing, `2LJ` two-level indexing with block jumps.
`--algorithm fb` selects forward–backward search, and `--dump DIR` writes the
MEMs of the first configuration. `build_twolevel --continuation frequent` uses the
most frequent successor of each phrase as its continuation instead of the
smallest, and `build_twolevel --text coll.txt --block-stats` compares the number
of blocks for four choices.

## The experiments

`twolevel/run_final.sh` runs the paper's comparison on chromosome 19 from the
Human Pangenome Reference Consortium (release 2), with collections of 16, 32 and
64 individuals: plain search, shortcuts, two-level indexing with and without
block jumps, KeBaB, ropebwt3 and Movi 2. `twolevel/run_continuation.sh` compares
the smallest and the most frequent successor as continuations, and
`twolevel/run_movi_final.sh` reruns only the Movi 2 part. The paths to the data
are set at the top of each script.

## License

GPL v3 or later (see `LICENSE`); libsais is under the Apache License 2.0.
