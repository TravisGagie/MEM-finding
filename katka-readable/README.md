# KATKA, readable version

This is a cleaned-up copy of the code behind *Taxonomic Classification with Tag Arrays*. It covers only what the paper's main pipeline uses:

- building the text, the BWT structures and the tag array;
- finding long MEMs with Boyer–Moore–Li;
- answering each MEM from the tag array.

It was written for someone who has read the paper. Every file starts with a comment saying which part of the paper it implements.

The programs are named `katka-*` instead of the original's `rz-*` (the rz-index of Appendix C is no longer used). Their options and file formats are the same, except as listed under *Differences from the original*, so the index files of the two versions can be mixed. On the test collection (40 genomes, 4 species, 3,000 reads), every build output is byte-identical to the original's. The classification outputs are identical for every option combination below (see *Checking against the original*).

## Map from the paper to the code

| Paper | File | What is there |
|---|---|---|
| §5 "The index", App. A "Reference text" | `src/prep.cpp` (`katka-prep`) | the text S = G₁ X rc(G₁) X … with both strands; `-d` datasets, `-s` genomes split by sample ID, `-U` homopolymer compression |
| §5 "The index" | `src/rlbwt.hpp` | the RLBWT (run starts, run heads, per-character run lengths), LF and the backward step |
| §5 "The index" | `src/rlcsa.hpp`, `src/elias_fano.hpp` | the RLCSA: per-character run lists, a backward step as a predecessor search (explicit lists with jump tables, or Elias–Fano) |
| §5 "The index", App. A "Lookup tables" | `src/search.hpp` | the BWT interval, one backward step on either structure, the t-mer lookup table |
| §5 "Boyer–Moore–Li", Fig. 1, App. A "BML without digests" | `src/bml.hpp` | BML on one segment of a read (`bml_segment`) |
| §5 "Long, noisy reads" | `src/bml.hpp` (`distinct_3mers`), `src/classify.cpp` (`-U`) | the low-complexity filter τ (`-D`) and homopolymer compression of reads |
| §5 "Longest MEMs" | `src/bml.hpp` (`find_mems`, `segments`), `src/kmer_filter.hpp` | adaptive L (`-X`); KeBaB pseudo-MEMs searched longest first (`-K`) |
| §6 "Tag arrays", "Counting" | `src/tag_array.hpp` | the run-length tag array; `count` (proportional credit); `list` with Muthukrishnan/Sadakane (even credit, range); sampled tags via LF (s > 1) |
| §6 "Construction", App. A "Indexes" | `src/build.cpp` (`katka-build -A`), `src/tagbuild.cpp` | the LF pass that gives the SA samples at run boundaries; φ streaming the SA bottom-up into tag runs, without a suffix array |
| §2 "Matches and credit", §4, Table 1 | `scripts/credit_schemes.py` (Zymo), `scripts/score_silva.py` (SILVA), `scripts/norm_credit.py` | the credit schemes, applied to the listings |
| §5 "Longest MEMs", Zymo prefix test | `scripts/maxmem_roc.py` | sensitivity, specificity and accuracy per threshold on the longest MEMs |
| App. A "Parallelism" | `src/classify.cpp` (`classify_all`) | threads take chunks of 256 reads and write them back in input order |

The order in which to read the code: `src/rlbwt.hpp` → `src/search.hpp` → `src/bml.hpp` → `src/tag_array.hpp` → `src/classify.cpp`. The builders can wait.

## Building

The code needs sdsl-lite (with libdivsufsort) and zlib.

```
make SDSL=/path/to/sdsl        # the prefix given to sdsl-lite's install.sh; default ../deps/sdsl
```

This puts `katka-prep`, `katka-build`, `katka-tagbuild`, `katka-csabuild`, `katka-kbbuild` and `katka-classify` in `bin/`.

## Running

The default SILVA configuration builds one BWT with Big-BWT. `list.tsv` has one line per genus, `<tag>\t<fasta>`, in tree order.

```
katka-prep list.tsv bac                           # bac.S, bac.tbl
bigbwt bac.S && mv bac.S.bwt bac.bwt              # BWT(S$), one terminator
katka-build -G -A bac.S bac.bwt bac.tbl bac       # bac.rlbwt (RLBWT), bac.rix (SA samples, for katka-tagbuild)
katka-tagbuild bac 1 bac.s1.tag                   # tag = genome (here: genus) index; bac.s1.tag + .rmq
katka-csabuild bac.bwt bac.csa                    # explicit RLCSA
katka-classify -j 8 -L 30 -c -C bac.csa -T bac.s1.tag -F 10 bac.rlbwt reads.fq out.listings
```

For Zymo, the tags are species and the BWT is built per dataset (`katka-prep -d` writes one file per species; Big-BWT each, then pfp-merge). Homopolymer compression is on.

```
katka-prep -d -s -U list.tsv zy
# Big-BWT on each zy.ds<j>.S, pfp-merge into zy.bwt (one terminator per dataset)
katka-build -G -A zy.S zy.bwt zy.tbl zy
katka-tagbuild zy 1 zy.tag -S -n                  # tag = species (column 2 of zy.tbl); -n: counting only, no RMQ
katka-csabuild zy.bwt zy.csa
katka-kbbuild -k 20 zy.S zy.k20.kbf
katka-classify -U -L 30 -D 14 -c -C zy.csa -T zy.tag -j 12 zy.rlbwt reads.fa out.listings   # threshold L = 30
katka-classify -U -X -L 20 -K zy.k20.kbf -c -C zy.csa -T zy.tag -j 12 zy.rlbwt reads.fa out.listings   # longest MEMs
```

The output has Cliffy's listing format, with counts: one record per read, `>name`, then `[a,b] {tag:count,...}` for each MEM, where `[a,b]` is the MEM's span in the read. Then score it:

```
python3 scripts/credit_schemes.py out.listings truth.txt zy.tbl         # Zymo, every credit scheme
python3 scripts/maxmem_roc.py out.listings truth.txt zy.tbl reads.fq    # -X runs: per-threshold sens./spec./accuracy
python3 scripts/score_silva.py --mode freq ...                          # SILVA (needs cliffy-experiments); freq = proportional
```

`katka-build` also accepts the original's six-argument form `katka-build -G S bwt - - S.tbl prefix`; the two `-` are ignored.

## What is left out

The original repository keeps everything we tried, mostly from Appendices B and C. None of it is used by the paper's main results, and it is left out here:

- minimizer digests;
- the rz-index's LZ77 grids (LCAs without tags);
- the r-index and sr-index as classifiers;
- forward–backward search;
- hybrid LCA/listing;
- verification of digested matches;
- grammar-compressed tags;
- the phrase-trimming experiments.

The options for these give an error saying so. Counting with rank queries (Appendix C) is not implemented in either version.

## Differences from the original

- **Program names:** `katka-prep`, `katka-build`, `katka-tagbuild`, `katka-csabuild`, `katka-kbbuild` and `katka-classify` replace `katka-prep`, `katka-build`, `katka-tagbuild`, `katka-csabuild`, `katka-kbbuild` and `katka-classify`.
- **The RLBWT file:** `katka-build -G` writes `prefix.rlbwt` instead of `prefix.rz`, and `katka-tagbuild` reads `prefix.rlbwt` (for s > 1). The format is the same, so an original `.rz` file can be renamed or linked. `katka-classify` takes the file by name and reads either, including an original `.rz` with the rz-index's grids, which it ignores.
- **Counting:** `katka-classify -c` replaces `rz-classify -l` with `RZ_COUNTS=1` in the environment. `-l` alone still lists the distinct tags, and neither option gives the range answers.
- **Kept:** the format tags inside the files (`RZTAG002`, `RZCSA001` and so on) are unchanged, so both versions read each other's `.tag`, `.rmq`, `.csa`, `.kbf` and `.rix` files.

## Checking against the original

On the test collection (`test-small`, 40 genomes of 4 species), the original and this version give byte-identical results.

- **Build outputs:** `katka-prep` (plain, `-d`, `-s`, `-U`, with empty and multi-record FASTA files), `katka-build -G -A` (one terminator and one per dataset), `katka-tagbuild` (s = 1, 2, 3, 4, 16; with and without `-S` and `-n`), `katka-csabuild` (explicit and Elias–Fano) and `katka-kbbuild`.
- **Classification outputs:**
  - reads in FASTQ, multi-line FASTA, lower case and with N;
  - answers by count, list and range, from sampled and unsampled tags;
  - search on the RLBWT and on both RLCSAs, with and without lookup tables;
  - with `-X`, `-K`, `-D`, `-U` and their combinations;
  - with 1, 3 and 4 threads.

On the full homopolymer-compressed Zymo index (581,802 reads; the original's `.rz`, `.tag`, `.csa` and `.kbf` files, read directly), `katka-classify` gives output byte-identical to the original's for:

- L = 30 with `-D 14`, on all reads and on a 1% subsample;
- the longest MEMs (`-X -L 1 -D 14`);
- the longest MEMs with the KeBaB filter (`-X -L 20 -K ... -D 14`).

With 12 threads on the laptop these took 675 s, 199 s and 110 s, against about 780 s, 211 s and 118 s for the original.
