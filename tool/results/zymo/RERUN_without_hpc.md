# Rerunning the Zymo results without homopolymer compression (not run; commands kept ready, 2026-11-04)

The Zymo results above use homopolymer compression (`-U`) and θ = 14. Without compression, the paper's
low-complexity filter is θ = 16 of the 64 3-mers. KATKA's laptop build of the uncompressed index is on the
external drive, under `/run/media/travis/Travis backup/spu2_plain/` (`$X` below):
- `zp.S` (53,741,194,206 characters) and `zp.bwt` (Big-BWT per dataset, then pfp-merge, 7 datasets);
- `zp.csa`, in KATKA's format;
- the per-dataset texts and BWTs;
- its reference split into `ref/*.fa` (`list.tsv`).

The new tool can reuse `zp.bwt` once its text is shown to be identical to a fresh `prep` without `-U`, as was done
for `zy.bwt`.

```bash
K=~/research/mem-finding/tool/bin/katka; X="/run/media/travis/Travis backup/spu2_plain"; W="/run/media/travis/Travis backup/zymo_plain_new"
cap() { systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "%e s wall, peak %M KB" "$@"; }
mkdir -p "$W" && cd "$W"
# 1. The text without -U, checked against KATKA's zp.S without writing it twice. The per-dataset files go to
#    /dev/null, and the text goes to standard output, through a pipe into cmp. (~/zymo/list.tsv lists the same 7
#    species files as gzipped FASTA; if the table differs from $X/zp.tbl in more than the paths, use $X/list.tsv.)
#    No -U here: prep writes zp.info with hpc 0, and find then reads uncompressed (its default follows zp.info).
#    Only if the texts are identical are KATKA's zp.S and zp.bwt linked in; otherwise stop here.
for j in 0 1 2 3 4 5 6; do ln -sf /dev/null zp.ds$j.S; done
ln -sf /dev/stdout zp.S
if cap $K prep -d -s ~/zymo/list.tsv zp | cmp - "$X/zp.S"; then echo "text identical to KATKA's zp.S"
  rm zp.S; ln -sf "$X/zp.S" zp.S; ln -sf "$X/zp.bwt" zp.bwt
else echo "TEXT DIFFERS (or prep failed): do not reuse $X/zp.bwt"; fi
diff <(cut -f1,2,4,5 zp.tbl) <(cut -f1,2,4,5 "$X/zp.tbl") && echo "table identical (paths aside)"
cat zp.info                                   # must say hpc 0
# 2. The index, in separate steps (as for zy): RLBWT and RLCSA, tags (with KATKA's converted tag runs if they exist
#    for zp, else computed: an LF pass and a phi pass over 53.7 G rows), KeBaB filter
cap $K index --csa zp
cap $K index --tags --by-tag --no-rmq zp
cap $K index --kebab 20 zp
# 3. The runs of results/zymo/ without -U, theta = 16 (do not use --preset zymo here: it adds -U and -D 14;
#    find warns "-U, but the index was built without homopolymer compression" if -U slips in)
R=~/zymo/all.fq
cap $K find -L 30 -D 16 -c --csa -j 12 zp $R L30.plain                    # threshold L = 30
cap $K find -X -L 20 -K 20 -D 16 -c --csa -j 12 zp $R X.plain             # longest MEMs, KeBaB
cap $K find -X -L 1 -D 16 -c --csa -j 12 zp $R X1.plain                   # longest MEMs, no filter (exact)
cap $K find -L 30 -D 16 -c -j 12 zp ~/zymo/pfx/p720.fq p720.L30.plain     # prefix test (paper: 97.24% / 97.98%)
# the KeBaB rules (results/zymo_kebab/): the same commands with Z="-D 16 -c --csa -j 12 --stats" (no -U), and
#    H and B pointing at this zp and zp.tbl: "$W" has a space, so quote $H and $B in those scripts or use a
#    symlink without one (e.g. ln -s "$W" ~/zymo/plain; H=~/zymo/plain/zp); and see the scoring notes below
# 4. Scoring (as in results/zymo/, with the paths of this run)
KR=~/research/mem-finding/katka-readable/scripts; T=~/zymo/all.truth; B=zp.tbl    # in "$W"; its path has a space
python3 -I $KR/credit_schemes.py L30.plain $T $B
python3 -I ~/research/mem-finding/tool/results/zymo/run_scripts/score_spu2.py L30.plain $T $B
python3 -I $KR/maxmem_roc.py X.plain $T $B  | head -8    # no reads argument: MEM lengths in bases
python3 -I $KR/maxmem_roc.py X1.plain $T $B | head -8
python3 -I - ~/zymo/pfx/truth p720.L30.plain <<'PY'
import sys
truth = dict(l.split()[:2] for l in open(sys.argv[1]))
has = {}; name = None
for line in open(sys.argv[2]):
    if line[0] == '>': name = line[1:].strip(); has[name] = False
    elif '[' in line: has[name] = True
pos = [k for k in has if truth[k] != 'yeast']; neg = [k for k in has if truth[k] == 'yeast']
print('prefix test, 720 bases, L = 30, -D 16: %d bacterial, %d yeast prefixes' % (len(pos), len(neg)))
print('  sensitivity %.2f%%  specificity %.2f%%' % (100 * sum(has[k] for k in pos) / len(pos), 100 * sum(not has[k] for k in neg) / len(neg)))
PY
```

Scoring notes:
- `credit_schemes.py` and `score_spu2.py` take MEM lengths from the listings (`[a,b]`, in read bases), so they
  are right with or without `-U`.
- Do not use `tool/scripts/maxmem_roc_stream.py` here: it always scores each MEM by its homopolymer-compressed
  length (it exists only to stream the reads for that). KATKA's `maxmem_roc.py` without the reads argument scores
  by length in bases and keeps only one number per read in memory; it prints thresholds in steps of 5, not 2.
- The same holds for `tool/scripts/same_longest_mems.py`, `tool/scripts/compare_longest_t.py` and
  `results/zymo_kebab/topt_lengths.py`, used by the KeBaB-rules comparisons: they compare compressed lengths
  and need a plain-length mode before those comparisons are rerun on `zp`.

Disk: about 60 GB on the external drive for the new index files and outputs. The text and BWT are reused there.

Not yet checked:
- whether `zp.S` equals a fresh `prep -d -s` (step 1);
- how long the tag runs take without KATKA's cache: about 1.4 times Zymo's compressed size, in rows.
- memory: n grows by 53.74 / 39.28 = 1.37 over `zy`. If r grows about as much (it did on synthetic data, below),
  the RLCSA is about 7.6 GB instead of 5.56 GB, and `index --csa` and the `--csa` searches peak near 9-10 GB,
  inside the 11 GB cap but with little room; the size of KATKA's `$X/zp.csa` against its `zy.csa` gives the real
  growth of r before starting.

## Checked on synthetic data (2026-11-04)

The block above (steps 1-4) was run end to end in a cloud session on synthetic data from `tests/gen_data.py`
(two seeds: 8 species of 5 genomes each, 2 records per genome, as gzipped FASTA, one file per species; 7 in
`list.tsv`, tagged by species name, the 8th the negative "yeast"; 600 reads with a truth file and 311 720-base
prefixes). A separate `prep -d -s` plus `katka bwt -d` stood in for KATKA's `$X/zp.S` and `zp.bwt`. Only the
paths and `cap()` were substituted (`prlimit --as` instead of `systemd-run`, unavailable there), and `$W` had a
space in its name, as on the external drive. Every step ran, the scoring included; `find -U` on this index
printed the mismatch warning; a list file that is missing, or prep with `-U`, makes step 1 report "TEXT DIFFERS"
and link nothing.

Changed after the check:
- step 1: the fifo (`mkfifo`, `prep &`, `cmp`) hung for ever when prep failed before opening its output (cmp
  waited on the fifo); prep now writes the text to standard output (`zp.S -> /dev/stdout`), piped into cmp, and
  KATKA's `zp.S` and `zp.bwt` are linked in only if cmp succeeds;
- the scoring commands are written out (step 4), with `B=zp.tbl`: an unquoted `$W/zp.tbl` splits at the space of
  "Travis backup", and `maxmem_roc.py` then takes the second half as the reads file;
- `maxmem_roc_stream.py` replaced by KATKA's `maxmem_roc.py` without reads, since the former always scores
  compressed lengths (see the scoring notes); the prefix test's scoring added (from `run_scripts/runs2.sh`).

Not changed, but to know: the tool records compression in `PREFIX.info` at prep time; `index` has no `-U`, and
`find` defaults to what `PREFIX.info` says, so leaving `-U` out of prep is enough. `find` warns about `-U` on an
uncompressed index, but not about `--no-hpc` on a compressed one (and `find -h` lists `-U / --no-hpc` on one
line, as if they were the same option: `--no-hpc` turns compression off).

Synthetic sizes, the same 7 species with and without `-U` (random sequence, so about 3/4 of the bases remain
after compression, against 0.73 for Zymo): n 1,381,593 -> 1,838,619 (x 1.33), r 349,999 -> 466,170 (x 1.33),
tag runs 436,734 -> 532,458 (x 1.22).
