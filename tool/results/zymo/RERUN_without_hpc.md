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
K=~/rz/MEM-finding/tool/bin/katka; X="/run/media/travis/Travis backup/spu2_plain"; W="/run/media/travis/Travis backup/zymo_plain_new"
cap() { systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "%e s wall, peak %M KB" "$@"; }
mkdir -p "$W" && cd "$W"
# 1. The text without -U, checked against KATKA's zp.S without writing it twice. The per-dataset files go to
#    /dev/null, and the text goes through a pipe. (~/zymo/list.tsv lists the same 7 species files as gzipped FASTA;
#    if the table differs from $X/zp.tbl in more than the paths, use $X/list.tsv instead.)
for j in 0 1 2 3 4 5 6; do ln -sf /dev/null zp.ds$j.S; done
mkfifo zp.S.fifo; ln -sf zp.S.fifo zp.S
cap $K prep -d -s ~/zymo/list.tsv zp & cmp zp.S.fifo "$X/zp.S" && echo "text identical to KATKA's zp.S"; wait
rm zp.S zp.S.fifo; ln -sf "$X/zp.S" zp.S; ln -sf "$X/zp.bwt" zp.bwt
diff <(cut -f1,2,4,5 zp.tbl) <(cut -f1,2,4,5 "$X/zp.tbl") && echo "table identical (paths aside)"
# 2. The index, in separate steps (as for zy): RLBWT and RLCSA, tags (with KATKA's converted tag runs if they exist
#    for zp, else computed: an LF pass and a phi pass over 53.7 G rows), KeBaB filter
cap $K index --csa zp
cap $K index --tags --by-tag --no-rmq zp
cap $K index --kebab 20 zp
# 3. The runs of results/zymo/ without -U, theta = 16
R=~/zymo/all.fq
cap $K find -L 30 -D 16 -c --csa -j 12 zp $R L30.plain                    # threshold L = 30
cap $K find -X -L 20 -K 20 -D 16 -c --csa -j 12 zp $R X.plain             # longest MEMs, KeBaB
cap $K find -X -L 1 -D 16 -c --csa -j 12 zp $R X1.plain                   # longest MEMs, no filter (exact)
cap $K find -L 30 -D 16 -c -j 12 zp ~/zymo/pfx/p720.fq p720.L30.plain     # prefix test (paper: 97.24% / 97.98%)
# the KeBaB rules (results/zymo_kebab/): the same commands without -U and with -D 16
```

Scoring is as in `results/zymo/`, with lengths in bases. The scripts that take the reads score each MEM by its
homopolymer-compressed length unless given `--plain`, or `--info zp.info` (plain when the index's `hpc` line is 0).
With either option, a MEM's length is end - start + 1:
- `credit_schemes.py` and `score_spu2.py` for L = 30 (they do not use lengths);
- `scripts/maxmem_roc_stream.py --info zp.info` for the longest MEMs. Its output is then that of
  `katka-readable/scripts/maxmem_roc.py` without the reads (thresholds in steps of 5, not 2), with the reads
  streamed;
- `scripts/same_longest_mems.py --info zp.info X.plain X1.plain $R` for the check that KeBaB keeps the longest MEMs;
- for the KeBaB rules (`results/zymo_kebab/run.sh`, `check_run.sh`), the same scripts with `--info zp.info`
  (`maxmem_roc_stream.py`, `compare_longest_t.py`, and `results/zymo_kebab/topt_lengths.py`), and the tbl `zp.tbl`.

```bash
S=~/rz/MEM-finding/tool/scripts; T=~/zymo/all.truth; I="$W/zp.info"; R=~/zymo/all.fq
python3 -I $S/maxmem_roc_stream.py --info $I X.plain $T zp.tbl $R
python3 -I $S/same_longest_mems.py --info $I X.plain X1.plain $R
python3 -I $S/compare_longest_t.py --info $I c.3.plain d.3.plain $R 3     # e.g., for the KeBaB rules
```

Without either option these scripts would compress homopolymers that the index did not compress, and report the
wrong lengths. The plain mode was checked on `tests/gen_data.py` data, with indexes built with and without `-U`.
Without `-U`, `maxmem_roc_stream.py --plain` gives output byte-identical to `maxmem_roc.py` without the reads. With
`--info` on a `-U` index, all four scripts give their default output, which is unchanged.

Disk: about 60 GB on the external drive for the new index files and outputs. The text and BWT are reused there.

Not yet checked:
- whether `zp.S` equals a fresh `prep -d -s` (step 1);
- how long the tag runs take without KATKA's cache: about 1.4 times Zymo's compressed size, in rows.
