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

Scoring is as in `results/zymo/`:
- `credit_schemes.py` and `score_spu2.py` for L = 30;
- `scripts/maxmem_roc_stream.py` for the longest MEMs. Without `-U`, lengths are in bases.

Disk: about 60 GB on the external drive for the new index files and outputs. The text and BWT are reused there.

Not yet checked:
- whether `zp.S` equals a fresh `prep -d -s` (step 1);
- how long the tag runs take without KATKA's cache: about 1.4 times Zymo's compressed size, in rows.
