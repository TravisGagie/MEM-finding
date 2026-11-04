#!/bin/bash
# A uniform-abundance MicrobeMixer simulation of 16S reads over the genera of a reference, and its scoring with
# katka (tool/NOTES.md, "Uniform-abundance simulation").
#
#   uniform_microbemixer.sh MICROBEMIXER SILVA_FASTA SILVA_TAX DOC_TO_TRAVERSAL OUT [NREADS] [REGIONS]
#
# MICROBEMIXER: a clone of github.com/oma219/MicrobeMixer (its sim_16s_reads.py needs art_illumina, from ART,
# and the Python module regex); SILVA_FASTA, SILVA_TAX: the SILVA reference and taxonomy the reads are
# simulated from (as in katka/silva/run_silva.sh); DOC_TO_TRAVERSAL: the reference's document -> genus file
# (katka/silva/work/ref/rna/doc_to_traversal.txt), whose genera are the ones simulated.  NREADS (default
# 200000) read pairs per region; REGIONS (default "V1_V2 V3_V4 V4_V4 V4_V5").
#
# Uniform abundance: every genus of the reference is listed with a share of 0.  MicrobeMixer gives each genus
# its share of the reads still left, then spreads what is left evenly over the genera from which it could
# extract the region, so equal nonzero shares would give the first genus about e times the reads of the
# last; with shares of 0, all reads are spread evenly (to within one read) over those genera.  Genera whose
# sequences have no region between the primers get no reads, as in any MicrobeMixer simulation.
#
# Genera whose extracted sequences give ART too few reads (MicrobeMixer's assertion "not enough reads in file") are
# dropped from that region's table, round by round, and listed in OUT/<region>/dropped.txt.
#
# Output: OUT/uniform.tsv (the abundance table; OUT/<region>/table.tsv after drops), OUT/<region>/<region>_mate_{1,2}.fq and _seqtax.txt (the
# truth, MicrobeMixer's), OUT/<region>/counts.txt (reads per genus, checked to be uniform).
# With KATKA=PREFIX (an index built as for SILVA, katka index --preset silva) and EXP (cliffy-experiments/src)
# and TRAV_TO_LENGTH (katka/silva/work/ref/trav_to_length.txt), it also runs `katka find --preset silva` on
# each mate and scores the listings with scripts/score_silva.py under proportional credit (MODE, default
# freq), writing OUT/<region>/score.csv and OUT/summary.txt.
set -euo pipefail
MM=$(cd "$1" && pwd); FASTA=$2; TAX=$3; D2T=$4; mkdir -p "$5"; OUT=$(cd "$5" && pwd)
NREADS=${6:-200000}; REGIONS=${7:-"V1_V2 V3_V4 V4_V4 V4_V5"}; MODE=${MODE:-freq}
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
command -v art_illumina > /dev/null || { echo "MicrobeMixer needs art_illumina (ART; e.g. apt install art-nextgen-simulation-tools)"; exit 1; }
python3 -c "import regex" 2> /dev/null || { echo "MicrobeMixer needs the Python module regex (pip install regex)"; exit 1; }
# the abundance table: every genus of the reference, share 0 (MicrobeMixer's format: full_path, read_count, percent_total)
python3 - "$D2T" "$OUT/uniform.tsv" <<'PY'
import sys
genera = []
for line in open(sys.argv[1]):
    t = line.split(None, 1)
    if len(t) == 2:
        g = t[1].strip()
        if g not in genera: genera.append(g)
with open(sys.argv[2], 'w') as o:
    o.write('full_path\tread_count\tpercent_total\n')
    for g in genera: o.write('%s\t1\t0.0\n' % g)
print('%d genera in the reference' % len(genera), file=sys.stderr)
PY
SUM=$OUT/summary.txt; : > "$SUM"; : > "$OUT/dropped.txt"
for r in $REGIONS; do
  d=$OUT/$r; p=${r/V4_V4/V4}; P=$MM/data/${p}_primers.txt
  [ -s "$P" ] || { echo "no primers $P"; exit 1; }
  if [ ! -s "$d/done" ]; then
    # MicrobeMixer asks ART for ceil(needed / sequences) * 2 reads per extracted sequence and asserts that a genus
    # then has more reads than it needs; ART makes no reads from a sequence shorter than the reads, so a genus with
    # few long enough sequences can fail the assertion.  Such genera are dropped from this region's table and
    # MicrobeMixer is rerun, until it passes (at most 10 rounds); OUT/<region>/dropped.txt lists them.
    rm -rf "$d"; mkdir -p "$d"; cp "$OUT/uniform.tsv" "$d/table.tsv"; : > "$d/dropped.txt"
    for round in 1 2 3 4 5 6 7 8 9 10; do
      if python3 "$MM/src/sim_16s_reads.py" simulate --biome-abundance "$d/table.tsv" --primers "$P" --silva-taxonomy "$TAX" \
          --silva-ref "$FASTA" --num-reads "$NREADS" --output-name "$r" --temp-dir "$d/" > "$d/stdout.log" 2>&1; then break; fi
      grep -q "not enough reads in file" "$d/stdout.log" && [ $round -lt 10 ] \
        || { echo "MicrobeMixer failed for $r in round $round (see $d/stdout.log)"; exit 1; }
      python3 - "$d" "$NREADS" <<'PY'
import os, sys
d, N = sys.argv[1], int(sys.argv[2])
rows = open(d + '/table.tsv').read().splitlines()
head, rows = rows[0], rows[1:]
def nreads(i):
    f = '%s/genus_%d_reads_1.fq' % (d, i)
    return sum(1 for _ in open(f)) // 4 if os.path.exists(f) else 0
# the genera MicrobeMixer gives reads (those with extracted sequences), and what it asks of each: with every share
# 0, N // live each, plus one for the first N % live in table order (its calculate_read_count_per_genus)
seqs = lambda i: '%s/genus_%d_seqs.fna' % (d, i)
live = [i for i in range(1, len(rows) + 1) if os.path.exists(seqs(i)) and os.path.getsize(seqs(i)) > 0]
need = {i: N // len(live) + (k < N % len(live)) for k, i in enumerate(live)}
bad = [(i, nreads(i), need[i]) for i in live if nreads(i) <= need[i]]
if not bad:
    sys.exit('no genus is short of reads, but MicrobeMixer said "not enough reads in file"')
with open(d + '/dropped.txt', 'a') as o:
    for i, k, n in bad:
        o.write('%s\t%d reads, %d needed\n' % (rows[i - 1].split('\t')[0], k, n))
drop = {i for i, _, _ in bad}
with open(d + '/table.tsv', 'w') as o:
    o.write(head + '\n' + ''.join(row + '\n' for i, row in enumerate(rows, 1) if i not in drop))
print('%s: dropped %d genera short of reads, %d left' % (os.path.basename(d), len(bad), len(rows) - len(drop)), file=sys.stderr)
PY
      rm -f "$d"/genus_* "$d/${r}_mate_1.fq" "$d/${r}_mate_2.fq"
    done
    rm -f "$d"/genus_*_reads_1.fq "$d"/genus_*_reads_2.fq "$d"/genus_*_seqs.fna "$d"/genus_*.aln
    echo "$r: $(wc -l < "$d/dropped.txt") genera dropped for too few reads (listed in $d/dropped.txt)" | tee -a "$OUT/dropped.txt"
    date > "$d/done"
  fi
  # reads per genus: uniform to within one read over the genera that got any
  python3 - "$d/${r}_seqtax.txt" "$d/counts.txt" <<'PY' | tee -a "$SUM"
import sys
c = [int(l.rstrip('\n').split('\t')[2]) for l in open(sys.argv[1]) if l.strip()]
nz = [x for x in c if x > 0]
open(sys.argv[2], 'w').write('\n'.join(map(str, c)) + '\n')
ok = max(nz) - min(nz) <= 1 if nz else False
print('%s: %d read pairs over %d of %d genera (%d to %d each): %s' % (sys.argv[1].split('/')[-2], sum(c), len(nz), len(c),
      min(nz) if nz else 0, max(nz) if nz else 0, 'uniform' if ok else 'NOT UNIFORM'))
if not ok: sys.exit(1)
PY
  if [ -n "${KATKA:-}" ]; then
    for m in 1 2; do [ -s "$d/out_$m" ] || $K find --preset silva -j "${THREADS:-$(nproc)}" "$KATKA" "$d/${r}_mate_$m.fq" "$d/out_$m" 2> "$d/find_$m.log"; done
    python3 "$HERE/scripts/score_silva.py" --mode "$MODE" --mate1-listings "$d/out_1" --mate2-listings "$d/out_2" \
        --doc-id-to-traversal "$D2T" --silva-tax-ranks "$TAX" --readset-truthset "$d/${r}_seqtax.txt" \
        --trav-to-length "$TRAV_TO_LENGTH" --exp-src "$EXP" --output "$d/score.csv" > "$d/score.log" 2>&1 \
      || { echo "scoring failed for $r (see $d/score.log)"; exit 1; }
    echo "  $(grep ",genus," "$d/score.csv")" | tee -a "$SUM"
  fi
done
echo "done: $SUM"
