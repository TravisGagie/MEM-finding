#!/bin/bash
# The desktop's SILVA queue, one step at a time, writing its outputs where scripts/silva_tables.py reads them
# (tool/NOTES.md, "SILVA queue: runner and tables").
#
#   silva_queue.sh STEP [STEP ...]        STEP: timing credit blam grammar digests uniform (or all), collect
#
# Environment (defaults for the desktop):
#   Q       the output directory (default ~/tool/silva/queue).  The finds' outputs (*.out) are large; everything
#           else (logs, scores, summaries, sizes.txt) is small, and is all scripts/silva_tables.py reads:
#           the step `collect` copies them to OUT (default tool/results/silva_queue), to commit
#   S       katka/silva (the reference in $S/work/ref, the taxonomy in $S/exp1_data)
#   IDX     the undigested index (default ~/tool/silva/bac: katka index --preset silva, see NOTES "SILVA")
#   DG      the digested index prefix (default ~/tool/silva/dg; built by the digests step if missing)
#   READS   the reads: READS/<region>_{1,2}.fq (default $S/work/compare/reads_every50, 200,000 pairs per region)
#   TRUTH   the truth files: TRUTH/<region>/<region>_seqtax.txt (default $S/work/reads/aquatic)
#   EXP     cliffy-experiments/src (default katka/cliffy-experiments/src), for the scoring
#   REGIONS (default "V1_V2 V3_V4 V4_V4 V4_V5"); PIN (a CPU for single-threaded runs, default 9); THREADS (12)
#   BIGBWT (default katka/deps/Big-BWT/bigbwt; without it, katka bwt), BIGREPAIR (bigrepair, for the grammar step);
#   MM, FASTA (the uniform step: MicrobeMixer, default katka/MicrobeMixer, and the SILVA reference it simulates
#   from, default $S/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta, as in katka/silva/run_silva.sh); MEMMAX (11G)
#
# The timing protocol (every timed run): `katka find` with one thread pinned to PIN, then with THREADS threads,
# on each region and mate; the tool's `timing:` line gives the mean per-read time (one thread) and the total time
# over the number of reads; reads per second with THREADS threads is the number of reads over the total time.
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); TOP=$(cd "$HERE/.." && pwd); K=$HERE/bin/katka
Q=${Q:-$HOME/tool/silva/queue}; S=${S:-$TOP/katka/silva}; IDX=${IDX:-$HOME/tool/silva/bac}; DG=${DG:-$HOME/tool/silva/dg}
READS=${READS:-$S/work/compare/reads_every50}; TRUTH=${TRUTH:-$S/work/reads/aquatic}
EXP=${EXP:-$TOP/katka/cliffy-experiments/src}; TAX=${TAX:-$S/exp1_data/tax_slv_ssu_138.1.txt}
D2T=${D2T:-$S/work/ref/rna/doc_to_traversal.txt}; T2L=${T2L:-$S/work/ref/trav_to_length.txt}
MM=${MM:-$TOP/katka/MicrobeMixer}; FASTA=${FASTA:-$S/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta}
BIGBWT=${BIGBWT-$TOP/katka/deps/Big-BWT/bigbwt}; [ -x "$BIGBWT" ] || BIGBWT=""
REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; PIN=${PIN:-9}; THREADS=${THREADS:-12}
CAP=""; systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0 true 2> /dev/null && CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
TS=""; command -v taskset > /dev/null && TS="taskset -c $PIN"
mkdir -p "$Q"
stamp() { echo "[$(date '+%F %T')] $*"; }
need() { for f in "$@"; do [ -e "$f" ] || { echo "missing $f"; exit 1; }; done; }
peak() {   # log, command...: wall time and peak memory appended to the log
  local log=$1; shift
  python3 - "$log" "$@" <<'PY'
import resource, subprocess, sys, time
log, cmd = sys.argv[1], sys.argv[2:]
t = time.time()
with open(log, 'a') as f:
    r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    f.write('wall %.1f s, peak %.2f GB, exit %d\n' % (time.time() - t, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024 / 1e9, r.returncode))
sys.exit(r.returncode)
PY
}
sizes() {   # name, files...: their total size, appended to Q/sizes.txt (replacing an earlier line of that name)
  local name=$1; shift
  local b; b=$(du -cbL "$@" | tail -1 | cut -f1)
  touch "$Q/sizes.txt"; grep -v "^$name	" "$Q/sizes.txt" > "$Q/sizes.tmp" || true; mv "$Q/sizes.tmp" "$Q/sizes.txt"
  printf '%s\t%s\t%s\n' "$name" "$b" "$*" >> "$Q/sizes.txt"
}
score() {   # mode, mate-1 listings, mate-2 listings, region, output csv
  local extra=""; [ $1 = norm ] && extra="--nseq $Q/nseq.txt"
  python3 "$HERE/scripts/score_silva.py" --mode $1 --mate1-listings "$2" --mate2-listings "$3" --doc-id-to-traversal "$D2T" \
      --silva-tax-ranks "$TAX" --readset-truthset "$TRUTH/$4/$4_seqtax.txt" --trav-to-length "$T2L" --exp-src "$EXP" \
      --output "$5" $extra > "$5.log" 2>&1 || { echo "scoring $1 $4 failed (see $5.log)"; exit 1; }
}
nseq() {   # sequences per document (for the normalized scheme), from the reference's per-genus files
  [ -s "$Q/nseq.txt" ] && return 0
  local n; n=$(cut -d' ' -f1 "$S/work/ref/done")
  for i in $(seq 1 $n); do grep -c '^>' "$S/work/ref/dna/doc_${i}_seq.fa"; done > "$Q/nseq.txt"
}
timed() {   # dir, name, index prefix, find options...: the timing protocol on every region and mate
  local d=$1 name=$2 P=$3; shift 3
  mkdir -p "$d"
  for r in $REGIONS; do for m in 1 2; do
    local in=$READS/${r}_$m.fq
    [ -s "$d/$name.${r}_$m.t1.log" ] || { stamp "$name $r mate $m, 1 thread"; peak "$d/$name.${r}_$m.t1.log" $CAP $TS $K find "$@" -j 1 "$P" "$in" "$d/$name.${r}_$m.out"; }
    [ -s "$d/$name.${r}_$m.t$THREADS.log" ] || { stamp "$name $r mate $m, $THREADS threads"; peak "$d/$name.${r}_$m.t$THREADS.log" $CAP $K find "$@" -j $THREADS "$P" "$in" "$d/$name.${r}_$m.out$THREADS"
      cmp -s "$d/$name.${r}_$m.out" "$d/$name.${r}_$m.out$THREADS" || { echo "$name $r $m: outputs differ with $THREADS threads"; exit 1; }; rm -f "$d/$name.${r}_$m.out$THREADS"; }
  done; done
}

step_timing() {   # the default configuration (index --preset silva; find --preset silva)
  need "$IDX.rlbwt" "$IDX.csa" "$IDX.tag"
  timed "$Q/timing" default "$IDX" --preset silva
  sizes default "$IDX.rlbwt" "$IDX.csa" "$IDX.tag" $( [ -e "$IDX.tag.rmq" ] && echo "$IDX.tag.rmq")
}
step_credit() {   # the credit schemes on the same MEMs (the default run's counts), every region
  [ -s "$Q/timing/default.V1_V2_1.out" ] || step_timing
  nseq
  mkdir -p "$Q/credit"
  for r in $REGIONS; do for mode in lca one cliffs list norm freq; do
    [ -s "$Q/credit/$mode.$r.csv" ] || { stamp "credit $mode $r"; score $mode "$Q/timing/default.${r}_1.out" "$Q/timing/default.${r}_2.out" $r "$Q/credit/$mode.$r.csv"; }
  done; done
}
step_blam() {     # B_lambda and the lemma's bound (the per-genus tags of the default index)
  mkdir -p "$Q/blam"
  [ -s "$Q/blam/blam.log" ] || { stamp "B_lambda"; peak "$Q/blam/blam.log" $CAP $K index --tags --no-rmq --tag-out "$Q/blam/b.tag" --tag-lcp 0,10,30,100,1000 "$IDX"; rm -f "$Q/blam/b.tag"; }
}
step_grammar() {  # grammar tags: built-in RePair against bigrepair, then the timing protocol with each
  need "$IDX.tagruns"
  [ -s "$Q/grammar/summary.txt" ] || { stamp "grammar tags"; THREADS=1 bash "$HERE/experiments/grammar_tags.sh" "$IDX" "$READS/V1_V2_1.fq" "$Q/grammar" -L 30 --csa -F 10; }
  timed "$Q/timing" grammar "$IDX" --preset silva -T "$Q/grammar/g.tag"
  local G="grammar"; [ -s "$Q/grammar/bg.tag" ] && G="grammar bigrepair" || echo "no bigrepair tag array (see $Q/grammar/summary.txt): only the built-in RePair is timed"
  [ -s "$Q/grammar/bg.tag" ] && timed "$Q/timing" bigrepair "$IDX" --preset silva -T "$Q/grammar/bg.tag"
  sizes grammar "$IDX.rlbwt" "$IDX.csa" "$Q/grammar/g.tag"
  [ -s "$Q/grammar/bg.tag" ] && sizes bigrepair "$IDX.rlbwt" "$IDX.csa" "$Q/grammar/bg.tag"
  for t in $G; do for r in $REGIONS; do for m in 1 2; do
    cmp -s "$Q/timing/default.${r}_$m.out" "$Q/timing/$t.${r}_$m.out" || { echo "$t ${r}_$m: answers differ from the default's"; exit 1; }
  done; done; done
  for t in $G; do for r in $REGIONS; do
    [ -s "$Q/credit/freq.$t.$r.csv" ] || score freq "$Q/timing/$t.${r}_1.out" "$Q/timing/$t.${r}_2.out" $r "$Q/credit/freq.$t.$r.csv"
  done; done
}
step_digests() {  # the digested index (App. B): explicit and Elias-Fano RLCSA with tags, RLBWT with sampled tags
  if [ ! -s "$DG.tag" ]; then
    stamp "the digested index"
    local L=$(dirname "$DG")/list.tsv; need "$L"
    mkdir -p "$Q/digests"
    peak "$Q/digests/build.log" $CAP $K prep --digest 4,11 "$L" "$DG"
    if [ -n "${BIGBWT:-}" ]; then peak "$Q/digests/build.log" $CAP "$BIGBWT" "$DG.S"; mv "$DG.S.bwt" "$DG.bwt"; rm -f "$DG.S.log"
    else peak "$Q/digests/build.log" $CAP $K bwt "$DG"; fi
    peak "$Q/digests/build.log" $CAP $K index --csa --efcsa --tags "$DG"
    peak "$Q/digests/build.log" $CAP $K index --tags --sample 4 --tag-out "$DG.s4.tag" "$DG"
  fi
  sizes digests.csa "$DG.rlbwt" "$DG.csa" "$DG.tag" "$DG.map"
  sizes digests.efcsa "$DG.rlbwt" "$DG.efcsa" "$DG.tag" "$DG.map"
  sizes digests.rlbwt.s4 "$DG.rlbwt" "$DG.s4.tag" "$DG.map"
  for L in 15 20 30 40 50; do
    timed "$Q/digests" csa.L$L "$DG" -L $L -c --csa -F 2
    timed "$Q/digests" efcsa.L$L "$DG" -L $L -c --efcsa -F 2
    [ $L -ge 30 ] && timed "$Q/digests" rlbwt.s4.L$L "$DG" -L $L -c -F 2 -T "$DG.s4.tag"
    for r in $REGIONS; do for m in 1 2; do   # the same MEMs and answers whichever the locate structure
      for v in efcsa rlbwt.s4; do f=$Q/digests/$v.L$L.${r}_$m.out; [ -e "$f" ] || continue
        if cmp -s "$Q/digests/csa.L$L.${r}_$m.out" "$f"; then echo "$v.L$L.${r}_$m identical" >> "$Q/digests/identical.log"; rm -f "$f"
        else echo "$v.L$L.${r}_$m DIFFERS from csa.L$L.${r}_$m"; exit 1; fi
      done
    done; done
    for r in $REGIONS; do
      [ -s "$Q/digests/freq.L$L.$r.csv" ] || score freq "$Q/digests/csa.L$L.${r}_1.out" "$Q/digests/csa.L$L.${r}_2.out" $r "$Q/digests/freq.L$L.$r.csv"
    done
  done
  echo "yes, the outputs of all $(sort -u "$Q/digests/identical.log" | wc -l) of their runs were byte-identical to the explicit RLCSA's" > "$Q/digests/identical.txt"
  for L in 15 20 30 40; do   # the undigested default configuration at each L
    timed "$Q/digests" undigested.L$L "$IDX" --preset silva -L $L
    for r in $REGIONS; do
      [ -s "$Q/digests/freq.undigested.L$L.$r.csv" ] || score freq "$Q/digests/undigested.L$L.${r}_1.out" "$Q/digests/undigested.L$L.${r}_2.out" $r "$Q/digests/freq.undigested.L$L.$r.csv"
    done
  done
}
step_uniform() {  # the uniform-abundance simulation and its accuracy (proportional credit)
  need "$MM/src/sim_16s_reads.py" "$FASTA"
  KATKA=$IDX EXP=$EXP TRAV_TO_LENGTH=$T2L bash "$HERE/experiments/uniform_microbemixer.sh" "$MM" "$FASTA" "$TAX" "$D2T" "$Q/uniform/work" ${NREADS:-200000} "$REGIONS"
  cp "$Q/uniform/work/summary.txt" "$Q/uniform/summary.txt"
  for r in $REGIONS; do cp "$Q/uniform/work/$r/score.csv" "$Q/uniform/$r.score.csv"; cp "$Q/uniform/work/$r/counts.txt" "$Q/uniform/$r.counts.txt"; cp "$Q/uniform/work/$r/dropped.txt" "$Q/uniform/$r.dropped.txt"; done
}

step_collect() {  # the small outputs (all that scripts/silva_tables.py reads), to OUT
  local out=${OUT:-$HERE/results/silva_queue}; mkdir -p "$out"
  (cd "$Q" && find . -type f ! -name '*.out' ! -name '*.out[0-9]*' ! -path './uniform/work/*' -print0 | xargs -0 -r cp --parents -t "$out")
  stamp "copied $(cd "$out" && find . -type f | wc -l) files ($(du -sh "$out" | cut -f1)) to $out"
}

[ $# -gt 0 ] || { echo "usage: silva_queue.sh timing|credit|blam|grammar|digests|uniform|all|collect ..."; exit 1; }
for s in "$@"; do
  case $s in
    all) for t in timing credit blam grammar digests uniform; do step_$t; done ;;
    timing|credit|blam|grammar|digests|uniform|collect) step_$s ;;
    *) echo "unknown step $s"; exit 1 ;;
  esac
done
stamp "done: $*; now: python3 $HERE/scripts/silva_tables.py $Q"
