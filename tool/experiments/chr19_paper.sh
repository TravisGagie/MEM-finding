#!/bin/bash
# The chromosome 19 experiments of "More RLBWT Tricks" (arXiv:2609.40299v2, Section 5): Table 1, Table 2 and
# the locating paragraph, with this tool.  See tool/CHR19.md.
#
#   chr19_paper.sh build    DATA WORK [SIZES]   indexes for each size (default "33 65 129" sequences = 16, 32, 64
#                                               individuals): prep, BWT, index --preset chr19; at the smallest
#                                               size also the STPD, parse-oracle, RLZ and locating structures
#   chr19_paper.sh table1   DATA WORK [SIZES]   Table 1: BML (L = 1000 / 100) with plain search, two-level
#                                               without and with block steps, FB with blocks
#   chr19_paper.sh table2   DATA WORK           Table 2 (the smallest size): time, peak memory, files read
#   chr19_paper.sh locating DATA WORK           the locating paragraph: the samples' size, the times, and
#                                               positions compared with STPD over the parse, MEM by MEM
#   chr19_paper.sh all      DATA WORK [SIZES]   all four
#
# DATA holds collection.txt (one sequence per line: CHM13 first, then both haplotypes of each individual, in
# order) and queries.txt (the two query haplotypes, one per line).  WORK receives the indexes, outputs and
# logs; WORK/table1.md, WORK/table2.md and WORK/locating.md are the results.
#
# Environment:
#   BIGBWT   Big-BWT's bigbwt (for the real data: katka bwt needs about 6 bytes per character, 48 GB at 32
#            individuals); without it, katka bwt
#   PIN      a CPU to pin the searches to (taskset)
#   MEMMAX   the memory cap of each step (default 11G; needs systemd-run --user; ignored if unavailable)
#   STPD     a built clone of the published STPD-index (github.com/regindex/STPD-index, commit 8a432a4) for
#            Table 2's last two rows (experiments/stpd_index_chr19.sh); without it, those rows are skipped
set -euo pipefail
CMD=$1; DATA=$(cd "$2" && pwd); mkdir -p "$3"; W=$(cd "$3" && pwd); SIZES=${4:-"33 65 129"}
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
MEMMAX=${MEMMAX:-11G}
RUN=""; [ -n "${PIN:-}" ] && command -v taskset > /dev/null && RUN="taskset -c $PIN"
CAP=""; systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0 true 2> /dev/null && CAP="systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0"
SMALL=$(echo $SIZES | tr ' ' '\n' | sort -n | head -1)
stamp() { echo "[$(date '+%F %T')] $*"; }
peak() {   # log, command...: runs it (output to log), appends "wall W s, peak P GB, exit E"
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
mb() { python3 -c "import os,sys; print('%.0f' % (sum(os.path.getsize(f) for f in sys.argv[1:]) / 1e6))" "$@"; }
pk() { grep -o 'peak [0-9.]* GB' "$1" | tail -1 | grep -o '[0-9.]*'; }
st() { grep -o 'wall [0-9.]* s, [0-9.]* us/read' "$1" | grep -o '^wall [0-9.]*' | grep -o '[0-9.]*'; }   # search time (no loading)
[ -s "$DATA/collection.txt" ] && [ -s "$DATA/queries.txt" ] || { echo "DATA needs collection.txt and queries.txt"; exit 1; }
[ -x "$K" ] || { echo "build the tool first: make -C $HERE"; exit 1; }
Q=$W/queries.fa
[ -s "$Q" ] || tr -c 'ACGT\n' '\n' < "$DATA/queries.txt" | awk 'length($0) > 0 { print ">p" NR - 1; print }' > "$Q"

build() {   # the indexes of one size n (sequences)
  local n=$1 P=$W/c$1
  if [ ! -s "$P.lvl" ]; then
    stamp "c$n: prep, BWT, index --preset chr19"
    mkdir -p "$W/fa$n"
    head -n "$n" "$DATA/collection.txt" | awk -v d="$W/fa$n" '{ f = sprintf("%s/s%03d.fa", d, NR); print ">s" NR > f; print > f; close(f) }'
    ls "$W/fa$n"/*.fa | sed 's/^/-\t/' > "$W/list$n.tsv"
    peak "$W/build$n.log" $CAP $K prep "$W/list$n.tsv" "$P"; rm -rf "$W/fa$n"
    if [ -n "${BIGBWT:-}" ]; then peak "$W/build$n.log" $CAP "$BIGBWT" "$P.S"; mv "$P.S.bwt" "$P.bwt"; rm -f "$P.S.log"
    else peak "$W/build$n.log" $CAP $K bwt "$P"; fi
    : > "$W/index$n.log"
    peak "$W/index$n.log" $CAP $K index --preset chr19 "$P"
  fi
  if [ "$n" = "$SMALL" ] && [ ! -s "$P.tloc" ]; then
    stamp "c$n: the STPD sample (through the parse oracle), STPD over the parse, the RLZ parse oracle, two-level locating"
    peak "$W/extra$n.log" $CAP $K index --stpd --parse-oracle "$P"     # c.ptxt (2-bit dictionary) and c.stpd
    peak "$W/extra$n.log" $CAP $K index --stpd-parse "$P"              # c.pstpd
    peak "$W/extra$n.log" $CAP $K index --rlzap "$P"                   # c.rptxt
    : > "$W/tloc$n.log"
    peak "$W/tloc$n.log" $CAP $K index --twolevel-locate "$P"          # c.tloc
  fi
}

search() {   # name, prefix, options...: one single-threaded search, its MEMs compared with plain search's
  local name=$1 P=$2; shift 2
  rm -f "$W/$name.log"
  peak "$W/$name.log" $CAP $RUN $K find --format mems --stats "$@" "$P" "$Q" "$W/$name.tsv"
}
same_mems() { cmp -s <(cut -f1-3 "$W/$1.tsv") <(cut -f1-3 "$W/$2.tsv"); }

table1() {
  local out=$W/table1.md
  { echo "| | $(for n in $SIZES; do printf '%s individuals | ' $(( (n - 1) / 2 )); done)"
    echo "|---|$(for n in $SIZES; do printf -- '---|'; done)"; } > "$out.head"
  declare -A row
  for n in $SIZES; do
    local P=$W/c$n
    for L in 1000 100; do
      stamp "Table 1: c$n, L = $L"
      search t1.plain.$L.$n "$P" -L $L --csa
      search t1.nj.$L.$n "$P" -L $L --csa --twolevel --no-jumps
      search t1.blocks.$L.$n "$P" -L $L --csa --twolevel
      for c in nj blocks; do same_mems t1.$c.$L.$n t1.plain.$L.$n || { echo "c$n L = $L: $c differs from plain search"; exit 1; }; done
    done
    search t1.fb.$n "$P" -L 1000 --csa --twolevel --fb
    same_mems t1.fb.$n t1.plain.1000.$n || { echo "c$n: FB differs from plain search"; exit 1; }
    row[plain]+="$(st $W/t1.plain.1000.$n.log) / $(st $W/t1.plain.100.$n.log) | "
    row[nj]+="$(st $W/t1.nj.1000.$n.log) / $(st $W/t1.nj.100.$n.log) | "
    row[blocks]+="$(st $W/t1.blocks.1000.$n.log) / $(st $W/t1.blocks.100.$n.log) | "
    row[fb]+="$(st $W/t1.fb.$n.log) | "
    row[size]+="$(mb $P.lvl) | "
    row[peak]+="$(pk $W/index$n.log) | "
  done
  { cat "$out.head"
    echo "| plain backward search (RLCSA) | ${row[plain]}"
    echo "| two-level, no block steps | ${row[nj]}"
    echo "| two-level with blocks | ${row[blocks]}"
    echo "| forward-backward with blocks | ${row[fb]}"
    echo "| two-level index (MB) | ${row[size]}"
    echo "| peak memory of the index build (GB) | ${row[peak]}"
    echo
    echo "Search time in seconds for L = 1000 / L = 100 (single-threaded, excluding loading). MEMs at L = 1000 / 100: $(for n in $SIZES; do printf '%s / %s; ' $(wc -l < $W/t1.plain.1000.$n.tsv) $(wc -l < $W/t1.plain.100.$n.tsv); done)every configuration gave the same MEMs as plain search."; } > "$out"
  rm -f "$out.head"; cat "$out"
}

table2() {
  local n=$SMALL P=$W/c$SMALL out=$W/table2.md
  local cfgs=(
    "plain|plain backward search (RLCSA)|--csa"
    "blocks|two-level with blocks|--csa --twolevel"
    "blocks-locate|  with locating|--csa --twolevel --locate"
    "stpd-text|STPD over the text|--stpd --no-intervals"
    "stpd-text-locate|  with locating|--stpd --locate"
    "stpd-oracle|STPD over the text, read through the parse|--stpd --parse-oracle --no-intervals"
    "stpd-parse|STPD over the parse, plain parse|--csa --stpd-parse --no-intervals"
    "stpd-parse-locate|  with locating|--csa --stpd-parse --locate"
    "stpd-parse-rlz|STPD over the parse, RLZ parse|--csa --stpd-parse --rlzap-oracle --no-intervals")
  for L in 1000 100; do
    for c in "${cfgs[@]}"; do
      IFS='|' read -r key label opts <<< "$c"
      stamp "Table 2: $key, L = $L"
      search t2.$key.$L "$P" -L $L $opts
      [ $key = plain ] || same_mems t2.$key.$L t2.plain.$L || { echo "Table 2: $key differs from plain search at L = $L"; exit 1; }
    done
  done
  # files read (MB): whole files, or for STPD over the parse the components it loads (the tool lists them)
  comp() { python3 - "$1" "$2" <<'PY'
import re, sys
text = open(sys.argv[1]).read()
line = next(l for l in text.splitlines() if l.startswith('STPD over the parse, components'))
def f(pat):
    m = re.search(pat + r' ([0-9]+\.[0-9]+)', line)
    return float(m.group(1)) if m else 0.0
base = sum(f(p) for p in (r'plain parse', r'boundaries', r'dictionary strings', r'offsets', r'parse sample',
                          r'entry table', r'\+ C', r'character-level sample', r'right end: [A-Z]+', r', B', r'dictionary hash'))
loc = f(r'parse phi\^-1') + f(r'\+ phi\^-1')
rlz = re.search(r'RLZAP parse \(MB\).*: ([0-9.]+) in all', text)
mode = sys.argv[2]
r = base + (loc if mode == 'locate' else 0) + (float(rlz.group(1)) if mode == 'rlz' and rlz else 0)
print('%.0f' % r)
PY
  }
  local fp=$(mb $P.rlbwt $P.csa) fb=$(mb $P.rlbwt $P.csa $P.lvl) fl=$(mb $P.rlbwt $P.csa $P.lvl $P.tloc)
  local ft=$(mb $P.S $P.stpd $P.rlbwt) fo=$(mb $P.ptxt $P.stpd $P.rlbwt)
  local fs=$(comp $W/t2.stpd-parse.1000.log plain) fsl=$(comp $W/t2.stpd-parse-locate.1000.log locate) fsr=$(comp $W/t2.stpd-parse-rlz.1000.log rlz)
  declare -A files=([plain]=$fp [blocks]=$fb [blocks-locate]=$fl [stpd-text]=$ft [stpd-text-locate]=$ft [stpd-oracle]=$fo [stpd-parse]=$fs [stpd-parse-locate]=$fsl [stpd-parse-rlz]=$fsr)
  { echo "| configuration | time (s) | peak (GB) | files (MB) |"
    echo "|---|---|---|---|"
    for c in "${cfgs[@]}"; do
      IFS='|' read -r key label opts <<< "$c"
      echo "| $label | $(st $W/t2.$key.1000.log) / $(st $W/t2.$key.100.log) | $(pk $W/t2.$key.1000.log) | ${files[$key]} |"
    done
    if [ -n "${STPD:-}" ]; then
      stamp "Table 2: the published STPD-index (experiments/stpd_index_chr19.sh)"
      for L in 1000 100; do
        THEIRS_BML_ONLY=1 OURS=$P PIN=${PIN:-} BIGBWT=${BIGBWT:-} "$HERE/experiments/stpd_index_chr19.sh" "$STPD" "$W/list$n.tsv" "$Q" "$W/stpdindex" $L > "$W/stpdindex$L.log" 2>&1
        cp "$W/stpdindex/summary.txt" "$W/stpdindex/summary.L$L.txt"
      done
      python3 - "$W" <<'PY'
import re, sys, os
w = sys.argv[1]
def times(L):
    s = open(os.path.join(w, 'stpdindex', 'summary.L%d.txt' % L)).read()
    find = float(re.search(r'theirs bml: .*?finding them ([0-9.]+) s', s).group(1))
    m = re.search(r'theirs bml-locate: .*?finding them ([0-9.]+) s.*?locating [0-9]+ occurrences ([0-9.]+) s', s)
    return find, float(m.group(1)) + float(m.group(2))
size = os.path.getsize(os.path.join(w, 'stpdindex', 'st.stpdi')) / 1e6
a, b = times(1000), times(100)
print('| published STPD-index (RLZ text) | %.1f / %.1f | (not measured) | %.0f |' % (a[0], b[0], size))
print('|   with locating | %.1f / %.1f | (not measured) | %.0f |' % (a[1], b[1], size))
PY
    else
      echo "| published STPD-index (RLZ text) | (set STPD to run it) | | |"
    fi
    echo
    echo "$(( (n - 1) / 2 )) individuals. Search time for L = 1000 / L = 100 (excluding loading), the peak memory of the search at L = 1000, and the size of what the search reads. Every configuration gave the same MEMs as plain search."; } > "$out"
  cat "$out"
}

locating() {
  local n=$SMALL P=$W/c$SMALL out=$W/locating.md
  [ -s "$W/t2.blocks-locate.1000.log" ] || table2
  for L in 1000 100; do
    stamp "locating: positions, L = $L"
    search loc.blocks.$L "$P" -L $L --csa --twolevel --locate --locate-out "$W/loc.blocks.$L.pos"
    search loc.stpd-parse.$L "$P" -L $L --csa --stpd-parse --locate --locate-out "$W/loc.stpd-parse.$L.pos"
  done
  python3 - "$W" "$P.tloc" <<'PY' > "$out"
import os, re, sys
w, tloc = sys.argv[1], sys.argv[2]
def wall(name):
    return re.search(r'wall ([0-9.]+) s, [0-9.]+ us/read', open(os.path.join(w, name + '.log')).read()).group(1)
occ = {L: sum(int(l.split('\t')[5]) for l in open(os.path.join(w, 'loc.blocks.%d.tsv' % L))) for L in (1000, 100)}
plain_rows = {L: sum(int(l.split('\t')[4]) - int(l.split('\t')[3]) for l in open(os.path.join(w, 't2.plain.%d.tsv' % L))) for L in (1000, 100)}
same = {L: open(os.path.join(w, 'loc.blocks.%d.pos' % L), 'rb').read() == open(os.path.join(w, 'loc.stpd-parse.%d.pos' % L), 'rb').read() for L in (1000, 100)}
print('- the samples, phi^-1 and the two lists (PREFIX.tloc): %.0f MB' % (os.path.getsize(tloc) / 1e6))
print('- two-level with blocks: %s / %s s without locating, %s / %s s with (L = 1000 / 100; Table 2)' %
      (wall('t2.blocks.1000'), wall('t2.blocks.100'), wall('t2.blocks-locate.1000'), wall('t2.blocks-locate.100')))
print('- STPD over the parse with locating: %s / %s s' % (wall('t2.stpd-parse-locate.1000'), wall('t2.stpd-parse-locate.100')))
print('- occurrences of all MEMs: %d (L = 1000), %d (L = 100); equal to the sums of ep - sp of plain search: %s' %
      (occ[1000], occ[100], occ == plain_rows))
print('- positions (--locate-out) byte-identical to STPD over the parse, MEM by MEM: L = 1000 %s, L = 100 %s' % (same[1000], same[100]))
PY
  cat "$out"
}

case $CMD in
  build) for n in $SIZES; do build $n; done ;;
  table1) table1 ;;
  table2) table2 ;;
  locating) locating ;;
  all) for n in $SIZES; do build $n; done; table1; table2; locating ;;
  *) echo "usage: chr19_paper.sh build|table1|table2|locating|all DATA WORK [SIZES]"; exit 1 ;;
esac
stamp "done: $CMD"
