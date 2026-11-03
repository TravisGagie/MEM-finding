#!/bin/bash
# STPD over the parse against character-level STPD (parse oracle) and two-level search with blocks, on a
# katka chr19 index (NOTES.md, "STPD over the parse").
#
#   stpd_parse_chr19.sh PREFIX patterns.fa WORK
# PREFIX: a katka index built with `index --preset chr19` (PREFIX.S, .rlbwt, .csa, .lvl, .tbl).  Missing
# pieces are added next to it (nothing existing is changed):
#   PREFIX.stpd   the character-level sample: `index --stpd --parse-oracle`, which reads the text through
#                 PREFIX.ptxt instead of holding it (needs PREFIX.bwt)
#   PREFIX.ptxt, PREFIX.pstpd   `index --stpd-parse` (reads only PREFIX.S and PREFIX.lvl)
# Environment: PIN (a CPU to pin the searches to), LS (thresholds, default "1000 100"), MEMMAX (default
# 11G: each step runs in a systemd cgroup with that cap and no swap, if systemd-run works).
# Output: WORK/summary.txt; the MEMs (TSV) and logs in WORK.
set -euo pipefail
P=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); PAT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir -p "$3"; W=$(cd "$3" && pwd)
LS=${LS:-"1000 100"}; MEMMAX=${MEMMAX:-11G}
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
for f in S rlbwt csa lvl tbl; do [ -s "$P.$f" ] || { echo "missing $P.$f (katka index --preset chr19)"; exit 1; }; done
RUN=""; [ -n "${PIN:-}" ] && command -v taskset > /dev/null && RUN="taskset -c $PIN"
CAP=""; systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0 true 2> /dev/null && CAP="systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0"
SUM=$W/summary.txt
note() { echo "$*" | tee -a "$SUM"; }
stamp() { echo "[$(date '+%F %T')] $*"; }
peak() {   # runs a command (output to the log named first), appending wall time and peak memory
  local log=$1; shift
  python3 - "$log" "$@" <<'PY'
import os, resource, subprocess, sys, time
log, cmd = sys.argv[1], sys.argv[2:]
t = time.time()
with open(log, 'a') as f:
    r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    f.write('wall %.1f s, peak memory %.2f GB, exit %d\n' % (time.time() - t, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024 / 1e9, r.returncode))
sys.exit(r.returncode)
PY
}
: > "$SUM"
note "STPD over the parse on $P; katka $(git -C "$HERE" log --oneline -1 2>/dev/null | cut -c1-60); patterns $PAT; L in $LS"

# ---- the components
if [ ! -s "$P.stpd" ]; then
  [ -s "$P.bwt" ] || { echo "missing $P.stpd and $P.bwt (the character-level sample needs the BWT)"; exit 1; }
  stamp "character-level sample, through the parse oracle"
  peak "$W/build_stpd.log" $CAP $K index --stpd --parse-oracle "$P"
  note "index --stpd --parse-oracle: $(tail -1 "$W/build_stpd.log")"
fi
if [ ! -s "$P.pstpd" ] || [ "$(head -c 8 "$P.ptxt" 2> /dev/null)" != KTPTXT02 ]; then   # (or a .ptxt of the old format)
  stamp "STPD over the parse"
  peak "$W/build_pstpd.log" $CAP $K index --stpd-parse "$P"
  note "index --stpd-parse: $(grep '^parse STPD' "$W/build_pstpd.log" | sed 's/^parse STPD: //'); $(tail -1 "$W/build_pstpd.log")"
fi

# ---- searches (single-threaded)
run() {   # name, options...
  local name=$1; shift
  stamp "$name"
  rm -f "$W/$name.log"
  peak "$W/$name.log" $CAP $RUN $K find --format mems --stats "$@" "$P" "$PAT" "$W/$name.tsv"
  note "  $(printf '%-26s' "$name") $(grep -o ' MEMs=[0-9]*' "$W/$name.log" | head -1 | tr -d ' '), search $(grep -o 'wall [0-9.]* s' "$W/$name.log" | head -1); $(tail -1 "$W/$name.log" | sed 's/, exit 0//')"
  grep -E '^STPD over the parse: right-end|^stpd_jumps|^searches=' "$W/$name.log" | sed 's/^/      /' | tee -a "$SUM" > /dev/null
}
same() {  # name: the same MEMs (name, start, end) as plain search
  if cmp -s <(cut -f1-3 "$W/$1.tsv") <(cut -f1-3 "$W/$2.tsv"); then note "      $1: the same MEMs as $2"; else note "      $1: MEMS DIFFER from $2"; fi
}
first=1
for L in $LS; do
  note "L = $L:"
  run "L$L.plain" -L "$L" --csa
  run "L$L.stpd-parse" -L "$L" --csa --stpd-parse --no-intervals;     same "L$L.stpd-parse" "L$L.plain"
  if [ $first = 1 ]; then
    grep '^STPD over the parse, components' "$W/L$L.stpd-parse.log" | sed 's/^STPD over the parse, components (MB): /  components (MB): /' | tee -a "$SUM"
  fi
  run "L$L.stpd-parse-locate" -L "$L" --csa --stpd-parse --locate;    same "L$L.stpd-parse-locate" "L$L.plain"
  if paste "$W/L$L.plain.tsv" "$W/L$L.stpd-parse-locate.tsv" | awk -F'\t' '$5 - $4 != $11 { bad++ } END { exit bad > 0 }'; then
    note "      L$L.stpd-parse-locate: every occurrence count equals ep - sp"; else note "      L$L.stpd-parse-locate: OCCURRENCE COUNTS DIFFER from ep - sp"; fi
  run "L$L.stpd-oracle" -L "$L" --stpd --parse-oracle --no-intervals; same "L$L.stpd-oracle" "L$L.plain"
  run "L$L.twolevel-blocks" -L "$L" --csa --twolevel;                 same "L$L.twolevel-blocks" "L$L.plain"
  first=0
done
stamp "done: $SUM"
