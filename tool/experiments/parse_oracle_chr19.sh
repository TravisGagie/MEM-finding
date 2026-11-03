#!/bin/bash
# STPD search with the text read through the two-level parse (katka index/find --parse-oracle) against the
# plain text (PREFIX.S in memory): the same MEMs, times, sizes and peak memory.
#
#   parse_oracle_chr19.sh PREFIX patterns.fa WORK [k,s ...]
# PREFIX: an existing katka index with PREFIX.S and PREFIX.stpd (katka index --stpd); patterns.fa: the
# queries; k,s: the syncmer settings of the parses to try (default 51,12, the chr19 preset's).  Each parse
# oracle is built in WORK/k<k>s<s>/c.ptxt next to symbolic links to PREFIX's files (.S, .stpd, .rlbwt, ...),
# so PREFIX is not changed.  Environment: PIN (a CPU to pin the searches to), LS (the thresholds, default "1000 100").
#
# Runs, for each L: --stpd --no-intervals and --stpd --locate, each with the text and with every parse
# oracle; at the first L also --stpd --fb --no-intervals.  (With intervals, STPD search also runs backward
# search on the BWT, which does not use the oracle.)  Output: WORK/summary.txt.
set -euo pipefail
P=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); PAT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir -p "$3"; W=$(cd "$3" && pwd); shift 3
SETS=("$@"); [ ${#SETS[@]} -gt 0 ] || SETS=(51,12)
LS=${LS:-"1000 100"}
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
for f in "$P.S" "$P.stpd" "$P.rlbwt"; do [ -s "$f" ] || { echo "missing $f (katka index --stpd $1)"; exit 1; }; done
RUN=""; [ -n "${PIN:-}" ] && command -v taskset > /dev/null && RUN="taskset -c $PIN"
SUM=$W/summary.txt
note() { echo "$*" | tee -a "$SUM"; }
stamp() { echo "[$(date '+%F %T')] $*"; }
peak() {   # runs a command (stdout and stderr to the log named first), appending wall time and peak memory
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
note "STPD search: plain text against parse oracles; katka $(git -C "$HERE" log --oneline -1 2>/dev/null || echo '?'); index $P; patterns $PAT; L in $LS"
note "text: $(ls -lL "$P.S" | awk '{printf "%.1f MB", $5 / 1e6}') (held in memory); STPD sample: $(ls -lL "$P.stpd" | awk '{printf "%.1f MB", $5 / 1e6}')"

# ---- the parse oracles
for ks in "${SETS[@]}"; do
  d=$W/k${ks/,/s}; mkdir -p "$d"
  for f in "$P".*; do
    x=${f#"$P"}; case $x in .ptxt|.ptxt.tmp|.bwt) continue;; esac
    ln -sf "$f" "$d/c$x"
  done
  if [ ! -s "$d/c.ptxt" ]; then
    stamp "parse oracle k,s = $ks"
    rm -f "$d/build.log"
    peak "$d/build.log" $K index --parse-oracle --syncmer "$ks" "$d/c"
  fi
  note "oracle k,s = $ks: $(grep '^parse oracle: n=' "$d/build.log" | sed 's/^parse oracle: //'); build $(tail -1 "$d/build.log" | sed 's/, exit 0//')"
done

# ---- searches (single-threaded)
run() {   # name, prefix, options...
  local name=$1 pre=$2; shift 2
  stamp "$name"
  rm -f "$W/$name.log"
  peak "$W/$name.log" $RUN $K find --format mems --stats "$@" "$pre" "$PAT" "$W/$name.tsv"
  local seeks; seeks=$(grep '^parse oracle seeks' "$W/$name.log" | sed 's/parse oracle seeks: //' || true)
  note "  $(printf '%-34s' "$name") $(grep -o ' MEMs=[0-9]*' "$W/$name.log" | head -1 | tr -d ' '), search $(grep -o 'wall [0-9.]* s' "$W/$name.log" | head -1); $(grep -o 'stpd_jumps=[0-9]* stpd_compares=[0-9]*' "$W/$name.log"); $(tail -1 "$W/$name.log" | sed 's/, exit 0//')${seeks:+; seeks: $seeks}"
}
first=1
for L in $LS; do
  note "L = $L:"
  modes=("nointervals:--no-intervals" "locate:--locate")
  [ $first = 1 ] && modes+=("fb:--fb --no-intervals")
  for md in "${modes[@]}"; do
    m=${md%%:*}; opts=${md#*:}
    run "L$L.$m.text" "$P" -L "$L" --stpd $opts
    for ks in "${SETS[@]}"; do
      t=k${ks/,/s}
      run "L$L.$m.$t" "$W/$t/c" -L "$L" --stpd --parse-oracle $opts
      if cmp -s "$W/L$L.$m.text.tsv" "$W/L$L.$m.$t.tsv"; then note "    $t: the same output as the text"; else note "    $t: OUTPUT DIFFERS from the text"; fi
    done
  done
  first=0
done
stamp "done: $SUM"
