#!/bin/bash
# Grammar-compressed tags: the built-in RePair against bigrepair, on one index (tool/NOTES.md, "Grammar tags").
#
#   grammar_tags.sh PREFIX READS WORK [FIND OPTIONS]
#
# PREFIX: a katka index with tags (index --preset silva, or index --tags): PREFIX.tag and the tag-run cache
# PREFIX.tagruns (PREFIX.sp.tagruns with BYTAG=1, i.e. --by-tag).  READS: reads to answer.  FIND OPTIONS
# (default "-L 30 --csa -F 10", the silva preset without its answer mode) are given to every search.
# Environment: BIGREPAIR (path to bigrepair, from https://gitlab.com/manzai/bigrepair; required), THREADS
# (searches, default 1), MEMMAX (an 11G cgroup per step if systemd-run works; default 11G), BGMEM (bigrepair's
# -m: default 0, its plain RePair.  Its memory-bounded RePair, which it uses for any -m > 0 and by default with
# 95% of the machine's RAM, stalled here on 64 million tag runs (over 30 min against 57 s for the plain one, which
# peaked at 0.83 GB), though it took 30 s on either half of them),
# BGTIME (a time limit for bigrepair, default 8h: if it does not finish, the rest runs without its tag array).
#
# Steps, each with its time and peak memory in WORK/summary.txt:
#   1. index --tags --grammar --tag-out WORK/g.tag         the built-in approximate RePair
#   2. scripts/tagruns_to_int32.py, then bigrepair -i -m     bigrepair's grammar of the same tag sequence
#   3. index --tags --bigrepair WORK/tags.int32 --tag-out WORK/bg.tag
#   4. find with PREFIX.tag, WORK/g.tag and WORK/bg.tag, counting (-c) and listing (-l): the outputs must be
#      byte-identical (only the tag array's encoding differs)
set -euo pipefail
P=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); R=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir -p "$3"; W=$(cd "$3" && pwd); shift 3
FIND=${*:-"-L 30 --csa -F 10"}
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
[ -n "${BIGREPAIR:-}" ] && [ -x "$BIGREPAIR" ] || { echo "set BIGREPAIR to bigrepair (git clone https://gitlab.com/manzai/bigrepair && make -C bigrepair)"; exit 1; }
BT=""; CACHE=$P.tagruns; [ -n "${BYTAG:-}" ] && { BT=--by-tag; CACHE=$P.sp.tagruns; }
CAP=""; systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0 true 2> /dev/null && CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
SUM=$W/summary.txt; : > "$SUM"
note() { echo "$*" | tee -a "$SUM"; }
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
last() { tail -1 "$1" | sed 's/, exit 0//'; }
[ -s "$P.tag" ] || { echo "no $P.tag (katka index --preset silva $P, or --tags)"; exit 1; }
[ -s "$CACHE" ] || { echo "no tag-run cache $CACHE (it is written by katka index --tags$BT)"; exit 1; }
note "grammar tags on $P; reads $R; find $FIND; katka $(git -C "$HERE" log --oneline -1 | cut -c1-60)"
note "tag runs: $(python3 -c "import struct,sys; f=open(sys.argv[1],'rb'); f.read(8); n,k=struct.unpack('<QQ', f.read(16)); print('%d runs over %d rows' % (k, n))" "$CACHE")"

rm -f "$W"/*.log
peak "$W/g.log" $CAP $K index --tags $BT --grammar --tag-out "$W/g.tag" "$P"
note "1. built-in RePair: $(grep -o 'tag array: .*' "$W/g.log" | sed 's/ (.*//'); $(last "$W/g.log")"
peak "$W/int32.log" python3 "$HERE/scripts/tagruns_to_int32.py" "$CACHE" "$W/tags.int32"
BG=yes
if ! peak "$W/bigrepair.log" timeout ${BGTIME:-8h} $CAP "$BIGREPAIR" -i -r -m ${BGMEM:-0} "$W/tags.int32"; then
  BG=no; note "2. bigrepair -i: DID NOT FINISH within ${BGTIME:-8h} (BGTIME), or failed; $(last "$W/bigrepair.log"); see $W/bigrepair.log and $W/tags.int32.log"
else
  note "2. bigrepair -i: .R $(stat -c %s "$W/tags.int32.R") bytes, .C $(stat -c %s "$W/tags.int32.C") bytes; $(last "$W/bigrepair.log")"
  peak "$W/bg.log" $CAP $K index --tags $BT --bigrepair "$W/tags.int32" --tag-out "$W/bg.tag" "$P"
  note "3. from bigrepair's grammar: $(grep -o 'tag array: .*' "$W/bg.log" | sed 's/ (.*//'); $(last "$W/bg.log")"
fi
TAGS="plain g"; [ $BG = yes ] && TAGS="plain g bg"
note "   tag array files: plain $(stat -L -c %s "$P.tag") bytes (+ RMQ $( [ -e "$P.tag.rmq" ] && stat -L -c %s "$P.tag.rmq" || echo 0)), built-in RePair $(stat -c %s "$W/g.tag"), bigrepair $( [ $BG = yes ] && stat -c %s "$W/bg.tag" || echo none)"
for m in c l; do
  for t in $TAGS; do
    tf=$W/$t.tag; [ $t = plain ] && tf=$P.tag
    peak "$W/find.$m.$t.log" $K find $FIND -$m -T "$tf" -j ${THREADS:-1} "$P" "$R" "$W/out.$m.$t"
  done
  same=yes; for t in ${TAGS#plain }; do cmp -s "$W/out.$m.plain" "$W/out.$m.$t" || same=NO; done
  note "4. find -$m: outputs byte-identical with the tag arrays ($TAGS): $same"
  for t in $TAGS; do note "     $t: $(grep -o 'timing: .*' "$W/find.$m.$t.log"); $(last "$W/find.$m.$t.log")"; done
  [ $same = yes ] || { echo "OUTPUTS DIFFER (find -$m)"; exit 1; }
done
note "done"
