#!/bin/bash
# The published STPD-index (github.com/regindex/STPD-index) against katka on the same MEM question:
# both strands, threshold L (default 1000), the same patterns.
#
#   stpd_index_chr19.sh STPDIR list.tsv patterns.fa WORK [L]
# STPDIR: a clone of STPD-index, built (see NOTES.md for the cmake line); list.tsv: the collection, one
# FASTA per sequence (as for katka prep); patterns.fa: the queries (FASTA; split at runs of N by both
# tools).  Environment: BIGBWT (path to Big-BWT's bigbwt; without it, katka bwt, which needs about 6 bytes
# per character), OURS (the prefix of an existing katka index of the same collection, built with
# `katka index --preset chr19 --stpd`; without it, one is built in WORK), PIN (a CPU to pin queries to).
#
# Their index: their stpd-build computes its st-colex- sampling with a suffix tree (stpd-sampling, about
# 68 bytes per character: far too much for whole chromosome collections), so the sampling and the arrays
# it needs come from katka instead, in O(r) memory plus the text (katka stpd-prep / index --stpd /
# stpd-export: byte-identical to stpd-sampling -P on the texts tested); stpd-build then finds them and
# builds its index (RLZ text oracle with the first sequence, both strands, as the reference; r-index phi).
# Their text has no separators (it must be over A/C/G/T only): the sequences and their reverse complements,
# split at runs of N, are concatenated, so a match can run across a junction; such MEMs are flagged.
#
# Their tools locate whole patterns; experiments/stpd_index_mems.cpp runs BML and forward-backward search
# on their longest-prefix matching (see its header), and their locate_all_occurrences for every MEM.
#
# Output in WORK: their MEMs and ours (TSV), summary.txt with times, sizes and the comparison.
set -euo pipefail
STPD=$(cd "$1" && pwd); LIST=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); PAT=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
mkdir -p "$4"; W=$(cd "$4" && pwd); L=${5:-1000}   # (absolute: the script cds into it)
HERE=$(cd "$(dirname "$0")/.." && pwd)
K=$HERE/bin/katka; DRV=$HERE/bin/stpd-index-mems
mkdir -p "$W"; cd "$W"
SUM=$W/summary.txt
note() { echo "$*" | tee -a "$SUM"; }
stamp() { echo "[$(date '+%F %T')] $*"; }
RUN=""; [ -n "${PIN:-}" ] && command -v taskset > /dev/null && RUN="taskset -c $PIN"
peak() {   # runs a command, appending its time and peak memory to the log named first
  local log=$1; shift
  python3 - "$log" "$@" <<'PY'
import os, resource, subprocess, sys, time
log, cmd = sys.argv[1], sys.argv[2:]
t = time.time()
with open(log, 'a') as f:
    f.flush()
    r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    f.write('wall %.1f s, peak memory %.2f GB, exit %d: %s\n' % (time.time() - t, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024 / 1e9,
            r.returncode, ' '.join(os.path.basename(c) if i == 0 else c for i, c in enumerate(cmd[:3]))))
sys.exit(r.returncode)
PY
}
: > "$SUM"
note "STPD-index ($(git -C "$STPD" log --oneline -1 2>/dev/null || echo '?')) against katka ($(git -C "$HERE" log --oneline -1 2>/dev/null || echo '?')); L = $L; patterns $PAT"

[ -x "$DRV" ] || "$HERE/experiments/build_stpd_index_mems.sh" "$STPD" "$DRV"

# ---- their index, from katka's sampling and arrays
if [ ! -s st.stpdi ]; then
  stamp "their text T (both strands, A/C/G/T pieces) and S = reverse(T)"
  peak build_st.log "$K" stpd-prep "$LIST" st
  stamp "BWT of S"
  if [ -n "${BIGBWT:-}" ]; then peak build_st.log "$BIGBWT" st.S && mv st.S.bwt st.bwt; else peak build_st.log "$K" bwt st; fi
  stamp "katka's STPD sample of S (= their st-colex- sampling of T)"
  peak build_st.log "$K" index --stpd st
  stamp "their input files: st.T.{colex_m,rbwt,pa,lcs}"
  peak build_st.log "$K" stpd-export st
  stamp "their stpd-build (RLZ reference: $(cat st.reflen) characters)"
  peak build_st.log "$STPD/build/CLI/stpd-build" -i st.T -o st.stpdi -l "$(cat st.reflen)"
  grep -q "Suffix Tree path decomposition found" build_st.log || { echo "stpd-build did not use katka's files"; exit 1; }
  rm -f st.bwt st.rlbwt st.stpd st.S          # (st.T is their text oracle's source; stpd-build removed its inputs)
fi
note "their index: $(ls -l st.stpdi | awk '{printf "%.1f MiB", $5 / 1048576}'); build steps (wall time, peak memory):"
grep -E "^wall" build_st.log | sed 's/^/  /' | tee -a "$SUM"

# ---- katka's index of the same collection (with separators)
if [ -z "${OURS:-}" ]; then
  OURS=$W/ours
  if [ ! -s ours.lvl ]; then
    stamp "katka's index"
    peak build_ours.log "$K" prep "$LIST" ours
    if [ -n "${BIGBWT:-}" ]; then peak build_ours.log "$BIGBWT" ours.S && mv ours.S.bwt ours.bwt; else peak build_ours.log "$K" bwt ours; fi
    peak build_ours.log "$K" index --preset chr19 --stpd ours
  fi
fi

# ---- runs (single-threaded)
theirs() {   # name, options
  local name=$1; shift
  stamp "theirs: $name"
  $RUN "$DRV" st.stpdi "$PAT" "$L" "$@" --pieces st.pieces > "theirs.$name.tsv" 2> "theirs.$name.log"
  note "theirs $name: $(grep '^stpd-index-mems' "theirs.$name.log" | sed 's/^stpd-index-mems: //')"
}
ours() {     # name, options
  local name=$1; shift
  stamp "ours: $name"
  $RUN "$K" find -L "$L" --format mems --stats "$@" "$OURS" "$PAT" "ours.$name.tsv" 2> "ours.$name.log"
  note "ours $name: $(grep -o 'reads=[0-9]* MEMs=[0-9]*' "ours.$name.log"), $(grep -o 'wall [0-9.]* s' "ours.$name.log"); $(grep searches "ours.$name.log")"
}
theirs bml
theirs bml-locate --locate
theirs fb --fb
ours plain --csa
ours twolevel-blocks --csa --twolevel
ours stpd-nointervals --stpd --no-intervals
ours stpd --stpd
ours stpd-twolevel --csa --stpd --twolevel
ours stpd-locate --stpd --locate
ours fb-twolevel-blocks --csa --twolevel --fb

# ---- comparison
python3 - "$SUM" <<'PY' | tee -a "$SUM"
import sys
def mems(f, occ_col=None):
    d = {}
    for line in open(f):
        x = line.rstrip('\n').split('\t')
        d[(x[0], int(x[1]), int(x[2]))] = (int(x[occ_col]) if occ_col is not None and x[occ_col] not in ('-',) else None, x[-1] == 'junction')
    return d
ref = mems('ours.plain.tsv')
counts = {k: None for k in ref}
for line in open('ours.plain.tsv'):
    x = line.split('\t'); counts[(x[0], int(x[1]), int(x[2]))] = int(x[4]) - int(x[3])
print('comparison (reference: katka plain BML, %d MEMs):' % len(ref))
for name, occ in [('ours.twolevel-blocks', None), ('ours.stpd-nointervals', None), ('ours.stpd', None), ('ours.stpd-twolevel', None),
                  ('ours.stpd-locate', 5), ('theirs.bml', None), ('theirs.bml-locate', 3)]:
    d = mems(name + '.tsv', occ)
    only_t = [k for k in d if k not in ref]; only_r = [k for k in ref if k not in d]
    junc = sum(1 for k in only_t if d[k][1])
    line = '  %-24s %7d MEMs, %d only here (%d with the first occurrence across a junction), %d missing' % (name, len(d), len(only_t), junc, len(only_r))
    if occ is not None:
        diff = sum(1 for k in d if k in counts and d[k][0] != counts[k])
        line += '; occurrence counts different from ep - sp: %d' % diff
    print(line)
fbr = mems('ours.fb-twolevel-blocks.tsv'); fbt = mems('theirs.fb.tsv')
print('  forward-backward: ours %d MEMs, theirs %d, %d in both' % (len(fbr), len(fbt), len(set(fbr) & set(fbt))))
PY
stamp "done: $SUM"
