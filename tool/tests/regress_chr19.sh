#!/bin/bash
# Regression against More RLBWT Tricks (chr19/twolevel): the same MEMs, with the same numbers of
# occurrences, from the old code (forward strand T and the reversed collection T') and from katka find
# (both strands, two-level search with blocks), for BML and forward-backward search.
#   regress_chr19.sh OLDBUILD DATADIR WORKDIR
# OLDBUILD is chr19/twolevel's CMake build directory (with build_twolevel, query_twolevel and
# rs_build/build_index); DATADIR has list.tsv and queries.fa (gen_data.py's chr/).
#
# The old index holds the forward strand only, so it reports a MEM by its occurrences in T, and katka by
# its occurrences on both strands; and a reverse-strand match can make katka's MEMs longer.  On the
# random sequences of gen_data.py, reverse-strand matches of 100 or more bases do not occur, so the
# outputs must agree.  Queries are split at runs of N, which the old code would match against N.
set -euo pipefail
OLD=$(cd "$1" && pwd); DATA=$(cd "$2" && pwd); W=$3
NEW=$(cd "$(dirname "$0")/.." && pwd)/bin/katka
mkdir -p "$W"; cd "$W"
# the collection, one sequence per line, and the patterns (A/C/G/T runs of the forward queries)
python3 - "$DATA" <<'EOF'
import sys, re
d = sys.argv[1]
def fasta(p):
    recs, cur = [], None
    for line in open(p):
        line = line.strip()
        if line.startswith('>'): cur = []; recs.append(cur)
        else: cur.append(line.upper())
    return [''.join(r) for r in recs]
with open('collection.txt', 'w') as f:
    for line in open(d + '/list.tsv'):
        f.write(''.join(fasta(line.rstrip('\n').split('\t')[-1])) + '\n')
pats = []
for q in fasta(d + '/queries.fa')[:2]:
    pats += [p for p in re.split('[^ACGT]+', q) if len(p) >= 100]
open('patterns.txt', 'w').write(''.join(p + '\n' for p in pats))
open('patterns.fa', 'w').write(''.join('>p%d\n%s\n' % (i, p) for i, p in enumerate(pats)))
EOF
ref=$(head -1 collection.txt | tr -d '\n' | wc -c)
"$OLD/rs_build/build_index" --text collection.txt --reference-length $ref --out idx --ell 25 --theta 0 > /dev/null 2>&1
"$OLD/rs_build/build_index" --text collection.txt --reference-length $ref --out idxr --ell 25 --theta 0 --reversed > /dev/null 2>&1
"$OLD/build_twolevel" --text collection.txt --rlbwt idx/rlbwt.bin --reversed-rlbwt idxr/rlbwt.bin --out lv 2> /dev/null
rm -rf dump; mkdir dump
"$OLD/query_twolevel" --index idx --reversed-index idxr --level lv.text.lvl --reversed-level lv.rev.lvl --patterns patterns.txt \
  --algorithm fb,bml --min-length 100,1000 --repeat 1 --configs plain,2LJ --dump dump > /dev/null 2>&1
"$NEW" prep "$DATA/list.tsv" new 2> /dev/null
"$NEW" bwt new 2> /dev/null
"$NEW" index --preset chr19 new 2> /dev/null
fails=0
np=$(grep -c . patterns.txt)
for alg in bml fb; do
  for L in 100 1000; do
    F=""; [ $alg = fb ] && F="--fb"
    "$NEW" find --preset chr19 -L $L $F new patterns.fa new.tsv 2> /dev/null
    for ((q = 0; q < np; ++q)); do
      awk -v p=p$q '$1 == p { print $2 "\t" $3 "\t" $5 - $4 }' new.tsv | sort -n > a.tsv
      sort -n dump/${alg}_L${L}_q$q.tsv > b.tsv
      if ! cmp -s a.tsv b.tsv; then echo "  DIFFERENT: $alg L=$L pattern $q"; fails=$((fails + 1)); fi
    done
    echo "  $alg L=$L: $(wc -l < new.tsv) MEMs in $np patterns"
  done
done
if [ $fails = 0 ]; then echo "regression against chr19/twolevel: passed"; else echo "regression against chr19/twolevel: $fails differences"; exit 1; fi
