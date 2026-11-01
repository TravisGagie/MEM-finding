#!/bin/bash
# Builds the tool and checks it on small synthetic data (tests/gen_data.py):
#   1. unit tests of the bit structures, the RMQ and the grammar;
#   2. brute force (tests/test.cpp, against a suffix array): MEMs and their intervals for BML and
#      forward-backward search with L = 1, 12, 25, 40, adaptive L, the t longest MEMs, KeBaB pseudo-MEMs
#      and the inline KeBaB window check, on reads and on random patterns of both strands; every search
#      configuration (RLBWT, RLCSAs, lookup table, two-level with and without block jumps) giving the
#      same output; the tag answers (counting and listing) with every kind of tag array; and two-level
#      longest-suffix searches against plain backward search on both strands;
#   3. the same with homopolymer compression; BML over minimizer digests (spans against brute force);
#      and on a chromosome-19-like collection with both continuation rules;
#   4. the command line: threads, two-level search with tags, a BWT with one terminator per dataset, the
#      presets;
#   5. if KATKA_OLD (katka-readable's bin/), RZ_OLD (a build of katka/rz-index) or CHR19_OLD
#      (chr19/twolevel's build/) is set, the regressions against the old code (tests/regress_katka.sh,
#      tests/regress_digest.sh, tests/regress_chr19.sh).
# Everything is written to selftest/ (ignored by git).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
make -C "$HERE" -j4 > /dev/null
K=$HERE/bin/katka
KT=$HERE/bin/katka-test
W=$HERE/selftest
rm -rf "$W"; mkdir -p "$W"; cd "$W"
step() { echo "== $*"; }
quiet() { "$@" > log.txt 2>&1 || { cat log.txt; echo "FAILED: $*"; exit 1; }; }

step "unit tests"
$KT unit

step "data"
python3 "$HERE/tests/gen_data.py" "$W/data" 1
D=$W/data

step "index of both strands (20 genomes, 4 species)"
(cd "$D" && quiet $K prep -s list.tsv "$W/t")
quiet $K bwt t
quiet $K index --csa --efcsa --tags --kebab 20 --twolevel --syncmer 21,8 t
quiet $K index --tags --sample 4 --tag-out t.s4.tag t
quiet $K index --tags --grammar --tag-out t.g.tag t
quiet $K index --tags --by-tag --tag-out t.sp.tag t
step "brute force: MEMs, intervals, KeBaB, two-level, tags (s = 1)"
$KT mems t "$D/reads.fq" --tags t.tag --kebab 20 --twolevel --random 40
for v in s4 g sp; do
  step "brute force: tags ($v)"
  $KT mems t "$D/reads.fa" --tags t.$v.tag
done

step "homopolymer compression"
(cd "$D" && quiet $K prep -s -U list.tsv "$W/h")
quiet $K bwt h
quiet $K index --csa --efcsa --tags --kebab 20 --twolevel --syncmer 21,8 h
$KT mems h "$D/reads.fq" --tags h.tag --kebab 20 --twolevel --hpc

step "minimizer digests (k = 4, w = 11)"
(cd "$D" && quiet $K prep -s --digest 4,11 list.tsv "$W/g")
quiet $K bwt g
quiet $K index --csa --efcsa --tags g
$KT mems g "$D/reads.fq" --digest --tags g.tag

step "chromosome-19-like collection (k = 51, s = 12)"
(cd "$D/chr" && quiet $K prep list.tsv "$W/c")
quiet $K bwt c
quiet $K index --csa --efcsa --twolevel c
$KT mems c "$D/chr/queries.fa" --twolevel --random 20
for f in c.*; do ln -sf "$f" "cf${f#c}"; done
rm cf.lvl; quiet $K index --csa --efcsa --twolevel --continuation frequent cf
step "  with the most frequent successor"
$KT mems cf "$D/chr/queries.fa" --twolevel --random 20

step "command line"
same() { if cmp -s "$1" "$2"; then echo "  same: $3"; else echo "FAILED: $3"; exit 1; fi; }
quiet $K find -L 20 -c --csa t "$D/reads.fq" a.out
quiet $K find -L 20 -c --csa -j 4 t "$D/reads.fq" b.out
same a.out b.out "4 threads"
quiet $K find -L 20 -c --csa --twolevel t "$D/reads.fq" b.out
same a.out b.out "two-level search with tags"
quiet $K find -L 20 -c --csa --twolevel --no-jumps -T t.g.tag t "$D/reads.fq" b.out
same a.out b.out "two-level search without jumps, grammar-compressed tags"
quiet $K find -L 20 -c --efcsa -F 8 -T t.s4.tag t "$D/reads.fa" b.out
same a.out b.out "Elias-Fano RLCSA, lookup table, sampled tags, FASTA"
gzip -c "$D/reads.fq" > reads.fq.gz
quiet $K find -L 20 -c --csa t reads.fq.gz b.out
same a.out b.out "gzipped reads"
quiet $K find -L 20 -c --csa --fb t "$D/reads.fq" b.out
same a.out b.out "forward-backward search (all MEMs of at least L)"
quiet $K find -L 1000 --format mems --csa c "$D/chr/queries.fa" a.out
quiet $K find --preset chr19 c "$D/chr/queries.fa" b.out
same a.out b.out "preset chr19 (two-level, blocks) = plain search"
quiet $K find --preset novel -L 200 t "$D/reads.fq" a.out
quiet $K find -L 200 --format mems --csa t "$D/reads.fq" b.out
same a.out b.out "preset novel (inline KeBaB) = plain BML"
quiet $K find --preset silva t "$D/reads.fq" a.out
quiet $K find -L 30 -c t "$D/reads.fq" b.out
same a.out b.out "preset silva"

step "preset zymo (one BWT string per dataset, species tags)"
(cd "$D" && quiet $K prep --preset zymo list.tsv "$W/z")
quiet $K bwt -d z
quiet $K index --preset zymo z
quiet $K find --preset zymo z "$D/reads.fq" a.out
quiet $K find -U -X -L 20 -K 20 -D 14 -c -T z.tag z "$D/reads.fq" b.out
same a.out b.out "preset zymo"
grep -q "{0:" a.out && echo "  species tags present"

if [ -n "${KATKA_OLD:-}" ]; then step "regression against katka-readable"; "$HERE/tests/regress_katka.sh" "$KATKA_OLD" "$D" "$W/regress_katka"; fi
if [ -n "${RZ_OLD:-}" ]; then step "regression against rz-index digests"; "$HERE/tests/regress_digest.sh" "$RZ_OLD" "$D" "$W/regress_digest"; fi
if [ -n "${CHR19_OLD:-}" ]; then step "regression against chr19/twolevel"; "$HERE/tests/regress_chr19.sh" "$CHR19_OLD" "$D/chr" "$W/regress_chr19"; fi
echo "self-test passed"
