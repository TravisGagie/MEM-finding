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
#   4. STPD search (stpd.hpp): the same MEMs and intervals as plain search, every search ending at the
#      lexicographically smallest occurrence, and the sample's size against the paper's definition;
#   5. the command line: threads, two-level search with tags, a BWT with one terminator per dataset, the
#      presets;
#   6. if KATKA_OLD (katka-readable's bin/), RZ_OLD (a build of katka/rz-index), CHR19_OLD
#      (chr19/twolevel's build/) or STPD_OLD (a built clone of regindex/STPD-index) is set, the
#      regressions against that code (tests/regress_katka.sh, tests/regress_digest.sh,
#      tests/regress_chr19.sh, tests/regress_stpd.sh).
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
quiet $K index --csa --efcsa --tags --kebab 20 --twolevel --syncmer 21,8 --stpd t
quiet $K index --tags --sample 4 --tag-out t.s4.tag t
quiet $K index --tags --grammar --tag-out t.g.tag t
quiet $K index --tags --by-tag --tag-out t.sp.tag --tag-lcp 0,10,30,100,1000 t
grep -E "rho|lambda" log.txt | sed 's/; r = .*//; s/, bound .*//' > blam.out
$KT blam t 1 0,10,30,100,1000 > blam.ref
same() { if cmp -s "$1" "$2"; then echo "  same: $3"; else echo "FAILED: $3"; exit 1; fi; }
same blam.out blam.ref "B_lambda of the tag runs against brute force"
step "brute force: MEMs, intervals, KeBaB, two-level, tags (s = 1)"
$KT mems t "$D/reads.fq" --tags t.tag --kebab 20 --twolevel --stpd --random 40
$KT stpd-size t
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
quiet $K index --csa --efcsa --twolevel --stpd c
$KT mems c "$D/chr/queries.fa" --twolevel --stpd --random 20
$KT stpd-size c
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
quiet $K find -L 20 -c --csa --stpd t "$D/reads.fq" b.out
same a.out b.out "STPD search with tags (intervals of the reported MEMs by backward search)"
quiet $K find -L 1000 --format mems --csa c "$D/chr/queries.fa" a.out
quiet $K find -L 1000 --format mems --stpd c "$D/chr/queries.fa" b.out
same a.out b.out "STPD search, chromosome-like"
quiet $K find -L 1000 --format mems --stpd --no-intervals c "$D/chr/queries.fa" b.out
cut -f1-3 a.out > a3.out; cut -f1-3 b.out > b3.out
same a3.out b3.out "STPD search without intervals (the same MEMs)"
quiet $K index --parse-oracle c
quiet $K find -L 1000 --format mems --stpd --parse-oracle c "$D/chr/queries.fa" b.out
same a.out b.out "STPD search with the parse oracle (the text read through the two-level parse)"
quiet $K find -L 1000 --format mems --stpd --locate c "$D/chr/queries.fa" d.out
quiet $K find -L 1000 --format mems --stpd --locate --parse-oracle c "$D/chr/queries.fa" e.out
same d.out e.out "STPD locating with the parse oracle"
quiet $K index --stpd-parse c
for L in 1 25 100 1000; do
  quiet $K find -L $L --format mems --csa c "$D/chr/queries.fa" pa.out
  quiet $K find -L $L --format mems --csa --stpd-parse c "$D/chr/queries.fa" pb.out
  same pa.out pb.out "STPD over the parse, L = $L"
  quiet $K find -L $L --format mems --csa --fb c "$D/chr/queries.fa" pa.out
  quiet $K find -L $L --format mems --csa --fb --stpd-parse c "$D/chr/queries.fa" pb.out
  same pa.out pb.out "STPD over the parse, forward-backward, L = $L"
done
quiet $K find -L 100 --format mems --csa c "$D/chr/queries.fa" pa.out
quiet $K find -L 100 --format mems --csa --stpd-parse --locate c "$D/chr/queries.fa" pb.out
paste pa.out pb.out | awk -F'\t' '$2 != $7 || $3 != $8 || $5 - $4 != $11 { bad++ } END { exit bad > 0 }' || { echo "FAILED: STPD over the parse, --locate counts"; exit 1; }
echo "  same: STPD over the parse, --locate (occurrence counts = ep - sp)"
quiet $K index --rlzap c
for L in 1 100 1000; do
  quiet $K find -L $L --format mems --csa c "$D/chr/queries.fa" pa.out
  quiet $K find -L $L --format mems --csa --stpd-parse --rlzap-oracle c "$D/chr/queries.fa" pb.out
  same pa.out pb.out "STPD over the parse through the RLZAP parse oracle, L = $L"
done
quiet $K find -L 100 --format mems --csa --stpd-parse --locate c "$D/chr/queries.fa" pa.out
quiet $K find -L 100 --format mems --csa --stpd-parse --locate --rlzap-oracle c "$D/chr/queries.fa" pb.out
same pa.out pb.out "STPD over the parse through the RLZAP parse oracle, --locate"
quiet $K index --twolevel-locate c
for L in 100 1000; do
  quiet $K find -L $L --format mems --csa c "$D/chr/queries.fa" pa.out
  quiet $K find -L $L --format mems --csa --twolevel --locate c "$D/chr/queries.fa" pb.out
  paste pa.out pb.out | awk -F'\t' '$2 != $7 || $3 != $8 || $4 != $9 || $5 != $10 || $5 - $4 != $11 { bad++ } END { exit bad > 0 }' || { echo "FAILED: two-level locate, L = $L"; exit 1; }
  echo "  same: two-level search with blocks, --locate (the same MEMs and intervals, occurrence counts = ep - sp), L = $L"
done
quiet $K find --preset chr19 c "$D/chr/queries.fa" b.out
same a.out b.out "preset chr19 (two-level, blocks) = plain search"
quiet $K find --preset novel -L 200 t "$D/reads.fq" a.out
quiet $K find -L 200 --format mems --csa t "$D/reads.fq" b.out
same a.out b.out "preset novel (inline KeBaB) = plain BML"
quiet $K find --preset silva t "$D/reads.fq" a.out
quiet $K find -L 30 -c t "$D/reads.fq" b.out
same a.out b.out "preset silva"

step "order-k block feasibility (katka-blocks)"
quiet $K find -L 30 -c --csa t "$D/reads.fq" list30.out
"$HERE/bin/katka-blocks" t "$D/reads.fq" list30.out -k 8,12 --check > blocks.out 2>&1 || { cat blocks.out; echo "FAILED: katka-blocks"; exit 1; }
grep -q "0 reads whose MEMs differ" blocks.out && [ "$(grep -c 'agree on succ for all k-mers: yes' blocks.out)" = 2 ] \
  && echo "  katka-blocks: BML's MEMs equal the listing's, k-mer tables checked, every step classified" \
  || { cat blocks.out; echo "FAILED: katka-blocks"; exit 1; }

step "preset zymo (one BWT string per dataset, species tags)"
(cd "$D" && quiet $K prep --preset zymo list.tsv "$W/z")
quiet $K bwt -d z
quiet $K index --preset zymo z
quiet $K find --preset zymo z "$D/reads.fq" a.out
quiet $K find -U -X -L 20 -K 20 -D 14 -c -T z.tag z "$D/reads.fq" b.out
same a.out b.out "preset zymo"
grep -q "{0:" a.out && echo "  species tags present"
mkdir -p blamtmp
quiet $K index --tags --by-tag --no-rmq --tag-out z.blam.tag --tag-lcp 0,5,20,50 --blam-external --blam-mem 0.000001 --tmp blamtmp z
grep -E "rho|lambda" log.txt | sed 's/; r = .*//; s/, bound .*//' > blam.out
$KT blam z 1 0,5,20,50 > blam.ref
same blam.out blam.ref "B_lambda with one terminator per dataset (boundary records and fingerprints, one value of lambda per pass, temporary files in --tmp) against brute force"
same z.blam.tag z.tag "the tag array built with --tag-lcp"
[ -z "$(ls blamtmp)" ] && echo "  temporary files removed"

if [ -n "${KATKA_OLD:-}" ]; then step "regression against katka-readable"; "$HERE/tests/regress_katka.sh" "$KATKA_OLD" "$D" "$W/regress_katka"; fi
if [ -n "${RZ_OLD:-}" ]; then step "regression against rz-index digests"; "$HERE/tests/regress_digest.sh" "$RZ_OLD" "$D" "$W/regress_digest"; fi
if [ -n "${STPD_OLD:-}" ]; then            # a clone of https://github.com/regindex/STPD-index, built in build/
  step "regression against STPD-index (the paper's st-colex- sampling)"
  "$HERE/tests/regress_stpd.sh" "$STPD_OLD/build/CLI" "$STPD_OLD/data/yeast.txt" "$W/regress_stpd_yeast"
  python3 -c "import sys; print(''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).replace('N', ''))" "$D/chr/h0.fa" > chrT.txt
  "$HERE/tests/regress_stpd.sh" "$STPD_OLD/build/CLI" chrT.txt "$W/regress_stpd_chr"
  step "STPD-index's MEMs (experiments/stpd_index_chr19.sh), on the chromosome-like collection"
  "$HERE/experiments/stpd_index_chr19.sh" "$STPD_OLD" "$D/chr/list.tsv" "$D/chr/queries.fa" "$W/stpd_index_chr" 100 > /dev/null 2>&1
  grep -A9 "^comparison" "$W/stpd_index_chr/summary.txt" | sed 's/^/  /'
  if grep -E "^  (ours|theirs)" "$W/stpd_index_chr/summary.txt" | grep -vq " 0 only here (0 with the first occurrence across a junction), 0 missing"; then echo "FAILED: STPD-index MEMs"; exit 1; fi
  if grep -q "occurrence counts different from ep - sp: [1-9]" "$W/stpd_index_chr/summary.txt"; then echo "FAILED: STPD-index occurrence counts"; exit 1; fi
fi
if [ -n "${CHR19_OLD:-}" ]; then step "regression against chr19/twolevel"; "$HERE/tests/regress_chr19.sh" "$CHR19_OLD" "$D/chr" "$W/regress_chr19"; fi
echo "self-test passed"
