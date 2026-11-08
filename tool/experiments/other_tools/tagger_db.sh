#!/bin/bash
# Tagger's index (https://github.com/biointec/tagger, at deps/tagger) on KATKA's reference, as run_tagger.sh builds it:
# one FASTA whose headers carry the genus as a tag G00000..G09117 (KATKA's document numbers, from 0), its categories
# file, and Tagger's PFP build script (with Big-BWT).  Tagger handles reverse complements at query time, so one strand.
#   tagger_db.sh DIR [LIST.tsv]     LIST.tsv: genus i's FASTA in line i, second column (default ~/tool/silva/list.tsv,
#                                   the full reference; a held-out split's list.tsv for its reduced reference)
# Timed: DIR/build.time (wall s, user s, sys s, peak KB), in an 11 GB cgroup.  The index is DIR/silva*.
set -euo pipefail
DIR=$1; LIST=${2:-$HOME/tool/silva/list.tsv}; TG=${TG:-$HOME/rz/deps/tagger}
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
rm -rf "$DIR"; mkdir -p "$DIR"; DIR=$(cd "$DIR" && pwd)
i=0; cut -f2 "$LIST" | while read -r f; do
  awk -v t=$(printf 'G%05d' $i) '/^>/ { print ">" t "_s" ++k; next } { print }' "$f"; i=$((i + 1)); done > "$DIR/ref_tags.fa"
n=$(wc -l < "$LIST"); for i in $(seq 0 $((n - 1))); do printf 'G%05d\n' $i; done > "$DIR/categories.txt"
echo "$(grep -c '^>' $DIR/ref_tags.fa) sequences, $n genera"
(cd "$TG/build" && $CAP /usr/bin/time -f "%e %U %S %M" -o "$DIR/build.time" bash ./tagger_build_pfp.sh -r "$DIR/silva" -t "$DIR/categories.txt" -f "$DIR/ref_tags.fa" > "$DIR/build.log" 2>&1)
rm -f "$DIR/ref_tags.fa"
echo "build: $(cat $DIR/build.time) (wall s, user s, sys s, peak KB); $(du -cb $DIR/silva* | tail -1 | cut -f1) bytes"
