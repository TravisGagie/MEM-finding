#!/bin/bash
# A Kraken 2 database on SILVA 138.1 NR99, built as kraken2-build --special silva does (its taxonomy script, the
# accession-to-taxon map, U to T), but from local files, so that the library can be chosen:
#   kraken2_db.sh DB silva [EXCLUDE]     all 510,508 NR99 sequences (Kraken 2's standard SILVA build, the one Cliffy used)
#   kraken2_db.sh DB ours [EXCLUDE]      the 402,378 sequences KATKA indexes (work/ref/dna/doc_*_seq.fa)
#   kraken2_db.sh DB LIST.tsv [EXCLUDE]  the files in LIST.tsv's second column (e.g. a held-out split's list.tsv)
# EXCLUDE: a file of sequence IDs (or a FASTA whose headers give them) to leave out.  Only `kraken2-build --build` is
# timed (DB/build.time: wall s, user s, sys s, peak KB), with THREADS threads (default 12), in an 11 GB cgroup.
# K2DATA: SILVA's SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz, tax_slv_ssu_138.1.acc_taxid.gz, tax_slv_ssu_138.1.txt.gz.
set -euo pipefail
DB=$1; SRC=$2; EX=${3:-}; T=${THREADS:-12}; D=${K2DATA:-$HOME/tool/silva/k2/data}
S=$(cd "$(dirname "$0")/../../../katka/silva" && pwd)
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
rm -rf "$DB"; mkdir -p "$DB/taxonomy" "$DB/library"; DB=$(cd "$DB" && pwd)
(cd "$DB" && zcat "$D/tax_slv_ssu_138.1.txt.gz" > tax.txt && perl /usr/lib/kraken2/build_silva_taxonomy.pl tax.txt && mv names.dmp nodes.dmp taxonomy/ && rm tax.txt)
zcat "$D/tax_slv_ssu_138.1.acc_taxid.gz" > "$DB/seqid2taxid.map"
if [ -n "$EX" ]; then { grep '^>' "$EX" | sed 's/^>//' || cat "$EX"; } | awk '{print $1}' > "$DB/excluded.ids"; else : > "$DB/excluded.ids"; fi
case $SRC in
  silva) zcat "$D/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" | sed -e '/^>/!y/U/T/' ;;
  ours)  for i in $(seq 1 "$(cut -d' ' -f1 $S/work/ref/done)"); do cat $S/work/ref/dna/doc_${i}_seq.fa; done ;;
  *.tsv) cut -f2 "$SRC" | while read -r f; do cat "$f"; done ;;
  *) echo "unknown source $SRC"; exit 1 ;;
esac | awk -v ex="$DB/excluded.ids" 'BEGIN { while ((getline l < ex) > 0) skip[l] = 1 }
  /^>/ { split(substr($0, 2), a, " "); keep = !(a[1] in skip); if (keep) n++ } keep { print } END { print n " sequences" > "/dev/stderr" }' \
  > "$DB/library/silva.fna" 2> "$DB/library.count"
grep -v '^>' "$DB/library/silva.fna" | tr -d '\n' | wc -c | awk '{print $1 " bases"}' >> "$DB/library.count"
cat "$DB/library.count"
$CAP /usr/bin/time -f "%e %U %S %M" -o "$DB/build.time" kraken2-build --db "$DB" --build --threads $T > "$DB/build.log" 2>&1
rm -f "$DB/library/silva.fna"
echo "build: $(cat $DB/build.time) (wall s, user s, sys s, peak KB); $(du -cb $DB/*.k2d | tail -1 | cut -f1) bytes of .k2d"
