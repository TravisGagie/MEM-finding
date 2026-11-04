#!/bin/bash
# The size of a Cliffy index, counted three ways (tool/NOTES.md, "Cliffy's index sizes"):
#   1. as Cliffy's paper counts it: cliffy-experiments' rule combine_all_index_sizes_into_csv_file_exp1 sums
#      .bwt.cliffy, .F.cliffy, .runcnt, .ftab and the six .taxcomp.* profile files (with --taxcomp), or
#      .bwt.cliffy, .F.cliffy, .sdap, .edap, .runcnt (without), and reports GiB (2^30 bytes);
#   2. every file the build wrote (PREFIX.fna*), as our comparison (katka/silva/compare_lib.sh) did with du -cb;
#   3. the files `cliffy run` opens (with strace, if CLIFFY and READS are given).
#
#   cliffy_index_size.sh PREFIX [CLIFFY READS.fq [RUN OPTIONS]]
# PREFIX is the --output of `cliffy build` (the files are PREFIX.fna, PREFIX.fna.bwt.cliffy, ...).  RUN OPTIONS are
# those of `cliffy run` besides --ref, --pattern, --output, --taxcomp and --ftab: default --num-col 7; for the
# digested index, --num-col 5 --minimizers --small-window 4 --large-window 11.
set -euo pipefail
PFX=$1; P=$1.fna
[ -s "$P" ] || { echo "no $P"; exit 1; }
sz() { local t=0; for f in "$@"; do [ -e "$f" ] && t=$((t + $(stat -c %s "$f"))); done; echo $t; }
if [ -e "$P.taxcomp.sdap" ]; then
  theirs=$(sz $P.bwt.cliffy $P.F.cliffy $P.taxcomp.sdap $P.taxcomp.edap $P.taxcomp.of.sdap $P.taxcomp.of.edap \
              $P.taxcomp.ofptr.sdap $P.taxcomp.ofptr.edap $P.runcnt $P.ftab)
  kind="taxonomic compression"
else
  theirs=$(sz $P.bwt.cliffy $P.F.cliffy $P.sdap $P.edap $P.runcnt); kind="full profiles"
fi
all=$(sz $P $P.*)
printf "Cliffy index %s (%s)\n" "$1" "$kind"
printf "  as Cliffy's paper counts it: %d bytes = %.2f GiB = %.2f GB\n" $theirs $(echo "$theirs / 1073741824" | bc -l) $(echo "$theirs / 1000000000" | bc -l)
printf "  every file of the build:     %d bytes = %.2f GB (%.2f times the above)\n" $all $(echo "$all / 1000000000" | bc -l) $(echo "$all / $theirs" | bc -l)
echo "  the files besides the paper's, largest first:"
for f in $P $P.*; do
  case ${f#$P} in .bwt.cliffy|.F.cliffy|.runcnt|.ftab|.taxcomp.sdap|.taxcomp.edap|.taxcomp.of.sdap|.taxcomp.of.edap|.taxcomp.ofptr.sdap|.taxcomp.ofptr.edap) continue;; esac
  [ -e "$P.taxcomp.sdap" ] || case ${f#$P} in .sdap|.edap) continue;; esac
  printf "    %14d  %s\n" $(stat -c %s $f) "${f#$1}"
done | sort -k1,1nr
if [ $# -ge 3 ]; then
  command -v strace > /dev/null || { echo "  (no strace: skipping the files cliffy run opens)"; exit 0; }
  T=$(mktemp -d); trap 'rm -rf $T' EXIT
  CL=$2 RD=$3; shift 3; [ $# -gt 0 ] || set -- --num-col 7
  strace -f -e trace=openat -o $T/trace "$CL" run --ref "$PFX" --pattern "$RD" --output $T/q --taxcomp --ftab "$@" > $T/log 2>&1 ||
    { echo "  cliffy run failed:"; sed 's/\x1b\[[0-9;]*m//g' $T/log | tail -3; exit 1; }
  mapfile -t opened < <(grep -v ENOENT $T/trace | grep -o "\"$P[^\"]*\"" | tr -d '"' | sort -u)
  printf "  opened by cliffy run:        %d bytes = %.2f GB (%s)\n" $(sz "${opened[@]}") $(echo "$(sz "${opened[@]}") / 1000000000" | bc -l) \
    "$(for f in "${opened[@]}"; do printf '%s ' "${f#$PFX}"; done)"
fi
