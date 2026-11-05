#!/bin/bash
# Block length as genomes are added, under the smallest and the most frequent successor (tool/NOTES.md, "Block
# length as genomes are added").
#
#   block_scaling.sh DATA WORK "n1 n2 ..." [L]
#
# DATA holds collection.txt (one sequence per line, the reference first) and queries.txt (held-out sequences, one
# per line): scripts/pangenome_sim.py's output, or chr19's (experiments/hprc_chr19.py; CHM13 first, then both
# haplotypes of each individual).  n: the number of sequences of each collection, the first n lines (chr19: 33,
# 65 and 129 for 16, 32 and 64 individuals; synthetic: 1 + N for N haplotypes).  L: the MEM length (default 1000).
#
# For each n, under WORK/c<n>: prep, the BWT (Big-BWT with BIGBWT=..., else katka bwt), index --csa --twolevel
# (the smallest successor) -- or, with EXISTING=<prefix pattern, %d for n> (e.g. ~/tool/chr19/c%d), links to an
# index built already --; then the two-level index with the most frequent successor under WORK/c<n>f (links to
# the same text, BWT and RLCSA; only PREFIX.lvl is rebuilt).  For each rule, `find -L L --format mems --csa
# --twolevel --stats` on the queries, with and without block steps (--no-jumps), single-threaded, pinned to PIN
# (default 0) if taskset exists; all outputs must be identical, and, with CHECK_PLAIN=1, equal to plain search's.
# MEMMAX (11G) caps each step if systemd-run works.
#
# Writes WORK/blocks.tsv and prints the table: per n and rule, the parse (phrases, distinct), the blocks (number,
# phrases per block, the fraction of phrase positions that are continuations, symbols, block BWT runs), the
# two-level index (MB: .lvl, its parse and block parts), and for the queries: character, parse and block steps
# per pattern character, the phrases crossed per block step (the mean block length as the patterns see it:
# parse steps without blocks minus parse steps with blocks, over block steps), and the search times with and
# without block steps.
set -euo pipefail
DATA=$(cd "$1" && pwd); mkdir -p "$2"; W=$(cd "$2" && pwd); NS=$3; L=${4:-1000}
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
CAP=""; systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0 true 2> /dev/null && CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
TS=""; command -v taskset > /dev/null && TS="taskset -c ${PIN:-0}"
stamp() { echo "[$(date '+%F %T')] $*" >&2; }
Q=$W/queries.fa
[ -s "$Q" ] || tr -c 'ACGT\n' '\n' < "$DATA/queries.txt" | awk 'length($0) > 0 { print ">p" NR - 1; print }' > "$Q"
QLEN=$(grep -v '^>' "$Q" | tr -d '\n' | wc -c)
for n in $NS; do
  P=$W/c$n; F=$W/c${n}f
  if [ -n "${EXISTING:-}" ]; then
    E=$(printf "$EXISTING" "$n")
    for x in S bwt rlbwt csa info tbl lvl; do [ -e "$E.$x" ] && ln -sf "$E.$x" "$P.$x"; done
    [ -s "$P.lvl" ] || { echo "no $E.lvl"; exit 1; }
  fi
  if [ ! -s "$P.lvl" ]; then
    stamp "c$n: prep, BWT, index --csa --twolevel (smallest successor)"
    mkdir -p "$W/fa$n"
    head -n "$n" "$DATA/collection.txt" | awk -v d="$W/fa$n" '{ f = sprintf("%s/s%04d.fa", d, NR); print ">s" NR > f; print > f; close(f) }'
    for f in "$W/fa$n"/s*.fa; do printf 'seq\t%s\n' "$f"; done > "$W/list$n.tsv"
    $CAP $K prep "$W/list$n.tsv" "$P" > "$W/build$n.log" 2>&1; rm -rf "$W/fa$n"
    if [ -n "${BIGBWT:-}" ]; then $CAP "$BIGBWT" "$P.S" >> "$W/build$n.log" 2>&1; mv "$P.S.bwt" "$P.bwt"; rm -f "$P.S.log"
    else $CAP $K bwt "$P" >> "$W/build$n.log" 2>&1; fi
    $CAP $K index --csa --twolevel "$P" > "$W/index$n.log" 2>&1
  fi
  if [ ! -s "$F.lvl" ]; then
    stamp "c$n: the two-level index with the most frequent successor"
    for x in S bwt rlbwt csa info tbl; do [ -e "$P.$x" ] && ln -sf "$(readlink -f "$P.$x")" "$F.$x"; done
    $CAP $K index --twolevel --continuation frequent "$F" > "$W/index${n}f.log" 2>&1
  fi
  if [ -n "${EXISTING:-}" ] && [ ! -s "$W/index$n.log" ]; then
    # the statistics of an existing smallest-successor index: rebuild its .lvl under a probe prefix (the same file)
    stamp "c$n: the existing index's two-level statistics (rebuilding its .lvl under a probe prefix)"
    G=$W/c${n}s; for x in S bwt rlbwt csa info tbl; do [ -e "$P.$x" ] && ln -sf "$(readlink -f "$P.$x")" "$G.$x"; done
    $CAP $K index --twolevel "$G" > "$W/index$n.log" 2>&1; rm -f "$G".*
  fi
  for rule in s f; do
    X=$P; [ $rule = f ] && X=$F
    for mode in blocks nojumps; do
      o=$W/q$n$rule.$mode
      [ -s "$o.log" ] && continue
      stamp "c$n, rule $rule: find ($mode)"
      opt=""; [ $mode = nojumps ] && opt="--no-jumps"
      $CAP $TS $K find -L $L --format mems --csa --twolevel $opt --stats -j 1 "$X" "$Q" "$o.tsv" > "$o.log" 2>&1
    done
  done
  for o in "$W/q${n}s.nojumps" "$W/q${n}f.blocks" "$W/q${n}f.nojumps"; do
    cmp -s "$W/q${n}s.blocks.tsv" "$o.tsv" || { echo "c$n: $o.tsv differs from the smallest-successor blocks' MEMs"; exit 1; }
  done
  if [ -n "${CHECK_PLAIN:-}" ] && [ ! -s "$W/q$n.plain.log" ]; then
    stamp "c$n: plain search"
    $CAP $K find -L $L --format mems --csa -j 4 "$P" "$Q" "$W/q$n.plain.tsv" > "$W/q$n.plain.log" 2>&1
    cmp -s "$W/q${n}s.blocks.tsv" "$W/q$n.plain.tsv" || { echo "c$n: two-level search differs from plain search"; exit 1; }
  fi
done
python3 - "$W" "$NS" "$QLEN" <<'PY'
import os, re, sys
W, NS, QLEN = sys.argv[1], sys.argv[2].split(), int(sys.argv[3])
def grab(f, rx, cast=float):
    m = re.search(rx, open(f).read()) if os.path.exists(f) else None
    return tuple(cast(x) for x in m.groups()) if m else None
rows = []
for n in NS:
    for rule, idx, lvl in (("smallest", "index%s.log" % n, "c%s.lvl" % n), ("frequent", "index%sf.log" % n, "c%sf.lvl" % n)):
        I = os.path.join(W, idx); r = "s" if rule == "smallest" else "f"
        par = grab(I, r"(\d+) phrases, (\d+) distinct", int) or (0, 0)
        blk = grab(I, r"blocks: (\d+) \(([\d.]+) phrases each\), (\d+) symbols; block BWT (\d+) runs") or (0, 0, 0, 0)
        blk = (int(blk[0]), blk[1], int(blk[2]), int(blk[3]))
        mb = grab(I, r"two-level: parse ([\d.]+) MB, blocks ([\d.]+) MB") or (0, 0)
        st = {}
        for mode in ("blocks", "nojumps"):
            f = os.path.join(W, "q%s%s.%s.log" % (n, r, mode))
            st[mode] = (grab(f, r"char_steps=(\d+) parse_steps=(\d+) block_steps=(\d+)", int) or (0, 0, 0),
                        (grab(f, r"wall ([\d.]+) s,") or (0,))[0], (grab(f, r"MEMs=(\d+)", int) or (0,))[0])
        (c, p, b), t, mems = st["blocks"]; (_, pn, _), tn, _ = st["nojumps"]
        m, nb = par[0], blk[0]
        rows.append(dict(n=n, rule=rule, phrases=m, distinct=par[1], blocks=nb, per_block=blk[1],
                         cont=(m - nb) / m if m else 0, symbols=blk[2], bbwt=blk[3],
                         lvl=os.path.getsize(os.path.realpath(os.path.join(W, lvl))) / 1e6 if os.path.exists(os.path.join(W, lvl)) else 0,
                         parse_mb=mb[0], blocks_mb=mb[1], char=c / QLEN, parse=p / QLEN, block=b / QLEN,
                         seen=(pn - p) / b if b else 0, time=t, time_nj=tn, mems=mems))
keys = ["n", "rule", "phrases", "distinct", "blocks", "per_block", "cont", "symbols", "bbwt", "lvl", "parse_mb", "blocks_mb",
        "char", "parse", "block", "seen", "time", "time_nj", "mems"]
with open(os.path.join(W, "blocks.tsv"), "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows: o.write("\t".join(str(r[k]) for k in keys) + "\n")
print("pattern characters: %d (queries.fa)\n" % QLEN)
h = ["n", "rule", "phrases", "distinct", "blocks", "phr/block", "contin.", "lvl MB", "char/ch", "parse/ch", "block/ch", "phr/blockstep", "time s", "no-jumps s"]
fmt = lambda r: [r["n"], r["rule"], "{:,}".format(r["phrases"]), "{:,}".format(r["distinct"]), "{:,}".format(r["blocks"]),
                 "%.2f" % r["per_block"], "%.3f" % r["cont"], "%.1f" % r["lvl"], "%.4f" % r["char"], "%.4f" % r["parse"],
                 "%.4f" % r["block"], "%.2f" % r["seen"], "%.2f" % r["time"], "%.2f" % r["time_nj"]]
tab = [h] + [fmt(r) for r in rows]
w = [max(len(str(x[i])) for x in tab) for i in range(len(h))]
for k, x in enumerate(tab):
    print("  ".join(str(v).rjust(w[i]) for i, v in enumerate(x)))
    if k == 0: print("  ".join("-" * v for v in w))
PY
