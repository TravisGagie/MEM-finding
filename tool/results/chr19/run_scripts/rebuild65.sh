#!/bin/bash
# step 1: c65 with the new two-level build (tool/NOTES.md, "For the desktop (2026-11-03): the new two-level build")
cd ~/research/rlbwt-tricks/chr19; K=~/research/mem-finding/tool/bin/katka; BB=~/research/mem-finding/katka/deps/Big-BWT/bigbwt
echo "== bigbwt c65 (c65.bwt was deleted for disk space) $(date +%T)"
if [ ! -s c65.bwt ]; then systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "bigbwt: %e s wall, peak %M KB" $BB c65.S > bigbwt.65.err 2>&1 || { tail bigbwt.65.err; exit 1; }; mv c65.S.bwt c65.bwt; rm -f c65.S.log c65.S.bigbwt.log; tail -1 bigbwt.65.err; fi
echo "== index --preset chr19 c65 $(date +%T)"
mv c65.lvl c65.lvl.prev 2>/dev/null
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "index: %e s wall, peak %M KB" $K index --preset chr19 c65 > index.65.new.err 2>&1; echo "exit $?"
cat index.65.new.err; ls -l c65.*
cmp c65.lvl c65.lvl.old && echo "same .lvl as the old build" || echo "DIFFERENT .lvl"
rm -f c65.lvl.prev; df -h / | tail -1; echo "== done $(date +%T)"
