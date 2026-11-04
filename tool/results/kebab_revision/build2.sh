#!/bin/bash
# kebab build with one thread: with more, kebab's process_sequences strdup's every sequence and frees only the last
# (4bed76f), so memory grows with the input and the 11 GB cgroup killed it at 43%. One thread uses kseq's buffer.
cd ~/zymo/rb; X="/run/media/travis/Travis backup/zymo_plain"
echo "== kebab-build $(date +%T)"
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "kebab-build: %e s wall, %U s user, peak %M KB" bin/kebab build -k 20 -e 0.1 -f 1 -t 1 -o zymo "$X/ref.fa" > kebab-build.err 2>&1; rc=$?
tr '\r' '\n' < kebab-build.err | grep -vE 'Estimating Cardinality: +[0-9.]+%$|Indexing: +[0-9.]+%$' | tail -12; echo "exit $rc"; [ $rc = 0 ] || exit 1
ls -l zymo.*
cd ~/rz/MEM-finding && git pull -q origin tool; cd ~/zymo/rb && bash run.sh > run.log 2>&1; echo "run exit $?"
