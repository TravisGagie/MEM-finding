#!/bin/bash
# Compiles stpd-index-mems against a built clone of https://github.com/regindex/STPD-index, with the flags
# and libraries of its own stpd-locate (CMake fetches sdsl and malloc_count into STPDIR/build/_deps).
#   build_stpd_index_mems.sh STPDIR OUT_BINARY
set -euo pipefail
D=$(cd "$1" && pwd); OUT=$2
SRC=$(cd "$(dirname "$0")" && pwd)/stpd_index_mems.cpp
B=$D/build
[ -f "$B/CLI/stpd-locate" ] || { echo "build STPD-index first: cd $D && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=\"-I\$PWD/_deps/sdsl-build/external/libdivsufsort/include\" && make -j4"; exit 1; }
INC="-I$D/include -I$B/_deps/sdsl-src/include -I$B/_deps/malloc_count-src -I$D/external/sux"
for d in common bitvector elias_fano phi_function move text_oracle text_oracle/RLZ binary_search text_index; do INC="$INC -I$D/include/$d"; done
c++ -I$B/_deps/sdsl-build/external/libdivsufsort/include -O3 -DNDEBUG -march=native -funroll-loops -fomit-frame-pointer -w -std=c++20 \
  $INC "$SRC" -o "$OUT" "$B/external/libmalloc_count.a" -ldl "$B/_deps/sdsl-build/lib/libsdsl.a" \
  "$B/_deps/sdsl-build/external/libdivsufsort/lib/libdivsufsort.a" "$B/_deps/sdsl-build/external/libdivsufsort/lib/libdivsufsort64.a"
echo "built $OUT"
