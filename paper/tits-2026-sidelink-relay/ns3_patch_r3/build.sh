#!/bin/bash
# Standalone build of the Paper 8 binary; touches nothing under ns-3-dev.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
B=/home/uzair/6g_ntn_ns3/ns-3-dev/build
OUT="$HERE/../bin"
python3 "$HERE/apply_edits.py"
FLAGS="-DEIGEN_MPL2_ONLY -DHAVE_BOOST -DHAVE_BOOST_UNITS -DHAVE_EIGEN3 -DHAVE_GSL -DHAVE_LIBXML2 -DHAVE_PACKET_H -DHAVE_SQLITE3 -DNS3_BUILD_PROFILE_OPTIMIZED -D__LINUX__ -DP8_STANDALONE -DNtnSlUeMac=NtnSlUeMacP8 -DNtnSlChannel=NtnSlChannelP8 -I$HERE -I$B/include -iquote $B/include/ns3 -I/usr/include/eigen3 -I/usr/include/libxml2 -O3 -DNDEBUG -std=c++20 -fPIE -fno-semantic-interposition -march=native -mtune=native"
LIBS=$(ls $B/lib/libns3.43-*-optimized.so | grep -E "libns3.43-(applications|ntn-constellation|core|network|mobility|point-to-point|internet|stats|bridge|traffic-control)-optimized" | sort)
nice -n 10 g++ $FLAGS -c "$HERE/ntn-nr-sidelink.cc" -o "$OUT/sl.o"
nice -n 10 g++ $FLAGS -c "$HERE/ntn-v2x-canyon-relay.cc" -o "$OUT/ex.o"
nice -n 10 g++ -o "$OUT/p8-canyon-relay3" "$OUT/ex.o" "$OUT/sl.o" $LIBS -Wl,-rpath,$B/lib -Wl,-rpath-link,$B/lib
rm -f "${OUT:?}/sl.o" "${OUT:?}/ex.o"
sha256sum "$OUT/p8-canyon-relay3"
