#!/bin/sh
set -eu

if [ "${BLEND:-}" ]; then
    blend_program=$BLEND
elif [ -x "../../../build/blend" ]; then
    blend_program="../../../build/blend"
else
    blend_program="blend"
fi

rm -f ex02_window2d.txt

"$blend_program" window2d -R0/99/0/99 -I1 -Bisotoxal_star.blend | \
awk '
    NR == 1 {y = $2}
    NR > 1 && $2 != y {print ""; y = $2}
    {print}
' > ex02_window2d.txt

printf 'Wrote ex02_window2d.txt\n'
