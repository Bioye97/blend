#!/bin/sh
set -eu

if [ "${BLEND:-}" ]; then
    blend_program=$BLEND
elif [ -x "../../../build/blend" ]; then
    blend_program="../../../build/blend"
else
    blend_program="blend"
fi

rm -f ex03_window2d.txt ex03_window2d_grid.txt ex03_window2d_nonzero.txt

"$blend_program" window2d -R-85/-30/-60/15 -I0.1 -Bsouth_america.blend -Cp -ME -N | \
awk -v nonzero="ex03_window2d_nonzero.txt" '
    $3 > 0 {print > nonzero}
    NR == 1 {y = $2}
    NR > 1 && $2 != y {print ""; y = $2}
    {print}
' > ex03_window2d_grid.txt

printf 'Wrote ex03_window2d_grid.txt and ex03_window2d_nonzero.txt\n'
