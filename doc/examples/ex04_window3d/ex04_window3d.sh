#!/bin/sh
set -eu

if [ "${BLEND:-}" ]; then
    blend_program=$BLEND
elif [ -x "../../../build/blend" ]; then
    blend_program="../../../build/blend"
else
    blend_program="blend"
fi

rm -f ex04_window3d.txt ex04_window3d_nonzero.txt

"$blend_program" window3d -R-85/-30/-60/15/0/80 -I0.1/0.1/5 -Bsouth_america.blend -Cp -ME -N | \
awk '$4 > 0 {print}' > ex04_window3d_nonzero.txt

printf 'Wrote ex04_window3d_nonzero.txt\n'
