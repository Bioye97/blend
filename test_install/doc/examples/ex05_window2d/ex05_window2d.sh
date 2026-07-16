#!/bin/sh
set -eu

if [ "${BLEND:-}" ]; then
    blend_program=$BLEND
elif [ -x "../../../build/blend" ]; then
    blend_program="../../../build/blend"
else
    blend_program="blend"
fi

rm -f ex05_window2d_grid.txt ex05_window2d_polygons.txt

"$blend_program" window2d -R0/100/0/100 -I0.1 -Bpolygons.blend -Cp | \
awk '
    NR == 1 {y = $2}
    NR > 1 && $2 != y {print ""; y = $2}
    {print}
' > ex05_window2d_grid.txt

while read -r polygon function taper; do
    [ -n "$polygon" ] || continue
    first_line=$(sed -n '1p' "$polygon")
    sed -n 'p' "$polygon" >> ex05_window2d_polygons.txt
    printf '%s\n\n' "$first_line" >> ex05_window2d_polygons.txt
done < polygons.blend

printf 'Wrote ex05_window2d_grid.txt and ex05_window2d_polygons.txt\n'
