#!/bin/sh
set -eu

cc_program=${CC:-cc}
prefix=${BLEND_PREFIX:-}

rm -f ex04_api ex04_api_grid.txt ex04_api.png south_america_monotone.txt

if [ -n "$prefix" ]; then
    "$cc_program" ex04_api.c -I "$prefix/include/blend" -L "$prefix/lib" -lblend -o ex04_api
elif [ -f "../../../build/libblend.a" ]; then
    "$cc_program" ex04_api.c ../../../build/libblend.a -I ../../../src -lm -o ex04_api
else
    "$cc_program" ex04_api.c ../../../src/blend_boundary.c ../../../src/blend_contribution.c \
        ../../../src/blend_error.c ../../../src/blend_polygon.c ../../../src/blend_report.c \
        ../../../src/blend_utils.c ../../../src/blend_window.c -I ../../../src -lm -o ex04_api
fi

./ex04_api
if command -v gnuplot >/dev/null 2>&1; then
    gnuplot ex04_api.gp
fi
printf 'Wrote ex04_api_grid.txt, south_america_monotone.txt, and ex04_api.png\n'
