#!/bin/sh
set -eu

cc_program=${CC:-cc}
prefix=${BLEND_PREFIX:-}

if [ -n "$prefix" ]; then
    "$cc_program" ex01_api.c -I "$prefix/include/blend" -L "$prefix/lib" -lblend -o ex01_api
elif [ -f "../../../build/libblend.a" ]; then
    "$cc_program" ex01_api.c ../../../build/libblend.a -I ../../../src -lm -o ex01_api
else
    "$cc_program" ex01_api.c ../../../src/blend_boundary.c ../../../src/blend_contribution.c \
        ../../../src/blend_error.c ../../../src/blend_polygon.c ../../../src/blend_report.c \
        ../../../src/blend_utils.c ../../../src/blend_window.c -I ../../../src -lm -o ex01_api
fi

./ex01_api
if command -v gnuplot >/dev/null 2>&1; then
    gnuplot ex01_api.gp
fi
printf 'Wrote ex01_api.txt and ex01_api.png\n'
