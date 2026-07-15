#!/bin/sh
set -eu

if [ "${BLEND:-}" ]; then
    blend_program=$BLEND
elif [ -x "../../../build/blend" ]; then
    blend_program="../../../build/blend"
else
    blend_program="blend"
fi

"$blend_program" window1d -R0/10 -I0.1 -Fconnes -T0.3/0.3 > connes.txt
printf 'Wrote connes.txt\n'
