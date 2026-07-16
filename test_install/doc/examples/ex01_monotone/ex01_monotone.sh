#!/bin/sh
set -eu

if [ "${BLEND:-}" ]; then
    blend_program=$BLEND
elif [ -x "../../../build/blend" ]; then
    blend_program="../../../build/blend"
else
    blend_program="blend"
fi

input="south_america.txt"
status_output="ex01_status.txt"
envelope_output="south_america_envelope.txt"
best_output="south_america_monotone.txt"
best_plot_output="south_america_monotone_plot.txt"

"$blend_program" monotone "$input" > "$status_output"
"$blend_program" monotone -Me "$input" > "$envelope_output"
"$blend_program" monotone -Mb -G256/256 "$input" > "$best_output"
awk 'NR == 1 {first = $0} {print} END {if (NR > 0) print first}' "$best_output" > "$best_plot_output"

printf 'Wrote %s, %s, %s, and %s\n' "$status_output" "$envelope_output" "$best_output" "$best_plot_output"
