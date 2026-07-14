#!/bin/sh
set -eu

GMT=${GMT:-gmt}
output="south_america_gshhg_segments.txt"

"$GMT" coast -M -R-85/-30/-60/15 -Dc -A50000 -W > "$output"

printf 'Wrote %s\n' "$output"
printf 'This file contains GMT/GSHHG coastline segments, not one closed BLEND polygon.\n'
printf 'Use south_america.txt for the self-contained monotone example.\n'
