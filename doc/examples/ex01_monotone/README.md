# ex01_monotone: South America Polygon With `blend monotone`

This example checks a simplified South America polygon with the `monotone`
module, then writes xy-monotone polygons using the `-Me` envelope method and
the `-Mb` best-IoU piecewise envelope method.

The committed `south_america.txt` file is a compact lon/lat polygon fixture for
testing. It is simplified, but dense enough to preserve the main continental
outline. The helper script `fetch_gshhg_reference.sh` can dump a GMT/GSHHG
reference coastline for the same region when GMT and GSHHG data are available,
but that output is a multi-segment coastline, not a single closed BLEND polygon.

Run from this directory:

```sh
./ex01_monotone.sh
```

The script writes:

```text
ex01_status.txt
south_america_envelope.txt
south_america_monotone.txt
south_america_monotone_plot.txt
```

To plot the input polygon and all xy-monotone conversions:

```sh
gnuplot ex01_monotone.gp
```

To force a specific BLEND executable:

```sh
BLEND=/path/to/blend ./ex01_monotone.sh
```
