set terminal pngcairo size 1500,1200 enhanced font "Helvetica,18"
set output "ex01_monotone.png"
set title "South America support polygon and xy-monotone conversions"
set xlabel "Longitude"
set ylabel "Latitude"
set border lw 1.4
set tics out nomirror
set grid lw 0.8 lc rgb "#d9dee7"
set key outside right top spacing 1.25
set size ratio -1
set xrange [-86:-30]
set yrange [-60:15]

plot "south_america.txt" using 1:2 with filledcurves closed fs transparent solid 0.18 lc rgb "#2f6f8f" title "input polygon", \
     "south_america.txt" using 1:2 with lines lw 2.8 lc rgb "#2f6f8f" notitle, \
     "south_america_envelope.txt" using 1:2 with filledcurves closed fs transparent solid 0.12 lc rgb "#f28e2b" title "-Me envelope", \
     "south_america_envelope.txt" using 1:2 with lines lw 2.8 lc rgb "#f28e2b" notitle, \
     "south_america_monotone_plot.txt" using 1:2 with filledcurves closed fs transparent solid 0.14 lc rgb "#9467bd" title "-Mb best IoU envelope", \
     "south_america_monotone_plot.txt" using 1:2 with lines lw 3.2 lc rgb "#9467bd" notitle
