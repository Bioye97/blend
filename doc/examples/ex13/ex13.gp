set terminal pngcairo size 900,700 enhanced font "Arial,12"
set output "ex13.png"
set title "South America polygon and xy-monotone conversions"
set xlabel "Longitude"
set ylabel "Latitude"
set grid
set key outside
set size ratio -1
set xrange [-86:-30]
set yrange [-60:15]

plot "south_america.txt" using 1:2 with filledcurves closed fs transparent solid 0.18 lc rgb "#4c78a8" title "input polygon", \
     "south_america.txt" using 1:2 with linespoints lw 2 pt 7 ps 0.6 lc rgb "#4c78a8" notitle, \
     "ex13_monotone_envelope.txt" using 1:2 with filledcurves closed fs transparent solid 0.16 lc rgb "#f58518" title "-Me envelope", \
     "ex13_monotone_envelope.txt" using 1:2 with linespoints lw 2 pt 5 ps 0.7 lc rgb "#f58518" notitle, \
     "ex13_monotone_best_plot.txt" using 1:2 with filledcurves closed fs transparent solid 0.14 lc rgb "#b279a2" title "-Mb best IoU envelope", \
     "ex13_monotone_best_plot.txt" using 1:2 with linespoints lw 3 pt 11 ps 0.55 lc rgb "#b279a2" notitle
