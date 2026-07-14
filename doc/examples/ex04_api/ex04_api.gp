set terminal pngcairo size 1800,1150 enhanced font "Helvetica,18"
set output "ex04_api.png"
set title "BLEND C API: South America 2-D support"
set xlabel "Longitude"
set ylabel "Latitude"
set xrange [-85:-30]
set yrange [-60:15]
set cbrange [0:1]
set palette rgb 33,13,10
set colorbox vertical user origin 0.90,0.20 size 0.020,0.58
set cblabel "weight" offset 1.5,0
set border lw 1.4
set tics out nomirror
set grid lw 0.8 lc rgb "#d9dee7"
set size ratio -1
set key outside right top spacing 1.25
set view map
set pm3d map interpolate 1,1 corners2color mean

splot "ex04_api_grid.txt" using 1:2:3 with pm3d notitle, \
      "south_america.txt" using 1:2:(1.01) with lines lw 1.8 lc rgb "#6b6b6b" title "South America input", \
      "south_america_monotone.txt" using 1:2:(1.02) with lines lw 3.0 lc rgb "#762a83" title "South America support"
