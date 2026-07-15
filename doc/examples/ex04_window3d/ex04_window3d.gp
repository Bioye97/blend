set terminal pngcairo size 1700,1250 enhanced font "Helvetica,18"
set output "ex04_window3d.png"
set title "BLEND window3d: South America support"
set xlabel "Longitude"
set ylabel "Latitude"
set zlabel "z"
set xrange [-85:-30]
set yrange [-60:15]
set zrange [0:80]
set cbrange [0:1]
set palette rgb 33,13,10
set colorbox vertical user origin 0.88,0.20 size 0.020,0.56
set cblabel "weight" offset 1.5,0
set border lw 1.4
set tics out nomirror
set grid lw 0.8 lc rgb "#d9dee7"
set xyplane at 0
set view 63,34
set key outside right top spacing 1.25

splot "ex04_window3d_nonzero.txt" every 25 using 1:2:3:4 with points pt 7 ps 0.20 lc palette title "weights", \
      "south_america_monotone.txt" using 1:2:(0) with lines lw 2.8 lc rgb "#762a83" title "South America support"
