set terminal pngcairo size 1450,900 enhanced font "Helvetica,16"
set output "ex05_window2d.png"
set title "BLEND window2d: multiple polygon supports"
set xlabel "x"
set ylabel "y"
set xrange [0:100]
set yrange [0:100]
set cbrange [0:1]
set palette rgb 33,13,10
set colorbox vertical user origin 0.90,0.22 size 0.020,0.56
set cblabel "weight" offset 1.5,0
set border lw 1.4
set tics out nomirror
set grid lw 0.8 lc rgb "#d9dee7"
set size ratio -1
set key outside right top spacing 1.25
set view map
set pm3d map interpolate 1,1 corners2color mean

splot "ex05_window2d_grid.txt" using 1:2:3 with pm3d notitle, \
      "ex05_window2d_polygons.txt" using 1:2:(1.02) with lines lw 2.6 lc rgb "#762a83" title "polygon supports"
