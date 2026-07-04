set terminal pngcairo size 700,650 enhanced font "Arial,12"
set output "ex06.png"

set title "BLEND example 06: cosine taper in a hexagon"
set xlabel "x"
set ylabel "y"
set size ratio -1
set xrange [0:99]
set yrange [0:99]
set cbrange [0:1]
set cblabel "Weight"
set palette rgb 33,13,10

plot "ex06.txt" using 1:2:3 with image title ""
