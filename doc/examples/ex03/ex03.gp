set terminal pngcairo size 700,650 enhanced font "Arial,12"
set output "ex03.png"

set title "BLEND example 03: cosine taper in a triangle"
set xlabel "x"
set ylabel "y"
set size ratio -1
set xrange [0:99]
set yrange [0:99]
set cbrange [0:1]
set cblabel "Weight"
set palette rgb 33,13,10

plot "ex03.txt" using 1:2:3 with image title ""
