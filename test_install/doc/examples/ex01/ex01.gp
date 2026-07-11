set terminal pngcairo size 800,500 enhanced font "Arial,12"
set output "ex01.png"

set title "BLEND example 01: symmetric cosine taper"
set xlabel "Sample"
set ylabel "Weight"
set grid
set xrange [1:100]
set yrange [0:1.05]

plot "ex01.txt" using 1:2 with lines linewidth 2 title "cosine"
