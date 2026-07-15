set terminal pngcairo size 700,650 enhanced font "Arial,12"
set output "ex02_window2d.png"

set title "BLEND window2d: 2-D cosine taper in a 4-pointed isotoxal star"
set xlabel "x"
set ylabel "y"
set size ratio -1
set xrange [0:99]
set yrange [0:99]
set cbrange [0:1]
set cblabel "Weight"
set palette rgb 33,13,10

plot "ex02_window2d.txt" using 1:2:3 with image title ""
