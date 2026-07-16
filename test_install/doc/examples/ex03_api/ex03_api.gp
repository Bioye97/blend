set terminal pngcairo size 900,750 enhanced font "Arial,12"
set output "ex03_api.png"

set title "BLEND example 03: 3-D cosine taper in a 4-pointed isotoxal star"
set xlabel "x"
set ylabel "y"
set zlabel "z"
set xrange [0:99]
set yrange [0:99]
set zrange [0:24]
set cbrange [0:1]
set cblabel "Weight"
set palette rgb 33,13,10
set view 62,32
set xyplane 0
set grid
unset key

splot "< awk '$4 >= 0.05 {print $1, $2, $3, $4}' ex03_api.txt" \
  using 1:2:3:4 with points pointtype 7 pointsize 0.25 palette
