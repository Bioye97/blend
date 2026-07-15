set terminal pngcairo size 1200,720 enhanced font "Helvetica,16"
set output "bartlett.png"
set title "bartlett window, taper ratio 0.3/0.3"
set xlabel "x"
set ylabel "weight"
set xrange [0:10]
set yrange [0:1.05]
set border lw 1.4
set tics out nomirror
set grid lw 1 lc rgb "#d8dee9"
set key off
plot "bartlett.txt" using 1:2 with lines lw 4 lc rgb "#176d8f"
