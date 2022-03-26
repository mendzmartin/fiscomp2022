#!/usr/bin/gnuplot


set terminal png size 1024,1024
set output 'error.png'
set xlabel 'step'
set ylabel 'absolute error'
unset key
set log x
set log y
#set xrange [0.000013:1]
#set yrange [0:0.18]
set autoscale x
set autoscale y
set title "Diferenciación numérica: Análisis de error" offset 0,-0.5
plot 'result.dat' using 1:2 with lines
