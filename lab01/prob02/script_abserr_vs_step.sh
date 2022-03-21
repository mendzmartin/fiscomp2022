#!/usr/bin/gnuplot

# Potential only #

set terminal png size 1024,1024
set output 'error.png'
set xlabel 'step'
set ylabel 'absolute error'
unset key
set log x
set log y
set xrange [0:1]
set yrange [0:1]
set title "1 Potential of Harmonic Oscilator with Sine DVR" offset 0,-0.5
plot 'result.dat' using 2:5 with lines
