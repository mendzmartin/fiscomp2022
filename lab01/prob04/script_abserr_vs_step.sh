#!/usr/bin/gnuplot


set terminal png size 1024,1024
set output 'error.png'
set xlabel 'interval points (n)'
set ylabel 'absolute error [%]'
set key
set log x
set log y
set autoscale x
set autoscale y
set title 'Integración numérica: Análisis de error' offset 0,-0.5
plot	'result.dat' using 1:5  title 'trapezoidal method' with lines, \
		'result.dat' using 1:7  title 'simpson 1/3 method' with lines, \
		'result.dat' using 1:9  title 'simpson 3/8 method' with lines, \
		'result.dat' using 1:11 title 'euler-legendre quadrature method' with lines
