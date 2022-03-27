#!/usr/bin/gnuplot


set terminal png size 1024,1024

set output 'error.png'

set xlabel 'interval points (n) - logarithmic scale'
set ylabel 'absolute error [%] - logarithmic scale'

set key

set log x
set log y

#set xrange [2:32768]
set autoscale x
set autoscale y

set style line 1 linetype rgb "red"		linewidth 3
set style line 2 linetype rgb "green"	linewidth 3
set style line 3 linetype rgb "blue"	linewidth 3
set style line 4 linetype rgb "violet"	linewidth 3
set style line 5 linetype rgb "black"	linewidth 1 dashtype '_  ____  _'

#set style line 2 linetype rgb "orange" linewidth 1 dashtype '-- '
#set style line 3 linetype rgb "yellow" linewidth 3
#set style line 5 linetype rgb "cyan" linewidth 4
#set style line 8 linecolor rgb '#0060ad' linetype 1 linewidth 2


set title 'Integración numérica: Análisis de error' offset 0,-0.5

f(x) = (1/(x**2))*5
g(x) = (1/(x**4))*1
h(x) = (1/(x**1))*10
plot	'result.dat' using 1:5  title 'trapezoidal method' 					with lines linestyle 1, \
		'result.dat' using 1:7  title 'simpson 1/3 method' 					with lines linestyle 2, \
		'result.dat' using 1:9  title 'simpson 3/8 method' 					with lines linestyle 3, \
		'result.dat' using 1:11 title 'euler-legendre quadrature method' 	with lines linestyle 4, \
		f(x) title 'x^-2' with lines linestyle 5, \
		g(x) title 'x^-4' with lines linestyle 5, \
		h(x) title 'x^-1' with lines linestyle 5
