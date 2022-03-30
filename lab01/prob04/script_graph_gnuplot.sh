#!/usr/bin/gnuplot

# kind of output to generate
set terminal png size 1024,1024
show terminal
# redirect the output toa file or device
set output 'error.png'

# generate de axis titles
set xlabel 'interval points (n) - logarithmic scale' font 'Times Roman Bold Italic,20'
set ylabel 'absolute error [%] - logarithmic scale' font 'Times Roman Bold Italic,20'

# set logarithmic scale
set log x
set log y

#set xrange [2:32768]
set autoscale x
set yrange [10**-13:10]
#set autoscale y

# plots functions with linespoints
#set style function linespoints

# define linestyles
set style line 1 linetype rgb "red" linewidth 3
set style line 2 linetype rgb "green" linewidth 3
set style line 3 linetype rgb "blue" linewidth 3
set style line 4 linetype rgb "violet" linewidth 3
#set style line 5 pointtype 5 pointsize 0.1


set title 'Numerical integration: error analysis' offset 0,-0.5 font 'Times Roman Bold Italic,20'

# define functions
f(x) = (1/(x**2))*5
g(x) = (1/(x**4))*1
h(x) = (1/(x**1))*8
i(x) = sqrt(x)*(10**-13)

# create plot
plot	'result_odd_intervals.dat' using 1:5  title 'trapezoidal method' with lines linestyle 1,\
		'result_odd_intervals.dat' using 1:7  title 'simpson 1/3 method' with lines linestyle 2,\
		'result_odd_intervals.dat' using 1:9  title 'simpson 3/8 method' with lines linestyle 3,\
		'result_odd_intervals.dat' using 1:11 title 'euler-legendre quadrature method' with lines linestyle 4,\
		f(x) title 'x^-2' with points pointtype 1 pointsize 0.8 linecolor rgb "black",\
		g(x) title 'x^-4' with points pointtype 2 pointsize 0.8 linecolor rgb "black",\
		h(x) title 'x^-1' with points pointtype 3 pointsize 0.8 linecolor rgb "black",\
		i(x) title 'x^.5' with points pointtype 4 pointsize 0.8 linecolor rgb "black"

# reset all of graph characteristics to default values
reset
