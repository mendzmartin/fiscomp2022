#!/usr/bin/gnuplot

######################################################
# GRAFICAR TRANSFORMADA DE FOURIER
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'fourier_transform.png'

	# generate de axis titles
	set xlabel 'w' font 'Times Roman Bold Italic,20'
	set ylabel 'F(w)' font 'Times Roman Bold Italic,20'

	set xrange[0:50]
	set autoscale y

	# plots functions with linespoints
	#set style function linespoints

	# define linestyles
	set style line 1 linetype rgb "red" linewidth 3
	set style line 2 linetype rgb "green" linewidth 3
	set style line 3 linetype rgb "blue" linewidth 3
	set style line 4 linetype rgb "violet" linewidth 3
	set style line 5 linetype rgb "black" linewidth 3
	set style line 6 linetype rgb "sienna4" linewidth 3
	set style line 7 linetype rgb "brown" linewidth 3
	set style line 8 linetype rgb "grey" linewidth 3
	set style line 9 linetype rgb "cyan" linewidth 3

	set title 'Fourier transform of f(t) = sin(t*pi/2)+cos(t*20pi). N=2^{10}, T=4' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'result_logistic_map.dat' using 1:2 every::1::512 title 'x1' smooth unique with linespoints linestyle 1,\
			'result_logistic_map.dat' using 1:2 every::513::1024 title 'x2' smooth unique with linespoints linestyle 2,\
			'result_logistic_map.dat' using 1:2 every::1025::1536 title 'x3' smooth unique with linespoints linestyle 3,\
			'result_logistic_map.dat' using 1:2 every::1537::2048 title 'x4' smooth unique with linespoints linestyle 4,\

	# reset all of graph characteristics to default values
	reset
