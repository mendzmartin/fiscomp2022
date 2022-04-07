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

	set autoscale x
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

	set title 'Fourier transform of f(t) = sin(t*pi/2)+cos(t*20pi). N=2^{10}, T=4' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'tfx.d' using 1:2 title 'fftw3: real part' with points linestyle 1,\
			'tfx.d' using 1:3 title 'fftw3: imaginary part' with points linestyle 2

	# reset all of graph characteristics to default values
	reset
	
	######################################################
	# GRAFICAR SOLUCION (POSICIÓN)
	######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'function.png'

	# generate de axis titles
	set xlabel 't' font 'Times Roman Bold Italic,20'
	set ylabel 'f(t)' font 'Times Roman Bold Italic,20'

	set autoscale x
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

	set title 'function - time domain' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'fx.d' using 1:2 title 'f(t) = sin(t*pi/2)+cos(t*20pi). N=2^{10}, T=4' with linespoints linestyle 1,\
			'fx_back.d' using 1:2 title 'inverse fourier transform' with linespoints linestyle 2,\
			(sin(3.14*0.5*x)+cos(20*3.14*x)) title 'explicit expresion'

	# reset all of graph characteristics to default values
	reset
