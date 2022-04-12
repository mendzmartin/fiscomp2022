#!/usr/bin/gnuplot

######################################################
# GRAFICAR TRANSFORMADA DE FOURIER
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'flips.png'

	# generate de axis titles
	set xlabel 'theta1' offset 0,-0.5 font 'Times Roman Bold Italic,20'
	set ylabel 'theta2' offset -0.5,0 font 'Times Roman Bold Italic,20'

	set pm3d map 

	# define color palette
	set palette defined (0 "green", 10 "#000F00", 10 "#FF0000", \
	100 "#310000", 100 "purple",1000 "#54025C",1000 "#0000FF", \
	10000 "#000B70", 10000 "white",10001 "white")

	set cbrange [0:10000]
	set size ratio .5
	set xrange [-3:3]
	set yrange [0:3]
	#unset surf

	set title 'flip color map' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	splot 'result_flips.dat' matrix
	# reset all of graph characteristics to default values
	reset
