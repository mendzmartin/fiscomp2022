#!/usr/bin/gnuplot

	# kind of output to generate
	set terminal png size 1024,1024
	set terminal post color
	#show terminal
	# redirect the output toa file or device
	#set output 'flips.ps'
	set output 'flips.png'

	# generate de axis titles
	set xlabel 'theta1' offset 0,-0.5 font 'Times Roman Bold Italic,20'
	set ylabel 'theta2' offset -0.5,0 font 'Times Roman Bold Italic,20'

	#set view map
	#set dgrid3d
	set pm3d map

	#define color palette
	set palette defined (0 "green", 10 "#000F00", 10 "#FF0000", \
	100 "#310000", 100 "purple",1000 "#54025C",1000 "#0000FF", \
	10000 "#000B70", 10000 "white",10001 "white")

	#set palette defined (0 "green", 0.5 "#000F00", 0.5 "#FF0000", \
	#1 "#310000", 1 "purple",1.5 "#54025C",1.5 "#0000FF", \
	#2 "#000B70", 2 "white",2.5 "white")

	set cbrange [0:10000]
	#set cbrange [0:2.5]
#	set size ratio .5
	set xrange [-3:3]
	set yrange [0:3]
	#unset surf

	set title 'flip color map' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
#	splot 'result_flips.dat' using 1:2:3 with pm3d
	plot 'result_flips.dat' using 1:2:3 with points pt 5 ps 0.2 palette t ''

	# reset all of graph characteristics to default values
	reset
