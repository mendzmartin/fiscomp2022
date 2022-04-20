#!/usr/bin/gnuplot
	set terminal png size 1024,1024
	set output 'prob_01.png'
	set ylabel 'longitud de la barra [m]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
	set xlabel 'tiempo [s]' offset -0.5,0 font 'Times Roman Bold Italic,20'
	set pm3d map
	# set palette defined (0 "green", 25 "#000F00", 25 "#FF0000", \
	# 50 "#310000", 50 "purple",75 "#54025C",75 "#0000FF", \
	# 100 "#000B70")
	set cbrange [0:100]
	set xrange [0:153750000]
	set yrange [0:1]
	set title 'Distribución de temperaturas' offset 0,-0.5 font 'Times Roman Bold Italic,20'
	plot for [i=1:50] '../results/result_01.dat' using (300*(i-1)*0.1025E+05):1:2 every::(100*(i-1)+1)::(100*i) with points pt 5 ps 2 palette t ''
	reset
