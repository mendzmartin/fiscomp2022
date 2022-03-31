#!/usr/bin/gnuplot

######################################################
# GRAFICAR SOLUCION (POSICIÓN)
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'solution_position.png'

	# generate de axis titles
	set xlabel 'time [s]' font 'Times Roman Bold Italic,20'
	set ylabel 'position [m]' font 'Times Roman Bold Italic,20'

	set xrange [0:10]
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


	set title 'Soluciones EDO 1er orden para n = 131072' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'result_y1.dat' using 2:3 every ::131071::262142 title 'exact solution' with lines linestyle 1,\
			'result_y1.dat' using 2:4 every ::131071::262142 title 'euler' with lines linestyle 2,\
			'result_y1.dat' using 2:5 every ::131071::262142 title 'RK2_{hu}' with lines linestyle 3,\
			'result_y1.dat' using 2:6 every ::131071::262142 title 'RK2_{mp}' with lines linestyle 4,\
			'result_y1.dat' using 2:7 every ::131071::262142 title 'RK2_{ra}' with lines linestyle 5,\
			'result_y1.dat' using 2:8 every ::131071::262142 title 'RK4_{cl}' with lines linestyle 6

	# reset all of graph characteristics to default values
	reset

######################################################
# GRAFICAR SOLUCION (VELOCIDAD)
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'solution_velocity.png'

	# generate de axis titles
	set xlabel 'time [s]' font 'Times Roman Bold Italic,20'
	set ylabel 'velocity [m/s]' font 'Times Roman Bold Italic,20'

	set xrange [0:10]
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


	set title 'Soluciones EDO 1er orden para n = 131072' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'result_y2.dat' using 2:3 every ::131071::262142 title 'exact solution' with lines linestyle 1,\
			'result_y2.dat' using 2:4 every ::131071::262142 title 'euler' with lines linestyle 2,\
			'result_y2.dat' using 2:5 every ::131071::262142 title 'RK2_{hu}' with lines linestyle 3,\
			'result_y2.dat' using 2:6 every ::131071::262142 title 'RK2_{mp}' with lines linestyle 4,\
			'result_y2.dat' using 2:7 every ::131071::262142 title 'RK2_{ra}' with lines linestyle 5,\
			'result_y2.dat' using 2:8 every ::131071::262142 title 'RK4_{cl}' with lines linestyle 6

	# reset all of graph characteristics to default values
	reset

######################################################
# GRAFICAR ERRORES
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'solution_rel_err.png'

	# generate de axis titles
	set xlabel 'time [s]' font 'Times Roman Bold Italic,20'
	set ylabel 'relative error [%] - logarithmic scale' font 'Times Roman Bold Italic,20'

	set key left top
	
	set logscale y
	
	set xrange [0:10]
	set autoscale y

	# plots functions with linespoints
	#set style function linespoints

	# define linestyles
	set style line 1 linetype rgb "red" linewidth 3
	set style line 2 linetype rgb "green" linewidth 3
	set style line 3 linetype rgb "blue" linewidth 3
	set style line 4 linetype rgb "violet" linewidth 3
	set style line 5 linetype rgb "sienna4" linewidth 3


#	set title 'Soluciones EDO 1er orden para n = 1024 (h->O(10^{-3}))' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
#	plot	'result_y1.dat' using 2:9 every ::1023::2046 title 'euler' with lines linestyle 1,\
#			'result_y1.dat' using 2:10 every ::1023::2046 title 'RK2_{hu}' with lines linestyle 2,\
#			'result_y1.dat' using 2:11 every ::1023::2046 title 'RK2_{mp}' with lines linestyle 3,\
#			'result_y1.dat' using 2:12 every ::1023::2046 title 'RK2_{ra}' with lines linestyle 4,\
#			'result_y1.dat' using 2:13 every ::1023::2046 title 'RK4_{cl}' with lines linestyle 5
	
#	set pointsize 0.3	
#	set title 'Soluciones EDO 1er orden para n = 512 (h->O(10^{-2}))' offset 0,-0.5 font 'Times Roman Bold Italic,20'
#	# create plot
#	plot	'result_y1.dat' using 2:9 every ::511::1022 title 'euler' with points linestyle 1,\
#			'result_y1.dat' using 2:10 every ::511::1022 title 'RK2_{hu}' with points linestyle 2,\
#			'result_y1.dat' using 2:11 every ::511::1022 title 'RK2_{mp}' with points linestyle 3,\
#			'result_y1.dat' using 2:12 every ::511::1022 title 'RK2_{ra}' with points linestyle 4,\
#			'result_y1.dat' using 2:13 every ::511::1022 title 'RK4_{cl}' with points  linestyle 5

	set pointsize 0.3
	set title 'Soluciones EDO 1er orden para n = 131072 (h->O(10^{-5}))' offset 0,-0.5 font 'Times Roman Bold Italic,20'
#	# create plot
	plot	'result_y1.dat' using 2:9 every ::131071::262142 title 'euler' with points linestyle 1,\
			'result_y1.dat' using 2:10 every ::131071::262142 title 'RK2_{hu}' with points linestyle 2,\
			'result_y1.dat' using 2:11 every ::131071::262142 title 'RK2_{mp}' with points linestyle 3,\
			'result_y1.dat' using 2:12 every ::131071::262142 title 'RK2_{ra}' with points linestyle 4,\
			'result_y1.dat' using 2:13 every ::131071::262142 title 'RK4_{cl}' with points  linestyle 5

	# reset all of graph characteristics to default values
	reset

######################################################
# GRAFICAR ENERGIAS
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'solution_energy.png'

	# generate de axis titles
	set xlabel 'time [s]' font 'Times Roman Bold Italic,20'
	set ylabel 'total energy [J]' font 'Times Roman Bold Italic,20'

	set xrange [0:10]
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


	set title 'Soluciones EDO 1er orden para n = 131072' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'result_energies.dat' using 2:3 every ::131071::262142 title 'exact solution' with lines linestyle 1,\
			'result_energies.dat' using 2:4 every ::131071::262142 title 'euler' with lines linestyle 2,\
			'result_energies.dat' using 2:5 every ::131071::262142 title 'RK2_{hu}' with lines linestyle 3,\
			'result_energies.dat' using 2:6 every ::131071::262142 title 'RK2_{mp}' with lines linestyle 4,\
			'result_energies.dat' using 2:7 every ::131071::262142 title 'RK2_{ra}' with lines linestyle 5,\
			'result_energies.dat' using 2:8 every ::131071::262142 title 'RK4_{cl}' with lines linestyle 6
			
	

	# reset all of graph characteristics to default values
	reset
	
######################################################
# GRAFICAR ERRORES
######################################################
	# kind of output to generate
	set terminal png size 1024,1024
	show terminal
	# redirect the output toa file or device
	set output 'solution_rel_err_glob.png'

	# generate de axis titles
	set xlabel 'time [s]' font 'Times Roman Bold Italic,20'
	set ylabel 'global relative error [%] - logarithmic scale' font 'Times Roman Bold Italic,20'

	set key left top
	
	set logscale y
	
	set autoscale x
	set autoscale y

	# plots functions with linespoints
	#set style function linespoints

	# define linestyles
	set style line 1 linetype rgb "red" linewidth 3
	set style line 2 linetype rgb "green" linewidth 3
	set style line 3 linetype rgb "blue" linewidth 3
	set style line 4 linetype rgb "violet" linewidth 3
	set style line 5 linetype rgb "sienna4" linewidth 3


	set title 'Soluciones EDO 1er orden' offset 0,-0.5 font 'Times Roman Bold Italic,20'

	# create plot
	plot	'result_global_error.dat' using (($1-1)/10):9 title 'euler' with lines linestyle 1,\
			'result_global_error.dat' using (($1-1)/10):10 title 'RK2_{hu}' with lines linestyle 2,\
			'result_global_error.dat' using (($1-1)/10):11 title 'RK2_{mp}' with lines linestyle 3,\
			'result_global_error.dat' using (($1-1)/10):12 title 'RK2_{ra}' with lines linestyle 4,\
			'result_global_error.dat' using (($1-1)/10):13 title 'RK4_{cl}' with lines linestyle 5

	# reset all of graph characteristics to default values
	reset
