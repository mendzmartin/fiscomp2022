#!/usr/bin/gnuplot
set terminal png size 1024,1024
set ylabel 'longitud de la barra [m]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set xlabel 'tiempo [s]' offset -0.5,0 font 'Times Roman Bold Italic,20'
set pm3d map
set cbrange [0:100]
set xrange [0:153750000]
set yrange [0+1.00E-03:1+1.00E-03]

set title 'Distribución de temperaturas: Método Explícito' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'explicit.png'
splot for [i=1:50] '../results/result_01_explicit.dat' using (300*(i-1)*0.1025E+05):1:2 every::(100*(i-1)+1)::(100*i) with points pt 5 ps 2 palette t ''

set title 'Distribución de temperaturas: Método Implicito' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'implicit.png'
splot for [i=1:50] '../results/result_01_implicit.dat' using (300*(i-1)*0.1025E+05):1:2 every::(100*(i-1)+1)::(100*i) with points pt 5 ps 2 palette t ''

set title 'Distribución de temperaturas: Método Crank-Nicolson' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'crank_nicolson.png'
splot for [i=1:50] '../results/result_01_cranknicolson.dat' using (300*(i-1)*0.1025E+05):1:2 every::(100*(i-1)+1)::(100*i) with points pt 5 ps 2 palette t ''

reset

set terminal png size 1024,1024
set ylabel 'x[m]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set xlabel 't[s]' offset -0.5,0 font 'Times Roman Bold Italic,20'
set zlabel 'T[°C]' offset -0.5,0 font 'Times Roman Bold Italic,20'
set samples 20
set isosamples 21
set cbrange [0:100]
set xyplane at 0
set title 'Distribución de temperaturas: Método Explícito' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'explicit_superf.png'
unset key
set dgrid3d
set contour base
#set samples 20
#set isosamples 21
#set hidden3d
#set cntrlabel start 2 font ",7"
set cntrparam level incremental -3, 0.5, 3
splot for [i=1:50] '../results/result_01_explicit.dat' using (300*(i-1)*0.1025E+05):1:2 every::(100*(i-1)-1)::(100*i) w lines lw 1 palette t ''