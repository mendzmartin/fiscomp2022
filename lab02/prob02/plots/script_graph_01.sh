#!/usr/bin/gnuplot
set terminal png size 1024,1024
set ylabel 'longitud de la barra [m]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set xlabel 'tiempo [s]' offset -0.5,0 font 'Times Roman Bold Italic,20'
set pm3d map
set cbrange [-1:1]
set xrange [0:205000]
set yrange [0:1]

set title 'Distribución de temperaturas: Método Explícito' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'explicit_vn.png'
splot for [i=1:10] '../results/result_01_explicit_vn.dat' using (2*(i-1)*0.1025E+05):1:2 every::(21*(i-1))::(21*i) with l lw 8 palette t '' #points pt 5 ps 2 palette t ''
reset