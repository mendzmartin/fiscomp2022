#!/usr/bin/gnuplot
set terminal png size 1024,1024
set ylabel 'longitud de la barra [m]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set xlabel 'tiempo [s]' offset -0.5,0 font 'Times Roman Bold Italic,20'
set pm3d map
t_write=2;t_write_total=(10+1);param_t=0.1025E+05
set cbrange [-1:1];set xrange [0:(param_t*t_write*(t_write_total-1))];set yrange [0:1]
#set title 'Distribución de temperaturas: Método Explícito' offset 0,-0.5 font 'Times Roman Bold Italic,20'
#set output 'explicit_vn.png'
#splot for [i=1:t_write_total] '../results/result_01_explicit_vn.dat' using (t_write*(i-1)*param_t):1:2 every::(21*(i-1))::(21*i) with lp lw 8 palette t ''
#set output 'implicit_vn.png'
#set title 'Distribución de temperaturas: Método Implícito' offset 0,-0.5 font 'Times Roman Bold Italic,20'
#splot for [i=1:t_write_total] '../results/result_01_implicit_vn.dat' using (t_write*(i-1)*param_t):1:2 every::(21*(i-1))::(21*i) with lp lw 8 palette t ''
#set output 'crank_nicolson_vn.png'
#set title 'Distribución de temperaturas: Método de Crank-Nicolson' offset 0,-0.5 font 'Times Roman Bold Italic,20'
#splot for [i=1:t_write_total] '../results/result_01_cranknicolson_vn.dat' u (t_write*(i-1)*param_t):1:2 every::(22*(i-1))::(23*i) with lp lw 8 palette t ''
reset