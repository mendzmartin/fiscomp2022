#!/usr/bin/gnuplot
set terminal png size 1024,1024
set xlabel 't_{adm}' offset 0,-1 font 'Times Roman Bold Italic,20'
set ylabel 'E_{abs}' offset 1,0 font 'Times Roman Bold Italic,20'
set xrange[0:1];set yrange[0:0.0026]
set grid
set output 'errores_absolutos_deltaT0001.png'
set title 'Error absoluto vs tiempo adimensional' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set key font ",20"
plot '../results/result_02_explicit_vn_err.dat' u 1:2 w lp lc rgb 'red' pt 2 title 'explícito',\
     '../results/result_02_implicit_vn_err.dat' u 1:2 w lp lc rgb 'blue' pt 3 title 'implícito',\
     '../results/result_02_cranknicolson_vn_err.dat' u 1:2 w lp lc rgb 'green' pt 4 title 'cranknicolson'
reset