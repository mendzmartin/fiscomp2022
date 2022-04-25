#!/usr/bin/gnuplot
set terminal png size 1024,1024
set output 'comparison.png'
set xlabel 'position [m]' font 'Times Roman Bold Italic,20'
set ylabel 'temperature [K]' font 'Times Roman Bold Italic,20'
set grid
set xrange [0:1]
set yrange [0:101]
set style line 1 linetype rgb "red" linewidth 3
set style line 2 linetype rgb "green" linewidth 3
set style line 3 linetype rgb "blue" linewidth 3
set style line 4 linetype rgb "orange" linewidth 3
set style line 5 linetype rgb 'black' linewidth 3
set style line 6 linetype rgb 'magenta' linewidth 3
set style line 7 linetype rgb 'purple' linewidth 3
set style line 8 linetype rgb 'cyan' linewidth 3
set title 'Temperature vs position bar' offset 0,-0.5 font 'Times Roman Bold Italic,20'
#smooth unique
plot '../results/result_02_aprox_explicit.dat' using 1:2 every ::0::10 title 'T_{aprox}(x,180[s]-exp)' with linespoints linestyle 1 pointsize 4,\
     '../results/result_02_aprox_implicit.dat' u 1:2 every ::0::10 title 'T_{aprox}(x,180[s]-imp)' with linespoints linestyle 1 pointsize 3,\
     '../results/result_02_aprox_implicit.dat' u 1:2 every ::0::10 title 'T_{aprox}(x,180[s]-imp)' with linespoints linestyle 3 pointsize 2,\
     '../results/result_02_aprox_cranknicolson.dat' using 1:2 every ::0::10 title 'T_{exact}(x,180[s])' with linespoints linestyle 4,\
     '../results/result_02_aprox_explicit.dat' using 1:2 every ::10::20 title 'T_{aprox}(x,1800[s]-exp)' with linespoints linestyle 5 pointsize 4,\
     '../results/result_02_aprox_implicit.dat' using 1:2 every ::10::20 title 'T_{aprox}(x,1800[s]-imp)' with linespoints linestyle 6 pointsize 3,\
     '../results/result_02_aprox_cranknicolson.dat' using 1:2 every ::10::20 title 'T_{exact}(x,1800[s])' with linespoints linestyle 7 pointsize 2,\
     '../results/result_02_exact.dat' using 1:2 every ::10::20 title 'T_{exact}(x,1800[s])' with linespoints linestyle 8 pointsize 1
reset