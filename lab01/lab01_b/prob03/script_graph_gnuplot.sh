#!/usr/bin/gnuplot
set terminal png size 1024,1024
set autoscale x
set autoscale y
set style line 1 linetype rgb "red" linewidth 0.1
set style line 2 linetype rgb "green" linewidth 3
set style line 3 linetype rgb "blue" linewidth 3
set style line 4 linetype rgb "violet" linewidth 3
set style line 5 linetype rgb "black" linewidth 3
set style line 6 linetype rgb "sienna4" linewidth 3
set style line 7 linetype rgb "brown" linewidth 3
set style line 8 linetype rgb "grey" linewidth 3
set style line 9 linetype rgb "cyan" linewidth 3
set key outside
set key font ",20"
set grid 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set output 'trayectory_E06.png'
set xrange[0:150]
set xlabel 't [s]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'theta [m]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'trayectoria para E/(m1*L1^2)=-0.6[s^{-2}], N=2^{21}' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_01.dat' using 1:2 title 'theta 1' with points pointtype 1 pointsize 0.1 lc rgb "red",\
	 'result_01.dat' using 1:3 title 'theta 2' with points pointtype 2 pointsize 0.1 lc rgb "blue"
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set output 'trayectory_E0.png'
set title 'trayectoria para E/(m1*L1^2)=0[s^{-2}], N=2^{21}' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_02.dat' using 1:2 title 'theta 1' with points pointtype 1 pointsize 0.1 lc rgb "red",\
	 'result_02.dat' using 1:3 title 'theta 2' with points pointtype 2 pointsize 0.1 lc rgb "blue"
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set output 'energy_E06.png'
set xlabel 't [s]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'E [J]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'energía para condiciones iniciales 1, N=2^{21}' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_01.dat' using 1:6 notitle with points pointsize 0.1 pointtype 1 lc rgb "red"

set output 'energy_E0.png'
set title 'energía para condiciones iniciales 2, N=2^{21}' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_02.dat' using 1:6 notitle with points pointsize 0.1 pointtype 2 lc rgb "red"
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set logscale y
set autoscale x
set key inside
set xlabel 'omega [rad/s]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'sde' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'spectrum.png'
set yrange[1e-6:]
plot 'result_ft_01.dat' using 1:4 title 'theta 1-E/(m1*L1^2)=-.745[s^{-2}]' smooth unique with linespoints linestyle 1,\
	 'result_ft_02.dat' using 1:4 title 'theta 2-E/(m1*L1^2)=-.745[s^{-2}]' smooth unique with linespoints linestyle 2,\
	 'result_ft_01_E0.dat' using 1:4 title 'theta 1-E/(m1*L1^2)=0[s^{-2}]' smooth unique with linespoints linestyle 3,\
	 'result_ft_02_E0.dat' using 1:4 title 'theta 2-E/(m1*L1^2)=0[s^{-2}]' smooth unique with linespoints linestyle 4,\
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

reset
