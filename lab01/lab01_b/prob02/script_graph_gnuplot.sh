#!/usr/bin/gnuplot

set terminal png size 1024,1024

# generate de axis titles
set xlabel 't' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'x(t)' offset 0,-0.5 font 'Times Roman Bold Italic,20'

set xrange[0:50]
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
set style line 7 linetype rgb "brown" linewidth 3
set style line 8 linetype rgb "grey" linewidth 3
set style line 9 linetype rgb "cyan" linewidth 3

set key outside
set key font ",20"
set grid 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set output 'logistic_map_x06_1.png'
set title 'mapa logístico x(1)=0.6' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_logistic_map.dat' using 1:2 every::1::512 title 'r=1.5' smooth unique with linespoints linestyle 1,\
		'result_logistic_map.dat' using 1:2 every::2049::2560 title 'r=3.3' smooth unique with linespoints linestyle 2,\
		'result_logistic_map.dat' using 1:2 every::4097::4608 title 'r=3.5' smooth unique with linespoints linestyle 3,\
		'result_logistic_map.dat' using 1:2 every::6145::6656 title 'r=3.55' smooth unique with linespoints linestyle 4,\

set terminal png size 1024,1024
set output 'logistic_map_x06_2.png'
set xrange[0:200]
set title 'mapa logístico x(1)=0.6' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_logistic_map.dat' using 1:2 every::1::512 title 'r=1.5' smooth unique with linespoints linestyle 1,\
		'result_logistic_map.dat' using 1:2 every::8193::8704 title 'r=4' smooth unique with linespoints linestyle 2
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set terminal png size 1024,1024
set output 'fourier_transform_1.png'
set autoscale x
set logscale y
set yrange[1e-05:1]
set xlabel 'frecuencia angular [rad/s]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'transformada de Fourier - mapa logístico x(1)=0.6' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_ft_1.dat' using 1:4 title 'r=1.5' smooth unique with linespoints linestyle 1,\
	 'result_ft_2.dat' using 1:4 title 'r=3.3' smooth unique with linespoints linestyle 2,\
	 'result_ft_3.dat' using 1:4 title 'r=3.5' smooth unique with linespoints linestyle 3,\
	 'result_ft_4.dat' using 1:4 title 'r=3.55.' smooth unique with linespoints linestyle 4,\
	 'result_ft_5.dat' using 1:4 title 'r=4' smooth unique with linespoints linestyle 6

set output 'fourier_transform_01.png'
set ylabel 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'transformada de Fourier - mapa logístico x(1)=0.6 & r=1.5' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_ft_1.dat' using 1:4 notitle smooth unique with linespoints linestyle 1

set output 'fourier_transform_02.png'
set ylabel 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'transformada de Fourier - mapa logístico x(1)=0.6 & r=3.3' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_ft_2.dat' using 1:4 notitle smooth unique with linespoints linestyle 1

set output 'fourier_transform_03.png'
set ylabel 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'transformada de Fourier - mapa logístico x(1)=0.6 & r=3.5' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_ft_3.dat' using 1:4 notitle smooth unique with linespoints linestyle 1

set output 'fourier_transform_04.png'
set ylabel 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'transformada de Fourier - mapa logístico x(1)=0.6 & r=3.55' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_ft_4.dat' using 1:4 notitle smooth unique with linespoints linestyle 1

set output 'fourier_transform_05.png'
set ylabel 'densidad espectral de energía [Js/rad]' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'transformada de Fourier - mapa logístico x(1)=0.6 & r=4' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'result_ft_5.dat' using 1:4 notitle smooth unique with linespoints linestyle 1
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set autoscale x
set autoscale y
set xlabel 'x' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'counter' offset 0,-0.5 font 'Times Roman Bold Italic,20'

set title 'normalized histogram-logistic map x(1)=0.6 & r=4' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'histogram.png'
plot 'histogram_01.dat' using 1:2 title 'N° bins=10' with boxes,\
	 'histogram_02.dat' using 1:2 title 'N° bins=100' with boxes,\
	 'histogram_03.dat' using 1:2 title 'N° bins=1000' with boxes
set title 'histograma normalizado - mapa logístico x(1)=0.6 & r=3.3' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'histogram_r33.png'
set style fill solid 0.5
plot 'histogram.dat' using 1:2 title 'N° bins=100' with boxes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set autoscale x
set autoscale y
unset logscale y
set xlabel 'r' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'x' offset 0,-0.5 font 'Times Roman Bold Italic,20'

set title 'mapa logistico: mapa de bifurcación' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'bif_graph_01.png'
plot 'bif_graph_01.dat' using 1:2 notitle with points pointsize 0.1

set yrange[0.13:0.174]
set yrange[3.847:3.8568]
set output 'bif_graph_02.png'
plot 'bif_graph_02.dat' using 1:2 notitle with points pointsize 0.1

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set autoscale
set xlabel 'r^*' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'lambda' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'lyapunov_exp.png'
set title 'exponentes de lyapunov - mapa logístico x(1)=0.6' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'lyapunov_exp.dat' using 1:2 notitle smooth unique with linespoints linestyle 1

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# reset all of graph characteristics to default values
reset
