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
set autoscale
set xlabel 'q-coordenada' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set ylabel 'p-momento' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set title 'diagrama de fases E=5' offset 0,-0.5 font 'Times Roman Bold Italic,20'
set output 'trayectory_E1.png'
plot 'results_E1.dat' using 1:3 title 'orbita 1 (red)' with points pointtype 1 pointsize 0.1 lc rgb "red",\
	 'results_E1.dat' using 2:4 title 'orbita 2 (blue)' with points pointtype 2 pointsize 0.1 lc rgb "blue"
set output 'trayectory_E2.png'
set title 'diagrama de fases E=20' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'results_E2.dat' using 1:3 title 'orbita 1 (red)' with points pointtype 1 pointsize 0.1 lc rgb "red",\
	 'results_E2.dat' using 2:4 title 'orbita 2 (blue)' with points pointtype 2 pointsize 0.1 lc rgb "blue"
set output 'trayectory_E3.png'
set title 'diagrama de fases E=100' offset 0,-0.5 font 'Times Roman Bold Italic,20'
plot 'results_E3.dat' using 1:3 title 'orbita 1 (red)' with points pointtype 1 pointsize 0.1 lc rgb "red",\
	 'results_E3.dat' using 2:4 title 'orbita 2 (blue)' with points pointtype 2 pointsize 0.1 lc rgb "blue"
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
