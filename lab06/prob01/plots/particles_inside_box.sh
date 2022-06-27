#!/usr/bin/bash
index=0
for i in {10..100}
    do
        gnuplot -e "set terminal jpeg;
        set title 'brownian dynamics - steady state - L=8.5499';
        set view 65,35;
        set xrange[-5:5];set yrange[-5:5];set zrange[-5:5];
        set xyplane at -6; L=8.5499;
        set xlabel 'x coordinate' rotate parallel;
        set ylabel 'y coordinate' rotate parallel;
        set zlabel 'z coordinate' rotate parallel;
        splot '../results/picture${i}.dat' u 1:2:3 w p pt 7 ps 0.2 lc 'dark-red' notitle,\
        '../results/fcc.dat' u 1:(-L/2):(-L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u 1:(L/2):(L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u 1:(-L/2):(L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u 1:(L/2):(-L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (-L/2):2:(-L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (L/2):2:(L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (-L/2):2:(L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (L/2):2:(-L/2) w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (-L/2):(-L/2):3 w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (L/2):(L/2):3 w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (-L/2):(L/2):3 w l lw 0.2 lc 'black' notitle,\
    '../results/fcc.dat' u (L/2):(-L/2):3 w l lw 0.2 lc 'black' notitle" > ../results/movie${index}.jpeg
        index=$((index+1))
    done


ffmpeg -r 9 -i ../results/movie%d.jpeg movie.mp4

for i in {0..100}
    do 
    rm -f ../results/movie${i}.jpeg
    done