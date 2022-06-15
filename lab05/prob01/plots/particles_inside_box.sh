#!/usr/bin/bash
index=0
for i in {10..20}
    do
        gnuplot -e "set terminal jpeg;
        set title 'molecular dynamics - steady state';
        set xlabel 'x coordinate' rotate parallel;
        set ylabel 'y coordinate' rotate parallel;
        set zlabel 'z coordinate' rotate parallel;
        splot '../results/picture${i}.dat' u 1:2:3 w p pt 7 ps 0.2 lc 'red' notitle" > ../results/movie${index}.jpeg
        index=$((index+1))
    done

ffmpeg -r 9 -i ../results/movie%d.jpeg movie.mp4

for i in {0..10}
    do 
    rm -f ../results/movie${i}.jpeg
    done