#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 8    last modified 2019-12-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2019
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal qt 0 font "Sans,9"
# set output
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "% h" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics back
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key fixed right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set view azimuth 0
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels 5
set cntrparam levels auto
set cntrparam firstlinetype 0 unsorted
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set title "" 
set title  font "" textcolor lt -1 norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" textcolor lt -1 norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 21.5645 : 869.140 ] noreverse writeback
set x2range [ 22.5429 : 869.271 ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ -0.0543657 : 0.871606 ] noreverse writeback
set y2range [ 0.0369184 : 0.867239 ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "es_AR.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
GNUTERM = "qt"
## Last datafile plotted: "result_01e_20x20_autocorr_T2.0.dat"
set terminal pdf size 4,4;set output 'autocorrelation1.pdf'
set multiplot layout 1,1     
    set xlabel "{/Symbol t}_{corr}(MC_{step})"
    unset grid;set key font ",12";set xlabel  font ",16" ;set ylabel  font ",16"
    set ylabel "A|U_{adim}({/Symbol t}_{corr})/A|M_{adim}({/Symbol t}_{corr})"
    set yrange[0:1]
    set title "Autocorrelación para T_{adim}/T_{adim}^c=0.88"
    f(x)=a*exp(-x/b)
    fit [0:1000] f(x) '../results/result_01e_20x20_autocorr_T2.0.dat' u 1:3 via a,b
    p '../results/result_01e_10x10_autocorr_T2.0.dat' u 1:2 w l lw 1 lc rgb 'red' t 'A|U_{adim},laticce(l)=10x10',\
    '../results/result_01e_10x10_autocorr_T2.0.dat' u 1:3 w l lw 1 lc rgb 'magenta' t 'A|M_{adim},l=10x10',\
    '../results/result_01e_20x20_autocorr_T2.0.dat' u 1:2 w l lw 1 lc rgb 'green' t 'A|U_{adim},l=20x20',\
    '../results/result_01e_20x20_autocorr_T2.0.dat' u 1:3 w l lw 1 lc rgb 'brown' t 'A|M_{adim},l=20x20',\
    f(x) w l lw 1 dt 2 lc 'black' t 'fit:{/Symbol t}_{corr}=600',\
    '../results/result_01e_20x20_autocorr_T2.0.dat' u (600):0 w l lw 1 dt 3 lc 'black' notitle
    set origin 0.17,0.28;set size 0.5,0.5;set xrange [0:300]; set yrange [0:1]; set xtics 100
    set logscale y;set nokey;set grid
    unset title;set ylabel 'ZOOM:logscale(eje y)';unset xlabel
    replot '../results/result_01e_10x10_autocorr_T2.0.dat' every 4 u 1:2 w l lw 1 lc rgb 'red' smooth mcsplines notitle,\
    '../results/result_01e_10x10_autocorr_T2.0.dat' every 4 u 1:3 w l lw 1 lc rgb 'magenta' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T2.0.dat' every 4 u 1:2 w l lw 1 lc rgb 'green' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T2.0.dat' every 4 u 1:3 w l lw 1 lc rgb 'brown' smooth mcsplines notitle,\
    f(x) w l lw 1 dt 2 lc 'black'
    reset
unset multiplot

set terminal pdf size 4,4;set output 'autocorrelation2.pdf'
set multiplot layout 1,1     
    set xlabel "{/Symbol t}_{corr}(MC_{step})"
    unset grid;set key font ",16";set xlabel  font ",16" ;set ylabel  font ",16"
    set ylabel "A|U_{adim}({/Symbol t}_{corr})/A|M_{adim}({/Symbol t}_{corr})"
    set title "Autocorrelación para T_{adim}/T_{adim}^c=0.98"
    set yrange[0:2]
    f(x)=a*exp(-x/b)
    fit [0:1000] f(x) '../results/result_01e_20x20_autocorr_T2.22.dat' u 1:3 via a,b
    p '../results/result_01e_10x10_autocorr_T2.22.dat' u 1:2 w l lw 1 lc rgb 'red' t 'A|U_{adim},laticce=10x10',\
    '../results/result_01e_10x10_autocorr_T2.22.dat' u 1:3 w l lw 1 lc rgb 'magenta' t 'A|M_{adim},laticce=10x10',\
    '../results/result_01e_20x20_autocorr_T2.22.dat' u 1:2 w l lw 1 lc rgb 'green' t 'A|U_{adim},laticce=20x20',\
    '../results/result_01e_20x20_autocorr_T2.22.dat' u 1:3 w l lw 1 lc rgb 'brown' t 'A|M_{adim},laticce=20x20',\
    '../results/result_01e_20x20_autocorr_T2.22.dat' u (200):0 w l lw 1 dt 3 lc 'black' notitle,\
    f(x) w l lw 1 dt 2 lc 'black' t 'fit:{/Symbol t}_{corr}=200'
    set origin 0.21,0.2;set size 0.5,0.5;set xrange [0:100]; set yrange [0:1]
    set logscale y;set nokey;set grid
    unset title;set ylabel 'ZOOM:logscale(eje y)';unset xlabel
    replot '../results/result_01e_10x10_autocorr_T2.22.dat' every 4 u 1:2 w l lw 1 lc rgb 'red' smooth mcsplines notitle,\
    '../results/result_01e_10x10_autocorr_T2.22.dat' every 4 u 1:3 w l lw 1 lc rgb 'magenta' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T2.22.dat' every 4 u 1:2 w l lw 1 lc rgb 'green' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T2.22.dat' every 4 u 1:3 w l lw 1 lc rgb 'brown' smooth mcsplines notitle,\
    f(x) w l lw 1 dt 2 lc 'black'
    reset
unset multiplot

set terminal pdf size 4,4;set output 'autocorrelation3.pdf'
set multiplot layout 1,1     
    set xlabel "{/Symbol t}_{corr}(MC_{step})"
    unset grid;set key font ",12";set xlabel  font ",16" ;set ylabel  font ",16"
    set ylabel "A|U_{adim}({/Symbol t}_{corr})/A|M_{adim}({/Symbol t}_{corr})"
    set title "Autocorrelación para T_{adim}/T_{adim}^c=0.99"
    set yrange[0:1]
    f(x)=a*exp(-x/b)
    fit [0:1000] f(x) '../results/result_01e_40x40_autocorr_T2.2676.dat' u 1:3 via a,b
    p '../results/result_01e_10x10_autocorr_T2.2676.dat' u 1:2 w l lw 1 lc rgb 'red' t 'A|U_{adim},laticce=10x10',\
    '../results/result_01e_10x10_autocorr_T2.2676.dat' u 1:3 w l lw 1 lc rgb 'magenta' t 'A|M_{adim},laticce=10x10',\
    '../results/result_01e_20x20_autocorr_T2.2676.dat' u 1:2 w l lw 1 lc rgb 'green' t 'A|U_{adim},laticce=20x20',\
    '../results/result_01e_20x20_autocorr_T2.2676.dat' u 1:3 w l lw 1 lc rgb 'brown' t 'A|M_{adim},laticce=20x20',\
    '../results/result_01e_40x40_autocorr_T2.2676.dat' u 1:2 w l lw 1 lc rgb 'blue' t 'A|U_{adim},laticce=40x40',\
    '../results/result_01e_40x40_autocorr_T2.2676.dat' u 1:3 w l lw 1 lc rgb 'cyan' t 'A|M_{adim},laticce=40x40',\
    f(x) w l lw 1 dt 2 lc 'black' t 'fit:{/Symbol t}_{corr}>1000'
    set origin 0.35,0.19;set size 0.5,0.5;set xrange [0:400]; set yrange [0:1]; set xtics 100
    set logscale y;set nokey;set grid
    unset title;set ylabel 'ZOOM:logscale(eje y)';unset xlabel
    replot '../results/result_01e_10x10_autocorr_T2.2676.dat' every 4 u 1:2 w l lw 1 lc rgb 'red' smooth mcsplines notitle,\
    '../results/result_01e_10x10_autocorr_T2.2676.dat' every 4 u 1:3 w l lw 1 lc rgb 'magenta' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T2.2676.dat' every 4 u 1:2 w l lw 1 lc rgb 'green' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T2.2676.dat' every 4 u 1:3 w l lw 1 lc rgb 'brown' smooth mcsplines notitle,\
    '../results/result_01e_40x40_autocorr_T2.2676.dat' every 4 u 1:2 w l lw 1 lc rgb 'blue' smooth mcsplines notitle,\
    '../results/result_01e_40x40_autocorr_T2.2676.dat' every 4 u 1:3 w l lw 1 lc rgb 'cyan' smooth mcsplines notitle,\
    f(x) w l lw 1 dt 2 lc 'black'
    reset
unset multiplot

set terminal pdf size 4,4;set output 'autocorrelation4.pdf'
set multiplot layout 1,1     
    set xlabel "{/Symbol t}_{corr}(MC_{step})"
    unset grid;set key font ",12";set xlabel  font ",16" ;set ylabel  font ",16"
    set ylabel "A|U_{adim}({/Symbol t}_{corr})/A|M_{adim}({/Symbol t}_{corr})"
    set title "Autocorrelación para T_{adim}/T_{adim}^c=1.4"
    set yrange[0:1]
    f(x)=a*exp(-x/b)
    fit [0:1000] f(x) '../results/result_01e_20x20_autocorr_T3.3.dat' u 1:3 via a,b
    p '../results/result_01e_10x10_autocorr_T3.3.dat' u 1:2 w l lw 1 lc rgb 'red' t 'A|U_{adim},laticce=10x10',\
    '../results/result_01e_10x10_autocorr_T3.3.dat' u 1:3 w l lw 1 lc rgb 'magenta' t 'A|M_{adim},laticce=10x10',\
    '../results/result_01e_20x20_autocorr_T3.3.dat' u 1:2 w l lw 1 lc rgb 'green' t 'A|U_{adim},laticce=20x20',\
    '../results/result_01e_20x20_autocorr_T3.3.dat' u 1:3 w l lw 1 lc rgb 'brown' t 'A|M_{adim},laticce=20x20',\
    '../results/result_01e_20x20_autocorr_T3.3.dat' u (160):0 w l lw 1 dt 3 lc 'black' notitle,\
    f(x) w l lw 1 dt 2 lc 'black' t 'fit:{/Symbol t}_{corr}=160'
    set origin 0.35,0.19;set size 0.5,0.5;set xrange [0:60]; set yrange [0:1]
    set logscale y;set nokey;set grid
    unset title;set ylabel 'ZOOM:logscale(eje y)';unset xlabel
    replot '../results/result_01e_10x10_autocorr_T3.3.dat' every 4 u 1:2 w l lw 1 lc rgb 'red' smooth mcsplines notitle,\
    '../results/result_01e_10x10_autocorr_T3.3.dat' every 4 u 1:3 w l lw 1 lc rgb 'magenta' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T3.3.dat' every 4 u 1:2 w l lw 1 lc rgb 'green' smooth mcsplines notitle,\
    '../results/result_01e_20x20_autocorr_T3.3.dat' every 4 u 1:3 w l lw 1 lc rgb 'brown' smooth mcsplines notitle,\
    f(x) w l lw 1 dt 2 lc 'black'
    reset
unset multiplot
