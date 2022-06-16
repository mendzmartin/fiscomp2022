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
set xrange [ -2.35860 : 2.64140 ] noreverse writeback
set x2range [ -2.15389 : 2.29236 ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ -2.50491 : 34.9951 ] noreverse writeback
set y2range [ -2.32372 : 32.4638 ] noreverse writeback
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
set locale "en_GB.UTF-8"
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
x = 0.0
## Last datafile plotted: "components_velocities_histogram.dat"

set terminal pdf size 8,8;set output 'histograms.pdf'

set multiplot layout 2,2
    set grid
    set autoscale
    set title "n_{p}=256,{/Symbol r}t=.8,T_{adim}=1.1"
    set grid;set key font ",12";set xlabel  font ",12" ;set ylabel  font ",12"

    set xlabel "velocity components (v_{x},v_{y},v_{z})"
    set ylabel "probability (p(v_{x}),p(v_{y}),p(v_{z}))"
    p '../results/components_velocities_histogram.dat' u 1:($2/25) w boxes t 'p(v_{x})',\
    '../results/components_velocities_histogram.dat' u 3:($4/25) w boxes t 'p(v_{y})',\
    '../results/components_velocities_histogram.dat' u 5:($6/25) w boxes t 'p(v_{z})',\
    '../results/exact_velocities_distributions.dat' u 1:2 w l lc 'black' smooth mcsplines t 'exact'

    set xlabel "total velocity (v_{tot})"
    set ylabel "probability (p(v_{tot})"
    p '../results/total_velocities_histogram.dat' u 1:($2/20) w boxes t 'p(v_{tot})',\
    '../results/exact_velocities_distributions.dat' u 7:8 w l lc 'black' smooth mcsplines t 'exact'


    reset

    set title "n_{p}=256,{/Symbol r}t=.8,T_{adim}=1.1"
    set style data histogram
    binwidth=0.2
    bin(x,width)=width*floor(x/width)

    set xlabel "velocity components (v_{x},v_{y},v_{z})"
    set ylabel "probability (p(v_{x}),p(v_{y}),p(v_{z}))"
    p '../results/exact_velocities_distributions.dat' u (bin($1,binwidth)):(1.0) smooth freq with boxes t 'p(v_{x})',\
    '../results/exact_velocities_distributions.dat' u (bin($3,binwidth)):(1.0) smooth freq with boxes t 'p(v_{y})',\
    '../results/exact_velocities_distributions.dat' u (bin($5,binwidth)):(1.0) smooth freq with boxes t 'p(v_{z})'

    binwidth=0.07
    set xlabel "total velocity (v_{tot})"
    set ylabel "probability (p(v_{tot})"
    p '../results/exact_velocities_distributions.dat' u (bin($7,binwidth)):(1.0) smooth freq with boxes t 'p(v_{tot})'

unset multiplot
#    EOF
