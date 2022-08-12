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
set xtics  norangelimit logscale autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit logscale autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit logscale autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit logscale autofreq 
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
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
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
set logscale z 10
set logscale y 10
set logscale x 10
set logscale cb 10
set logscale y2 10
set logscale x2 10
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
x = 0.0
unset logscale
## Last datafile plotted: "md_diffusion_vs_time_T0.75_20.dat"
sizex=8;sizey=8;set terminal pdf size sizex,sizey;set output 'bd_diffusion_vs_time.pdf'
rows=2;columns=2;set multiplot layout rows,columns
    #set xrange[0.8:1.2];set yrange[:3];set xtics 0.05;set ytics 0.25
    set autoscale;set grid
    set title "n_{p}=256,r_{cutoff}=2.5,FCC structure;T_{adim}=0.75\n\
    Brownian Dynamic Simulations\n\
    BD -> t_{eq}=100000,t_{run}=50000,{/Symbol D}t=.001"
    set xlabel "Correlation time (log({/Symbol t}/{t_{0}}))"
    set ylabel "Mean squared displacement (log(MSD))"
    # set xlabel "Correlation time ({/Symbol t})"
    # set ylabel "Mean squared displacement (MSD)"
    set key left
    set grid;set key font ",12";set xlabel  font ",12" ;set ylabel  font ",12"
    set logscale
    D0=0.75/27.0491
    x1label=0.001;y1label=0.0001;set label "diffusion" at x1label,y1label center
    x2label=0.1;y2label=1;set label "diffusion" at x2label,y2label center
    x3label=0.01;y3label=0.01;set label "subdiffusion" at x3label,y3label center
    p '../results/bd_diffusion_vs_time_T0.75_11.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'red' t '{/Symbol r}=0.83',\
    '../results/bd_diffusion_vs_time_T0.75_12.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'blue' t '{/Symbol r}=0.86',\
    '../results/bd_diffusion_vs_time_T0.75_13.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'green' t '{/Symbol r}=0.89',\
    '../results/bd_diffusion_vs_time_T0.75_14.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'magenta' t '{/Symbol r}=0.92',\
    '../results/bd_diffusion_vs_time_T0.75_15.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'brown' t '{/Symbol r}=0.95',\
    '../results/bd_diffusion_vs_time_T0.75_16.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'cyan' t '{/Symbol r}=0.98',\
    '../results/bd_diffusion_vs_time_T0.75_17.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'black' t '{/Symbol r}=1.01',\
    '../results/bd_diffusion_vs_time_T0.75_18.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'dark-blue' t '{/Symbol r}=1.04',\
    '../results/bd_diffusion_vs_time_T0.75_19.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'dark-green' t '{/Symbol r}=1.07',\
    '../results/bd_diffusion_vs_time_T0.75_20.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'orange' t '{/Symbol r}=1.1'

    unset label
    set ylabel "Diffusion constant (log(D)=log(MSD/({/Symbol t}{D_{0}})))"
    set key right top
    p '../results/bd_diffusion_vs_time_T0.75_11.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'red' t '{/Symbol r}=0.83',\
    '../results/bd_diffusion_vs_time_T0.75_12.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'blue' t '{/Symbol r}=0.86',\
    '../results/bd_diffusion_vs_time_T0.75_13.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'green' t '{/Symbol r}=0.89',\
    '../results/bd_diffusion_vs_time_T0.75_14.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'magenta' t '{/Symbol r}=0.92',\
    '../results/bd_diffusion_vs_time_T0.75_15.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'brown' t '{/Symbol r}=0.95',\
    '../results/bd_diffusion_vs_time_T0.75_16.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'cyan' t '{/Symbol r}=0.98',\
    '../results/bd_diffusion_vs_time_T0.75_17.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'black' t '{/Symbol r}=1.01',\
    '../results/bd_diffusion_vs_time_T0.75_18.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'dark-blue' t '{/Symbol r}=1.04',\
    '../results/bd_diffusion_vs_time_T0.75_19.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'dark-green' t '{/Symbol r}=1.07',\
    '../results/bd_diffusion_vs_time_T0.75_20.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'orange' t '{/Symbol r}=1.1'

    set title "n_{p}=256,r_{cutoff}=2.5,FCC structure;T_{adim}=1.35\n\
    Brownian Dynamic Simulations\n\
    BD -> t_{eq}=100000,t_{run}=50000,{/Symbol D}t=.001"
    D0=0.75/27.0491
    x1label=0.001;y1label=0.001;set label "diffusion" at x1label,y1label center
    x2label=0.1;y2label=1;set label "diffusion" at x2label,y2label center
    x3label=0.01;y3label=0.01;set label "subdiffusion" at x3label,y3label center
    set ylabel "Mean squared displacement (log(MSD))"
    set logscale
    set key left top
    p '../results/bd_diffusion_vs_time_T1.35_11.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'red' t '{/Symbol r}=0.83',\
    '../results/bd_diffusion_vs_time_T1.35_12.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'blue' t '{/Symbol r}=0.86',\
    '../results/bd_diffusion_vs_time_T1.35_13.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'green' t '{/Symbol r}=0.89',\
    '../results/bd_diffusion_vs_time_T1.35_14.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'magenta' t '{/Symbol r}=0.92',\
    '../results/bd_diffusion_vs_time_T1.35_15.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'brown' t '{/Symbol r}=0.95',\
    '../results/bd_diffusion_vs_time_T1.35_16.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'cyan' t '{/Symbol r}=0.98',\
    '../results/bd_diffusion_vs_time_T1.35_17.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'black' t '{/Symbol r}=1.01',\
    '../results/bd_diffusion_vs_time_T1.35_18.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'dark-blue' t '{/Symbol r}=1.04',\
    '../results/bd_diffusion_vs_time_T1.35_19.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'dark-green' t '{/Symbol r}=1.07',\
    '../results/bd_diffusion_vs_time_T1.35_20.dat' u ($1*D0):3 w lp pt 7 ps 0.5 lw 0.1 lc 'orange' t '{/Symbol r}=1.1'

    unset label
    set ylabel "Diffusion constant (log(D)=log(MSD/({/Symbol t}{D_{0}})))"
    set key right top
    p '../results/bd_diffusion_vs_time_T1.35_11.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'red' t '{/Symbol r}=0.83',\
    '../results/bd_diffusion_vs_time_T1.35_12.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'blue' t '{/Symbol r}=0.86',\
    '../results/bd_diffusion_vs_time_T1.35_13.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'green' t '{/Symbol r}=0.89',\
    '../results/bd_diffusion_vs_time_T1.35_14.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'magenta' t '{/Symbol r}=0.92',\
    '../results/bd_diffusion_vs_time_T1.35_15.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'brown' t '{/Symbol r}=0.95',\
    '../results/bd_diffusion_vs_time_T1.35_16.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'cyan' t '{/Symbol r}=0.98',\
    '../results/bd_diffusion_vs_time_T1.35_17.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'black' t '{/Symbol r}=1.01',\
    '../results/bd_diffusion_vs_time_T1.35_18.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'dark-blue' t '{/Symbol r}=1.04',\
    '../results/bd_diffusion_vs_time_T1.35_19.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'dark-green' t '{/Symbol r}=1.07',\
    '../results/bd_diffusion_vs_time_T1.35_20.dat' u ($1*D0):(1/($2*6)) w lp pt 7 ps 0.5 lw 0.1 lc 'orange' t '{/Symbol r}=1.1'

unset multiplot
#    EOF
