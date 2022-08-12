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
set xrange [ -3.66291 : 33.8371 ] noreverse writeback
set x2range [ -3.66291 : 33.8371 ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ -6.60604 : -1.60604 ] noreverse writeback
set y2range [ -6.48713 : -1.85579 ] noreverse writeback
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
## Last datafile plotted: "md_energies_without_linked_list.dat"
sizey=12;sizex=4;set terminal pdf size sizex,sizey;set output 'energies_vs_time.pdf'
rows=3;columns=1;set multiplot layout rows,columns
    set ylabel "Potential adimensional energy (U_{adim}/{n_{p}})"
    set key right
    set grid;set key font ",12";set xlabel  font ",12" ;set ylabel  font ",12"
    set parametric
    t1=3000;set trange [-6:-5.3]

    set autoscale;set xrange[0:10000];set yrange[-6:-5.2]
    set xlabel "MD_{step} (Molecular Dynamic)"
    set title "n_{p}=256,r_{cutoff}=2.5,FCC structure,T_{adim}=0.75\n\
    {/Symbol r}=.8,{/Symbol D}t=.005;MD_{step}=10000"
    x1label=1500;y1label=-5.6;set label "transitory" at x1label,y1label center
    x2label=6000;y2label=-5.6;set label "estationary" at x2label,y2label center
    p '../results/md_energies_with_linkedlist_rho08.dat' u 0:(($2/256)-0.4198) w l lw 2 lc 'red' t 'with linked-list Uadim_{med}=-5.34+-0.14' smooth mcsplines,\
    '../results/md_energies_without_linkedlist_rho08.dat' u 0:(($2/256)-0.4198) w l lw 2 lc 'blue' t 'without linked-list Uadim_{med}=-5.35+-0.15' smooth mcsplines,\
    t1,t w l lw 3 lc 'black' dt 2 notitle

    set autoscale;set xrange[0:100000];set yrange[-6:-1.5]
    t2=80000;set trange [-6:-5]
    set xlabel "BD_{step} (Brownian Dynamic)"
    set title "n_{p}=256,r_{cutoff}=2.5,FCC structure,T_{adim}=0.75\n\
    {/Symbol r}=.8,{/Symbol D}t=.001;BD_{step}=100000;{BD_{step}}^{ens}=10"
    unset label
    x1label=40000;y1label=-4;set label "transitory" at x1label,y1label center
    x2label=90000;y2label=-4;set label "estationary" at x2label,y2label center
    p '../results/bd_energies_with_linkedlist_rho08.dat' u 1:(($2/256)-0.4198) w l lw 2 lc 'blue' t 'with linked-list Uadim_{med}=-5.34+-0.64' smooth mcsplines,\
    '../results/bd_energies_without_linkedlist_rho08.dat' u 1:(($2/256)-0.4198) w l lw 2 lc 'red' t 'without linked-list Uadim_{med}=-5.37+-0.52' smooth mcsplines,\
    t2,t w l lw 3 lc 'black' dt 2 notitle

    set autoscale;set xrange[0:10000];set yrange[-6:-5.2]
    t3=6000;set trange [-6:-5.3]
    set xlabel "MCD_{step} (Monte Carlo Dynamic)"
    set title "n_{p}=256,r_{cutoff}=2.5,FCC structure,T_{adim}=0.75\n\
    {/Symbol r}=.8,MCD_{step}=10000"
    unset label
    x1label=3000;y1label=-5.6;set label "transitory" at x1label,y1label center
    x2label=8000;y2label=-5.6;set label "estationary" at x2label,y2label center
    p '../results/mcd_energies_without_linkedlist_rho08.dat' u 1:(($2/256)-0.4198) w l lw 2 lc 'red' t 'without linked-list Uadim_{med}=-5.36+-0.19' smooth mcsplines,\
    t3,t w l lw 3 lc 'black' dt 2 notitle
unset multiplot


# results
# MD
#   +++++++++++ without linked-list +++++++++++
#     cpu_time      delta_t     r_cutoff        T_ref          n_p         t_eq        t_run  tau_max_corr
#   0.1682E+03   0.5000E-02   0.2500E+01   0.7500E+00          256            0        10000          100
#   Uadim_med=  -4.9272859003801894      +-  0.14657077219363179
#   +++++++++++ with linked-list +++++++++++ 
#     cpu_time      delta_t     r_cutoff        T_ref          n_p         t_eq        t_run  tau_max_corr
#   0.1521E+03   0.5000E-02   0.2500E+01   0.7500E+00          256            0        10000          100
#   Uadim_med=  -4.9212789648311270      +-  0.14319123660881800 
# BD
#   +++++++++++ without linked-list +++++++++++
#     cpu_time      delta_t     r_cutoff        T_ref          n_p         t_eq        t_run tau_max_corr
#   0.2473E+04   0.1000E-02   0.2500E+01   0.7500E+00          256            0       100000          100
#    Uadim_med=  -4.9473788928258715      +-  0.52574288609920783
#   +++++++++++ with linked-list +++++++++++
#     cpu_time      delta_t     r_cutoff        T_ref          n_p         t_eq        t_run tau_max_corr         t_ens
#   0.2134E+04   0.1000E-02   0.2500E+01   0.7500E+00          256            0       100000          100
#    Uadim_med=  -4.9167994132768884      +-  0.64471896882402180    
# MCD
#   +++++++++++ without linked-list +++++++++++
#     cpu_time      delta_t     r_cutoff        T_ref          n_p   MC_step_eq  MC_step_run  tau_max_corr
#   0.2267E+03   0.5000E-01   0.2500E+01   0.7500E+00          256            0        10000          100
#   Uadim_med=  -4.9441419509924955      +-  0.19483034310410302
#   +++++++++++ with linked-list +++++++++++
#         
