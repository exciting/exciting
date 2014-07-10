#!/bin/bash
#-------------------------------------------------------------------------------

dpisize=100
len=${#DPIPNG}
if [ "$len" -gt 0 ]; then dpisize=$DPIPNG ; fi

dollar="$"
if [ -f 'gnu-input' ] ; then rm gnu-input ; fi
inpfile=convergence-test

#-------------------------------------------------------------------------------

sed '/^$/d' $inpfile > gdum
tail -n1 gdum > idum
awk '{ print $3 }' idum > rdum
eta=$(cat rdum)
awk '{printf "%20.3f\n",sqrt($1*$1)}' rdum > gdum
awk '{ print $1 }' gdum > edum
ene=$(cat edum)
awk -v eta="$eta" \
        '/ / {printf "%4i  %6.2f  %20.10f\n",$1,$2,$3-eta}' $inpfile > idum

#-------------------------------------------------------------------------------
cat>>gnu-input<<***
 set ter pos landscape enhanced color solid lw 2
 set out 'PLOT.ps'
 set sty data l
 set multiplot 
 set origin  0.05, -0.01

 set format z "%10.3f"
 set xtics 2
 set xtics 2 offset -0.5, -0.5 font "Helvetica,18"
 set ytics 1 offset  0.5, -0.5 font "Helvetica,18"
 set ticslevel 0.1

# set ztick 0.05
# set palette rgbformulae 22,13,-31
# set pm3d at b

 set title "Convergence test for $1" offset -3.0, -0.8 font 'Helvetica, 25'
 set xlabel "ngridk" offset 0.0, -0.8 font 'Helvetica, 22'
 set ylabel "rgkmax" offset 0.0, -0.6 font 'Helvetica, 22'

 set label "Energy+$ene [Ha]" at screen 0.08, screen 0.45 \
           font 'Helvetica, 22' rotate by +90 center

 set view 60, 120, 1, 0.94

#set xrange [3:20]
#set yrange [6.0:9.0]
#set zrange [-5312.345:-5312.435]

 unset key
 set dgrid3d $2,$3
 set hidden3d

# set contour base
# set isosample 20

 splot "idum" u 1:2:3 

# set multiplot   
# set sty data lp
#set xtics 0.01
# set mxtics 0
# set grid 
# set key o r
# set key spacing 1.5 
# set key box
#===============================================
# set size    0.86, 0.83
# set origin  0.05, 0.05
# set yr [$yyymin:$yyymax]
# set xr [:]
#set logscale y
#set logscale x
# set title  "$ttitle" offset 0,-0.8  font "Helvetica,18"
# set xlabel ""        offset 0.0,0.0 font "Helvetica,18"     
# set ylabel "$ylabel" offset 0.0,0.0 font "Helvetica,18"
# set xtics                           font "Helvetica,16"
# set ytics                           font "Helvetica,16"
# plot '$1' u 1:$coloumn pt 7 ps $pointsize lt 1 lw 2 title ' 1'    
#===============================================
# unset multiplot  
 quit
***
gnuplot < gnu-input
rm -f gnu-input idum rdum edum gdum
#-------------------------------------------------------------------------------
convert -density $dpisize -rotate 90 PLOT.ps PLOT.png
#-------------------------------------------------------------------------------
