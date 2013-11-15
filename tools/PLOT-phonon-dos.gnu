#!/bin/bash
#_______________________________________________________________________________

pointsize=0.01

dpisize=100
len=${#DPIPNG}
if [ "$len" -gt 0 ]; then dpisize=$DPIPNG ; fi

inpf='PHDOS.OUT'
xlab='Frequency [cm^-^1]'
title=''

awk \
    'NF==2 {printf "%20.10f   %20.10f\n",$1*2.194746313705e5,$2/2.194746313705e5};\
     NF==0 {printf("\n")}' $inpf > gnu-file

if [ -f 'gnu-input' ]; then rm gnu-input; fi

cat>>gnu-input<<***
 set terminal postscript landscape enhanced color solid linewidth 2
 set out 'PLOT.ps'
 set multiplot   
 set style data linespoints
#set xtics 0.02
 set mxtics 0
 set grid 
 set key o r
 set key spacing 1.5 
 set key box
 set nokey
#===============================================
 set size    0.86, 0.83
 set origin  0.05, 0.05
 set yr [:]
 set xr [:]
 set format y "%6.3f"
 set title "$title" offset 0,-0.8 font "Helvetica,25"
 set xlabel "$xlab" offset 0.0,-0.7 font "Helvetica,22"     
 set ylabel "Phonon DOS [states/cm^-^1]" offset 0.0,0.0 font "Helvetica,22"
 set xtics font "Helvetica,16"
 set ytics font "Helvetica,16"
 plot 'gnu-file' u 1:2 pt 7 ps $pointsize lt 1 lw 2 title ''    
#===============================================
 unset multiplot  
 quit
***
gnuplot < gnu-input
rm gnu-input gnu-file
#-------------------------------------------------------------------------------
convert -density $dpisize -rotate 90 PLOT.ps PLOT.png
#-------------------------------------------------------------------------------
