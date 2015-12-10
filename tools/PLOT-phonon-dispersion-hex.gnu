#!/bin/bash
#_______________________________________________________________________________

pointsize=0.01

dpisize=100
len=${#DPIPNG}
if [ "$len" -gt 0 ]; then dpisize=$DPIPNG ; fi

inpf='PHDISP.OUT'

alat=`awk -F"\"" '/scale/{print $2}' input.xml`

#echo 
#read -p 'Enter lattice parameter in Bohr >>>> ' alat
#echo

awk -v a="$alat"\
    'NF==2 {printf "%11.8f   %20.10f\n",$1*a,$2*2.194746313705e5};\
     NF==0 {printf("\n")}' $inpf > gnu-file

if [ -f 'gnu-input' ]; then rm gnu-input; fi

xk=4.18879016 
xx=6.28318531
xg=9.91078404

y0=$1
yy=$2
if [ "$y0" = "" ]; then y0=0;   fi
if [ "$yy" = "" ]; then yy=2600; fi

cat>>gnu-input<<***
 set ter pos landscape enhanced color solid lw 2
 set out 'PLOT.ps'
 set multiplot   
 set sty data lp
#set xtics 0.02
 set mxtics 0
 set grid 
 set xtics ("{/Symbol G}" 0, "K" $xk, "M" $xx,"{/Symbol G}" $xg)
 set arrow from $xk,$y0 to $xk,$yy nohead lt 0 lw 2
 set arrow from $xx,$y0 to $xx,$yy nohead lt 7 lw 2
 set arrow from $xg,$y0 to $xg,$yy nohead lt 7 lw 2
 set key o r
 set key spacing 1.5 
 set key box
 set nokey

#===============================================
 set size    0.86, 0.83
 set origin  0.05, 0.05
 set yr [$1:$2]
 set xr [0:$xg]
 set title "Lattice parameter = $alat Bohr"  offset 0,-0.3  font "Helvetica,22"
 set xlabel ""   offset 0.0,0.0 font "Helvetica,1"     
 set ylabel "Frequency [cm^-^1]" offset -0.5,0.0 font "Helvetica,22"
 set xtics font "Helvetica,18"
 set ytics font "Helvetica,18"
 plot 'gnu-file' u 1:2 pt 7 ps $pointsize lt 1 lw 2 title '' 
#===============================================
 unset multiplot  
 quit
***
gnuplot < gnu-input
rm gnu-input #gnu-file
#-----------------------------------------------------------------------
convert -density $dpisize -rotate 90 PLOT.ps PLOT.png
#-----------------------------------------------------------------------

