#!/bin/bash
#_______________________________________________________________________________

ul='_'

dpisize=100
len=${#DPIPNG}
if [ "$len" -gt 0 ]; then dpisize=$DPIPNG ; fi

if [ -f 'gnu-input' ]; then rm gnu-input; fi

input_list=`ls -d RAMAN_SPEC_*.OUT`

for input in $input_list ; do
   test_entry=`head -n 100 $input | tail -1 | awk '{print $2}'`
   result=$( echo "$test_entry == 0.00000" | bc -l)
   if [ $result -eq 1 ]; then continue ; fi
   suffixoc=$(echo $input | cut -c12-15)
   suffixmod=$(echo $input | cut -c17-22)
   cat>>gnu-input<<***
    set terminal postscript landscape enhanced color solid lw 2 font 22
    set output 'PLOT_SPEC_$suffixoc$ul$suffixmod.ps'
    set xlabel 'Raman shift [cm^{-1}]'
    set ylabel 'S [10^{-5} m^{-1} sr^{-1}]' offset 2.0
    set key spacing 2.0 box top right
    plot 'RAMAN_SPEC_$suffixoc$ul$suffixmod.OUT' u 1:2 w l t '  $suffixoc{\137}$suffixmod'
    quit
***
   gnuplot < gnu-input
   rm gnu-input
   #-----------------------------------------------------------------------
   convert -density $dpisize -rotate 90 PLOT_SPEC_$suffixoc$ul$suffixmod.ps PLOT_SPEC_$suffixoc$ul$suffixmod.png
done
#-----------------------------------------------------------------------

