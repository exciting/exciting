#!/bin/bash
#_______________________________________________________________________________

ul='_'

dpisize=100
len=${#DPIPNG}
if [ "$len" -gt 0 ]; then dpisize=$DPIPNG ; fi

if [ -f 'gnu-input' ]; then rm gnu-input; fi

input_list=`ls -d RAMAN_RESONANCE_*.OUT`

for input in $input_list ; do
   test_entry=`head -n 20 $input | tail -1 | awk '{print $5}'`
   result=$( echo "$test_entry == 0.0" | bc -l)
   if [ $result -eq 1 ]; then continue ; fi
   suffixoc=$(echo $input | cut -c17-20)
   suffixmod=$(echo $input | cut -c22-27)
   cat>>gnu-input<<***
    set terminal postscript landscape enhanced color solid lw 2 font 22
    set output 'PLOT_RESONANCE_$suffixoc$ul$suffixmod.ps'
    set xlabel 'laser energy [eV]'
    set xrange [1.0:4.1]
    set xtics format '%.2f'
    set logscale y 10
    set ylabel '|d{/Symbol c}{\251}/dQ|^2 [au]' offset 2.0
    set key spacing 2.0 box top left
    plot 'RAMAN_RESONANCE_$suffixoc$ul$suffixmod.OUT' u 2:5 w l t '  $suffixoc{\137}$suffixmod'
    quit
***
   gnuplot < gnu-input
   rm gnu-input
   #-----------------------------------------------------------------------
   convert -density $dpisize -rotate 90 PLOT_RESONANCE_$suffixoc$ul$suffixmod.ps PLOT_RESONANCE_$suffixoc$ul$suffixmod.png
done
#-----------------------------------------------------------------------

