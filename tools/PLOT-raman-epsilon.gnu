#!/bin/bash
#_______________________________________________________________________________


dpisize=100
len=${#DPIPNG}
if [ "$len" -gt 0 ]; then dpisize=$DPIPNG ; fi

if [ -f 'gnu-input' ]; then rm gnu-input; fi

echo 
read -p 'Enter range of plotted displacements in Bohr >>>> ' umin umax
echo

input_list=`ls -d RAMAN_EPSILON_*.OUT`

for input in $input_list ; do
   test_entry=`head -n 20 $input | tail -1 | awk '{print $2}'`
   result=$( echo "$test_entry == 0.0" | bc -l)
   if [ $result -eq 1 ]; then continue ; fi
   suffix=$(echo $input | cut -c15-25)
   oc=$(echo $input | cut -c17-18)
   cat>>gnu-input<<***
    set terminal postscript enhanced color solid lw 3 size 14,7 font 22
    set output 'PLOT_EPSILON_$suffix.ps'
    set multiplot layout 1,2
    set xlabel 'abs. displacement [Bohr]'
    set xrange [$umin:$umax]
    set ylabel 'Re {/Symbol e}^{$oc}({/Symbol w}_L)' offset 0.0
    #set ytics format '%.2f'
    set key spacing 2.0 box top left
    plot 'RAMAN_EPSILON_$suffix.OUT' u 1:2 w l t 'fit', '' u 4:5 w p pt 6 lc 3 t 'data'
    # ---
    set xlabel 'abs. displacement [Bohr]'
    #set xrange [$umin:$umax]
    set ylabel 'Im {/Symbol e}^{$oc}({/Symbol w}_L)' offset 0.0
    #set ytics format '%.3f'
    set key spacing 2.0 box top left
    plot 'RAMAN_EPSILON_$suffix.OUT' u 1:3 w l t 'fit', '' u 4:6 w p pt 6 lc 3 t 'data'
    # ---
    unset multiplot
    quit
***
   gnuplot < gnu-input
   rm gnu-input
   #-----------------------------------------------------------------------
   convert -density $dpisize -rotate 90 PLOT_EPSILON_$suffix.ps PLOT_EPSILON_$suffix.png
done
#-----------------------------------------------------------------------

