#!/bin/bash
printf "\n################################################"
printf "\n#            WANNIER FUNCTION PLOT             #"
printf "\n################################################\n\n"

read -p "Which Wannier function do you want to plot? Enter index here. " i
n=$(printf "%04d" $i)

s='a'
while [ "$s" != "s" -a "$s" != "r" -a "$s" != "i" ]; do
    read -p "What do you want to plot, square modulus (s), real part (r) or imaginary part (i)? " s
done

if [ "$s" == "s" ]; then
    bash wannierplot_smod.sh $n > plot.xcrysden
elif [ "$s" == "i" ]; then
    bash wannierplot_imag.sh $n > plot.xcrysden
else
    bash wannierplot_real.sh $n > plot.xcrysden
fi
xcrysden -s plot.xcrysden >/dev/null 2>&1 &
