#!/bin/bash
#_______________________________________________________________________________
#
tempfile="temp"
workdir="workdir"
prefix="rundir-"
number_of_leading_zeros=1
#
current_number_of_dirs=$(ls -d $prefix* 2>>$tempfile|wc -l)
number_of_new_dir=$(echo "$current_number_of_dirs + 1"|bc)
#
for i in $(seq $number_of_leading_zeros) ; do
    if [ $number_of_new_dir -lt $(echo "10^$i"|bc) ] ; then
       number_of_new_dir="0$number_of_new_dir"
    fi
done
#
eta=$number_of_new_dir
tot=workdir/INFO.OUT
OUT=energy-vs-step
touch $OUT
awk -v eta="$eta" '/Total energy    / {printf "%5i   %20.10f\n",eta,$4}' $tot | tail -n1>>$OUT
energy_list=$(cat $OUT | tail -n1)
echo $energy_list > dum ; energy=$(cat dum | awk '{ printf "%20.8f", $2 }') ; rm -f dum
echo $energy_list > dum ;  niter=$(cat dum | awk '{ printf "%3i", $1 }')    ; rm -f dum
#
parameters=$(LATTICE-parameters.sh $workdir/input.xml | tail -n2 | head -n1)
#
PLOT-energy.py >/dev/null 2>&1
#
if [ -f initial-step ] ; then
   rm -f initial-step
   energy_zero=$energy
   echo "$energy_zero" > energy-zero
fi
#
energy_zero=$(cat energy-zero | awk '{ printf "%20.8f", $1 }')
echo "$energy" "$energy_zero" > dum
delta=$(awk '/ / {printf "%10.6f", $1-$2}' dum) ; rm -f dum
#
echo "$parameters" "$delta"
#
cp $workdir/input.xml ./input_opt.xml
mv $workdir $prefix$number_of_new_dir 2>>$tempfile
rm $tempfile