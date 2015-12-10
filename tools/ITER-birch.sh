#!/bin/bash    
#_______________________________________________________________________________
#
if [ ${#1} -ne 0 ] ; then
    col=$1
else 
    col="2"
fi 
#
emin=`head -n1 energy-vs-strain`
i=0
for xmin in $emin ; do
    if [ $i -eq 0 ] ; then a=$xmin ; fi
    i=1
done
#
emax=`tail -n1 energy-vs-strain`
i=0
for xmax in $emax ; do
    if [ $i -eq 0 ] ; then b=$xmax ; fi
    i=1
done
#
a=`echo "scale=8;0 - $a" | bc -l` 
a=`echo "scale=8;$a * 1.000001" | bc -l` 
b=`echo "scale=8;$b * 1.000001" | bc -l` 
#
list_maxstrain=`seq -s " " 1 50`
#
rm -rf birch-iter
echo
for maxstrain in $list_maxstrain ; do
    x=`echo "scale=2;$maxstrain / 100." | bc -l`
    ra=`echo "scale=0;$x / $a" | bc -l`
    rb=`echo "scale=0;$x / $b" | bc -l`
    #echo $ra, $rb, $x, $maxstrain
    if [[ "$ra" -eq "0"  ||  "$rb" -eq "0" ]] ; then 
              echo " maxstrain = "$maxstrain"%"
              PLOT-birch.py -maxstrain $x > birch-dum 2>&1
              parameters=$(tail -n3 birch-dum | head -n1)
              slen=${#parameters}
              if [ "$slen" -ne "0" ] ; then echo $x $parameters >> birch-iter ; fi 
    fi 
done
#
PLOT-column.py birch-iter 1 $col
rm -rf birch-dum
echo
