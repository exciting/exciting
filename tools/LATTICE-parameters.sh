#!/bin/bash
#_______________________________________________________________________________
#
IXML=input.xml
TMP=$1
len=${#TMP}
if [ "$len" -gt 0 ] ; then IXML=$TMP ; fi
#
$EXCITINGTOOLS/exciting2sgroup.py $IXML sgroup.in 
sgroup sgroup.in 1>sgroup.out 2>sgroup.err 
#
sym=$(cat sgroup.out| head -n1) 
abc=$(cat sgroup.out| head -n4 | tail -n1) 
ABG=$(cat sgroup.out| head -n6 | tail -n1)
#
echo $ABG > dum ; A=$(cat dum | awk '{ printf "%10.5f", $1 }') ; rm -f dum
echo $ABG > dum ; B=$(cat dum | awk '{ printf "%10.5f", $2 }') ; rm -f dum
echo $ABG > dum ; G=$(cat dum | awk '{ printf "%10.5f", $3 }') ; rm -f dum
#
echo $abc > dum ; a=$(cat dum | awk '{ printf "%12.5f", $1 }') ; rm -f dum
echo $abc > dum ; b=$(cat dum | awk '{ printf "%10.5f", $2 }') ; rm -f dum
echo $abc > dum ; c=$(cat dum | awk '{ printf "%10.5f", $3 }') ; rm -f dum
#
echo
echo "    " $sym
echo "     a          b          c          alpha      beta       gamma"
echo "$a" "$b" "$c" "$A" "$B" "$G" 
echo
rm -f sgroup.* 
#-------------------------------------------------------------------------------