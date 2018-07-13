#!/bin/bash    
#_______________________________________________________________________________
#
EXECUTABLE=$EXCITINGROOT/bin/excitingser
#
CURRENT=$PWD
#-------------------------------------------------------------------------------

WORKDIRS=`ls -d */`
for wdir in $WORKDIRS ; do
    #echo $wdir
    if [ -f $wdir/phrun ] ; then
       cd $wdir 
       echo $EXECUTABLE > ./exciting
#
       echo
       echo "*  Running exciting phdisp & phdos in" \"$wdir\"
#
       time $EXECUTABLE > output.screen
       cd ../
    else
       echo
       echo "*  Skipping directory" \"$wdir\"
    fi
done

echo

#ThermoX.py 

#-------------------------------------------------------------------------------