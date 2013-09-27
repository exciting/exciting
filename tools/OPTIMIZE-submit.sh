#!/bin/bash
#
EXECUTABLE=$EXCITINGROOT/bin/excitingser

label=`ls -d *_??`
for dirn in $label ; do
    cd $dirn
    cp -f $dirn.xml input.xml
    echo
    echo 'SCF calculation of "'$dirn'" starts --------------------------------'
    time $EXECUTABLE | tee output.screen
    cd ../
done
echo 
