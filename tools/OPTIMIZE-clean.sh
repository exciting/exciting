#!/bin/bash
#

label=`ls -d *_??`
for dirn in $label ; do
    cd $dirn/
    mkdir tmp
    cd tmp/
    mv ../$dirn.xml .
    mv ../INFO.OUT .
    mv ../RMSDVEFF.OUT .
    mv ../GEOMETRY.OUT .
    mv ../TOTENERGY.OUT .
    mv ../info.xml .
    cd ../
    ls -l | grep ^- | awk '{print $9}' | xargs rm
    mv tmp/* .
    rm -r tmp
    cd ../
done
