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
    mv ../TOTENERGY.OUT .
    mv ../WARNINGS.OUT .
    mv ../info.xml .
    mv ../geometry.xml .
    mv -f ../geometry_opt.xml .
    cd ../
    ls -l | grep ^- | awk '{print $9}' | xargs rm
    mv tmp/* .
    rm -r tmp
    cd ../
done
