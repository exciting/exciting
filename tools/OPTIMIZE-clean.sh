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
    mv ../info.xml .
    mv ../geometry.xml .
    if [ -f ../DFSCFMAX.OUT ];     then mv ../DFSCFMAX.OUT .     ; fi
    if [ -f ../WARNINGS.OUT ];     then mv ../WARNINGS.OUT .     ; fi
    if [ -f ../geometry_opt.xml ]; then mv ../geometry_opt.xml . ; fi
    cd ../
    eta=$(ls -l | grep ^- | wc -l)
    if [ $eta != '0' ]; then 
       ls -l | grep ^- | awk '{print $9}' | xargs rm
    fi
    mv tmp/* .
    rm -r tmp
    cd ../
done
