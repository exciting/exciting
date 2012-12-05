#!/bin/bash
#

label=`ls -d Dst??`
for Dstn in $label ; do
    cd $Dstn
    Dstn_num_list=`ls -d ${Dstn}_??`
    for Dstn_num in $Dstn_num_list ; do
        cd $Dstn_num/
        mkdir tmp/
        cd tmp/
        mv ../INFO.OUT .
        mv ../info.xml .
        mv ../input.xml .
        mv ../RMSDVEFF.OUT .
        mv ../GEOMETRY.OUT .
        mv ../TOTENERGY.OUT .
        mv ../DTOTENERGY.OUT .
        cd ../
        ls -l | grep ^- | awk '{print $9}' | xargs rm
        mv tmp/* .
        rm -r tmp/
        cd ../
    done
    cd ../
done
