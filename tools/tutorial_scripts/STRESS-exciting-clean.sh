#!/bin/bash
#

label=`ls -d Dst??`
for Dstn in $label ; do
    cd $Dstn
    Dstn_num_list=`ls -d ${Dstn}_??`
    for Dstn_num in $Dstn_num_list ; do
        cd $Dstn_num/
        mkdir ../tmp
        if [ -f INFO.OUT         ]; then mv INFO.OUT         ../tmp ; fi
        if [ -f info.xml         ]; then mv info.xml         ../tmp ; fi
        if [ -f input.xml        ]; then mv input.xml        ../tmp ; fi
        if [ -f geometry.xml     ]; then mv geometry.xml     ../tmp ; fi
        if [ -f geometry_opt.xml ]; then mv geometry_opt.xml ../tmp ; fi
        if [ -f RMSDVEFF.OUT     ]; then mv RMSDVEFF.OUT     ../tmp ; fi
        if [ -f TOTENERGY.OUT    ]; then mv TOTENERGY.OUT    ../tmp ; fi
        if [ -f WARNINGS.OUT     ]; then mv WARNINGS.OUT     ../tmp ; fi
        rm -Rf ./*
        mv ../tmp/* ./
        rm -R ../tmp/
        cd ../
    done
    cd ../
done
