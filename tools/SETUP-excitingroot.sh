#!/bin/bash
#-------------------------------------------------------------------------------
IXML=input.xml
TMP=$1
len=${#TMP}
if [ "$len" -gt 0 ] ; then IXML=$TMP ; fi
#
if [ -f $IXML ]; then
    awk -v xcr=$EXCITINGROOT '{ gsub(/\$EXCITINGROOT/, xcr); print }' $IXML > tmpfile
    mv tmpfile $IXML
else
    echo 
    echo " File $IXML do not found!"
    echo 
fi
