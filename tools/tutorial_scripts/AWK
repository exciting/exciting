#!/bin/bash
#-------------------------------------------------------------------------------
IXML=input.xml
TMP=$1
len=${#TMP}
if [ "$len" -gt 0 ] ; then IXML=$TMP ; fi
awk -v xcr=$EXCITINGROOT '{ gsub(/\$EXCITINGROOT/, xcr); print }' $IXML > tmpfile
mv tmpfile $IXML

