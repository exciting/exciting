#!/bin/bash
#-------------------------------------------------------------------------------
IXML=input.xml
TMP=$1
len=${#TMP}
if [ "$len" -gt 0 ] ; then IXML=$TMP ; fi
#
if [ -f $IXML ]; then
    xmllint $IXML --schema $EXCITINGROOT/xml/excitinginput.xsd --noout
else
    echo 
    echo " File $IXML do not found!"
    echo 
fi
