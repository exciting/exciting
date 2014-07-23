#!/bin/bash
#-------------------------------------------------------------------------------
IXML=input.xml
TMP=$1
len=${#TMP}
if [ "$len" -gt 0 ] ; then IXML=$TMP ; fi
#
xmllint $IXML --schema $EXCITINGROOT/xml/excitinginput.xsd --noout
#