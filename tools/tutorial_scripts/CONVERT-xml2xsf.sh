#!/bin/bash
#-------------------------------------------------------------------------------
#
EXCITINGVISUAL=$EXCITINGROOT/xml/visualizationtemplates
EXCITINGCONVERT=$EXCITINGROOT/xml/inputfileconverter
#
if test "$#" -le 1 ; then
   echo
   echo " Illegal number of arguments!"
   echo
   exit
fi
#
if test "$#" -ge 2 ; then
   XMLFILE=$2
   FILEXMLROOT=$(basename $XMLFILE .xml)
   if test "$#" -ge 3 ; then
      XSFFILE=$3
      FILEXSFROOT=$(basename $XSFFILE .xsf)
   else 
      FILEXSFROOT=$FILEXMLROOT
   fi
fi
#
nd="$(tr [A-Z] [a-z] <<< "$1")"
ext=2xsf
PLOT=plot$nd$ext
#
if [ -f $XMLFILE ] || [-f $FILEXMLROOT.xml] ; then
    if [ "$FILEXMLROOT" == "input" ] ; then
        xsltproc $EXCITINGCONVERT/xmlinput2xsf.xsl input.xml > input.xsf
    else
        xsltproc $EXCITINGVISUAL/$PLOT.xsl $FILEXMLROOT.xml > $FILEXSFROOT.xsf
    fi
else
    echo 
    echo " File $FILEXROOT.xml do not found!"
    echo 
fi
#
#-------------------------------------------------------------------------------
