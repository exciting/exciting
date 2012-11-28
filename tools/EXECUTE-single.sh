#!/bin/bash    
#_______________________________________________________________________________
#
EXECUTABLE=$EXCITINGROOT/bin/excitingser
#
CURRENT=$PWD
#-------------------------------------------------------------------------------
WORKDIR=rundir-00
if [ ${#1} -gt 0 ]; then WORKDIR=$1 ; fi
if [ -d "$WORKDIR" ]; then rm -Rf $WORKDIR ; fi
#
RUNDIR=$CURRENT
#
if [ -d "$EXCITINGRUNDIR" ]; then
   RUNDIR=$EXCITINGRUNDIR
else 
   if [ ${#EXCITINGRUNDIR} -gt 0 ]; then
      mkdir $EXCITINGRUNDIR
      RUNDIR=$EXCITINGRUNDIR
   fi
fi 
#
XCRUNDIR=xc-rundir
i=0
while [ -d "$RUNDIR/$XCRUNDIR" ]; do 
   i=$(($i + 1))
   XCRUNDIR="$XCRUNDIR$i"
done
#-------------------------------------------------------------------------------
echo
echo "===> Output directory is \""$WORKDIR"\" <==="
echo
#-------------------------------------------------------------------------------
echo $EXECUTABLE > $CURRENT/exciting
#
cd $RUNDIR
#
echo
echo "Running exciting for file input.xml -------------------------------------"
echo
#
rm -Rf $XCRUNDIR
mkdir  $XCRUNDIR
cd     $XCRUNDIR
cp $CURRENT/input.xml input.xml
#
time $EXECUTABLE | tee output.screen
#
mv $RUNDIR/$XCRUNDIR $CURRENT/$WORKDIR
echo
echo "Run completed for file input.xml ----------------------------------------"
echo
echo
#

