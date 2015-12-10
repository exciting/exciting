#!/bin/bash    
#_______________________________________________________________________________
#
EXECUTABLE=$EXCITINGROOT/bin/excitingser
#
CURRENT=$PWD
#-------------------------------------------------------------------------------
WORKDIR=workdir
if [ ${#1} -gt 0 ]; then WORKDIR=$1 ; fi
if [ ! -d "$WORKDIR" ]; then 
echo ; echo "Working directory \""$WORKDIR"\" does NOT exist!!" ; echo ; exit
fi
#
RUNDIR=$CURRENT/$WORKDIR
#
if [ -d "$EXCITINGRUNDIR" ]; then
   RUNDIR=$EXCITINGRUNDIR
else 
   len=${#EXCITINGRUNDIR}
   if [ "$len" -gt 0 ]; then
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
echo $EXECUTABLE > $CURRENT/$WORKDIR/exciting
cp $CURRENT/input.xml $CURRENT/$WORKDIR/source.xml
#
if [ -f $CURRENT/planar ] ; then cp $CURRENT/planar $CURRENT/$WORKDIR ; fi
#
if [ -f $CURRENT/$WORKDIR/dft-0.5 ] ; then
    OUT=$CURRENT/$WORKDIR/bandgap-vs-rcut
else
    OUT=$CURRENT/$WORKDIR/energy-vs-strain
fi
#
if [ -f $OUT ] ; then mv $OUT $OUT.save ; fi
#
cd $WORKDIR 
input_list=`ls -d input-*`
cd $RUNDIR
#
aloop=0
#
for input in $input_list ; do
    echo
    echo "Running exciting for file" $input "----------------------------------"
    echo
#
    rm -Rf $XCRUNDIR
    mkdir  $XCRUNDIR
    cd     $XCRUNDIR
    cp $CURRENT/$WORKDIR/$input input.xml
#
    time $EXECUTABLE | tee output.screen
#
    suffix=$(echo $input | cut -c7-8)
#
    eta=$(cat $CURRENT/$WORKDIR/strain-$suffix)
    tot=INFO.OUT
#
    if [ -f $CURRENT/$WORKDIR/dft-0.5 ] ; then
        awk -v eta="$eta" \
        '/Estimated fundam/ {printf "%11.8f   %20.10f\n",eta,$5*27.211396132}' $tot | tail -n1>>$OUT
    else
        awk -v eta="$eta" \
        '/Total energy    / {printf "%11.8f   %20.10f\n",eta,$4}' $tot | tail -n1>>$OUT
    fi
#
    cd ../
#
    rm -Rf $CURRENT/$WORKDIR/rundir-$suffix
    mv $XCRUNDIR $CURRENT/$WORKDIR/rundir-$suffix
    cd $CURRENT/$WORKDIR
#
    if [ -f $CURRENT/$WORKDIR/dft-0.5 ] ; then
        PLOT-plot.py bandgap-vs-rcut "o-" "rcut [Bohr]" "Bandgap [eV]" 
    else    
        if [ $aloop = 1 ]; then PLOT-energy.py ; fi
    fi
#
    cd $RUNDIR
    echo
    echo "Run completed for file" $input "-------------------------------------"
    echo
    aloop=1
#
done
echo

