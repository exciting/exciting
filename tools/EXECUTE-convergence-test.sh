#!/bin/bash    
#_______________________________________________________________________________
#
EXECUTABLE=$EXCITINGROOT/bin/excitingser
#
CURRENT=$PWD
#
WORKDIR=workdir
if [ ${#1} -gt 0 ]; then WORKDIR=$1 ; fi
if [ ! -d "$WORKDIR" ]; then 
echo ; echo "Working directory \""$WORKDIR"\" does NOT exist!!" ; echo ; exit
fi
#-------------------------------------------------------------------------------
echo
echo "===> Output directory is \""$WORKDIR"\" <==="
echo
#-------------------------------------------------------------------------------
echo $EXECUTABLE > $CURRENT/$WORKDIR/exciting
#
OUT=$CURRENT/$WORKDIR/convergence-test
#
if [ -f $OUT ] ; then mv $OUT $OUT.save ; fi
#
touch $OUT
#
cd $WORKDIR
#
input_list=`ls -d ngridk*`
#
for input in $input_list ; do
    echo
    echo "Running exciting in directory" $input "-----------"
    echo
    cd $input
    time $EXECUTABLE | tee output.screen
    echo
    echo "Run completed for directory" $input "-------------"
    echo
    nk=$(cat ngridk)
    rg=$(cat rgkmax)
    awk -v nk="$nk" -v rg="$rg"\
        '/  / {printf "%4i  %6.2f  %18.8f\n",nk,rg,$1}' TOTENERGY.OUT |\
        tail -n1>>$OUT
    cd ..
done
echo