#!/bin/csh
if ( "x"$1 == "x-m" ) then
  set machinefile=$2
else
  set machinefile="~/.machines.rc"
endif

if ( ! -e $machinefile ) then
  echo Unable to read $machinefile
  if (! -e ~/.machines.rc) then
    cat - <<EOF
I have created a ~/.machines.rc file for you.
Please edit it to reflect you machine setup.
EOF
    cat - >! ~/.machines.rc <<EOF
#configuration file for chl, minload, and pollmach
#this line indicates the waiting time between the checks of machine availability
#in seconds
set waitbetweenpoll=60
#The remainder of this file lists the machines
#columns in this file are separated by a +
#each line corresponds to a machine
#the first column indicates the command to obtain the load on a remote machine
#the second column is the command prefix to lauch a job on that remote machine
# note that the command must cd into the same directory on the remote machine
# as on the local machine. The "node" command does that automatically.
# Type "node" on the command line for more info.
# 
#the first line sets the threshold load for not starting a job (here 0.5)
#do not remove the "none" keyword
echo 0.5 + none
#list the machines here
#example on a local network (e.g. beowolf) with shared disk
rsh machinename uptime | getvalue average + node machinename
#secure version of the above
ssh -n user@machinename uptime | getvalue average + node -s user@machinename  
#if the machines do not share the same disk
ssh -n user@machinename uptime | getvalue average + node -s user@machinename  
#example for queueing system
llq -u $USER | wc -l | awk '{print ($1-4)/2}' + poe 
EOF
  endif
  exit
endif

cat $machinefile | grep -v '^#' | grep -v '^set ' | awk -F + '{printf $2 " + "; system($1)}' | awk -F + '{print $2+0 " + " $1}'
