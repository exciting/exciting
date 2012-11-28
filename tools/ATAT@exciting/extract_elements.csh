#!/usr/bin/tcsh -f
# AUTHORS: 
# Juergen Spitaler
#
# DATE:
# Oct. 25, 2010
#
# LICENCE AND AGREEMENTS:
# All Copyright is with the 
# Chair of Atomistic Modelling and Design of Materials
# University of Leoben, 
# Franz-Josef-Strasse 18, 8700 Leoben Austria
#
# SYNTAX:
# extract_elements.csh lat.in-file
#
# EXPLANATION:
# Extracts all substitutional elements from the 
# lat.in-file given as arguments, and writes
# them to standard output line by line.

if (  "$1" == "-h" || "$1" == "-H" ) then
    goto help
endif

set latfile = $1

if (! -e $latfile) then
  printf '>> ERROR in extract_elements.csh: File %s not found!\n' $latfile
  goto help
endif

awk '/,/{\
   for (i=4;i<=NF;i++) {str=gensub(" ","","g",$i);str_tot=str_tot""str}\
   n_ele=split(str_tot,elements,",")\
   for(i=1;i<=n_ele;i++){print elements[i]}\
 }' $latfile


exit

help:
awk '/# *SYNTAX:/,$0 !~ / *#.*[a-zA-Z]+/;\
     /# *EXPLANATION:/,$0 !~ / *#.*[a-zA-Z]+/' < $0

