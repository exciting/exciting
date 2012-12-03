#!/usr/bin/tcsh -f
# AUTHORS: 
# Monodeep Chakraborty and Juergen Spitaler
#
# DATE:
# Jul. 2, 2010
#
# LICENCE AND AGREEMENTS:
# All Copyright is with the 
# Chair of Atomistic Modelling and Design of Materials
# University of Leoben, 
# Franz-Josef-Strasse 18, 8700 Leoben Austria

set bin = $0:h

set eldir = ../

set latfile = ../lat.in

# Extract elements from $latfile and check,
# if the corresponding .lat.in exist in $eldir

set elements =  ` $bin/extract_elements.csh ../lat.in` 

foreach i ($elements)
  set elfile = $eldir/$i.lat.in
  if (! -e $elfile) then
    printf '>> ERROR in vegard2new.csh: File %s not found!\n' $elfile
    exit
  endif
end

#Check if $latfile really  exists.
if (! -e $latfile ) then
    printf '>> ERROR in vegard2new.csh: File %s not found!\n' $latfile
    exit
endif

awk -v latfile=$latfile -v eldir=$eldir \
    'BEGIN{\
        while(getline<latfile){\
          # Check lat.in file for substitutional elements \
          # and save them in the array "elements". \
          if ($0 ~ /,/){\
          #  print $0; \
            for (i=4;i<=NF;i++) {str=gensub(" ","","g",$i);str_tot=str_tot""str}\
            n_ele=split(str_tot,elements,",")\
          #  for(i=1;i<=n_ele;i++){print elements[i]}\
          }\
        };\
        # Read lattice parameters from element_name.lat.in file corresponding \
        # to each element found in str.out.\
        for(i=1;i<=n_ele;i++){;\
          fname=eldir"/"elements[i]".lat.in"\
        #  print fname\
          getline<fname;ll[i,1,1]=$1;ll[i,1,2]=$2;ll[i,1,3]=$3;\
          getline<fname;ll[i,2,1]=$1;ll[i,2,2]=$2;ll[i,2,3]=$3;\
          getline<fname;ll[i,3,1]=$1;ll[i,3,2]=$2;ll[i,3,3]=$3;\
        }\
        for(i=1;i<=n_ele;i++){;\
        #  print elements[i]\
        #  printf "%10.6f %10.6f %10.6f\n", ll[i,1,1], ll[i,1,2], ll[i,1,3] \
        #  printf "%10.6f %10.6f %10.6f\n", ll[i,2,1], ll[i,2,2], ll[i,2,3] \
        #  printf "%10.6f %10.6f %10.6f\n", ll[i,3,1], ll[i,3,2], ll[i,3,3] \
        };\
        # Get concentration of each substitutional elements in str.out:\
        while(getline<"str.out"){\
          for(i=1;i<=n_ele;i++){\
             if ($0 ~ elements[i]){n[i]++; ntot++}\
          }\
        };\
        # for(i=1;i<=n_ele;i++){print elements[i],n[i]};\
        # print ntot\
        # Now calculate concentration-weighted lattice constants \
        # using Vegards law: \
        for(i=1;i<=n_ele;i++){\
          fac=n[i]/ntot; \
          # print elements[i], "fac = " fac\
          for(i_row=1;i_row<=3;i_row++){\
            for(i_col=1;i_col<=3;i_col++){\
              ll_vegard[i_row,i_col]=ll_vegard[i_row,i_col]+ll[i,i_row,i_col]*fac\
            }\
          }\
        }\
        #Print out results -- lattice vectors from Vegards law:\
        for(i_row=1;i_row<=3;i_row++){\
            printf "%12.8f %12.8f %12.8f\n", ll_vegard[i_row,1], ll_vegard[i_row,2], ll_vegard[i_row,3]\
        }\
      }'\
        < str.out
    


exit
help:
awk '/# *SYNTAX:/,$0 !~ / *#.*[a-zA-Z]+/;\
     /# *EXPLANATION:/,$0 !~ / *#.*[a-zA-Z]+/' < $0

