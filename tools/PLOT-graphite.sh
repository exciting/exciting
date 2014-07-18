#!/bin/bash

cd output/
files=`ls INFO_*_PBE.OUT`

if [ -e .dummy ]; then 
    rm -f .dummy 
fi

for f in $files; do

    d=${f/INFO_d/}

    echo ${d/_PBE.OUT/} \
         `grep -s "Total energy" $f | tail -1 | awk '{print $4}'` \
         `grep -s "Exchange energy" $f | tail -1 | awk '{print $4}'` \
         `grep -s "Correlation energy" $f | tail -1 | awk '{print $4}'` \
         `grep -s "Exchange energy" ${f/PBE/revPBE} | tail -1 | awk '{print $4}'` \
         `grep -s "Correlation energy" ${f/PBE/LDA} | tail -1 | awk '{print $4}'` \
         `grep -s "Ec_NL=" noloco_d${d/_PBE.OUT/.out} | awk '{print $2}'` \
         >> .dummy

done

echo "# d[A]     E_pbe[meV]    E_c^NL[meV]    E_vdWDF[meV]"
sort -n .dummy | \
awk 'BEGIN{n=0}
     {
        n++; 
        d[n]=$1;
        pbe[n]=$2;
        pbe_x[n]=$3;
        pbe_c[n]=$4;
        rpbe_x[n]=$5;
        lda_c[n]=$6;
        enl[n]=$7;
     }
     END{
        const=27.2*1000/4; # [Ha] -> [meV/atom]
        for(i=1;i<=n;i++){
            printf "%f  %12.8f  %12.8f  %12.8f\n", \
                d[i],
                (pbe[i]-pbe[n])*const,
                (enl[i]-enl[n])*const,
                ((pbe[i]-pbe_x[i]-pbe_c[i]+rpbe_x[i]+lda_c[i]+enl[i])-\
                 (pbe[n]-pbe_x[n]-pbe_c[n]+rpbe_x[n]+lda_c[n]+enl[n]))*const;
        }
     }'

rm -rf .dummy
exit
