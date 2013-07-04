#!/bin/bash

if (($#<1)); then
    echo "USAGE: $0 <cif file>"
    exit
fi

HMs=`grep "_symmetry_space_group_name_H-M" $1 | \
     awk 'BEGIN{FS="\047"};{print $2}' | sed -e 's/ //g'`
a=`grep "_cell_length_a" $1 | awk '{printf "%.3f", $2}'`
b=`grep "_cell_length_b" $1 | awk '{printf "%.3f", $2}'`
c=`grep "_cell_length_c" $1 | awk '{printf "%.3f", $2}'`
alpha=`grep "_cell_angle_alpha" $1 | awk '{printf "%.3f", $2}'`
beta=`grep "_cell_angle_beta" $1 | awk '{printf "%.3f", $2}'`
gamma=`grep "_cell_angle_gamma" $1 | awk '{printf "%.3f", $2}'`

# Extract coordinates ( fancy thing )
awk 'BEGIN{n=0};
  /loop_/{getline;
          while($0~/_atom_site/){
              n++;
              if($0~/_atom_site_label/){na=n};
              if($0~/_atom_site_fract_x/){nx=n};
              if($0~/_atom_site_fract_y/){ny=n};
              if($0~/_atom_site_fract_z/){nz=n};
              getline;
          }
          if(n>0){
              while(NF>=n){
                  label=$(na);
                  gsub(/[0-9]/,"",label);
                  printf "%s  %12.10f  %12.10f  %12.10f\n", label, $(nx), $(ny), $(nz);
                  getline;
              }
              exit;
          }
  }' $1 > dummy

# Sort the list by atom name
sort -b -i -d dummy > coords
rm -f dummy

awk 'BEGIN{label="";first=1}
    { 
      if($1!=label){
        if(!first){print  "              </wspecies>"}; first=0;
        printf "              <wspecies speciesfile=\"%s.xml\">\n", $1;
        printf "                <wpos coord=\"%12.10f  %12.10f  %12.10f\"/>\n", \
              $2, $3, $4;
        label=$1;
      }
      else{
        printf "                <wpos coord=\"%12.10f  %12.10f  %12.10f\"/>\n", \
                $2, $3, $4;
      }
    }
    END{
      print  "              </wspecies>";
    }' coords > dummy

cat > spacegroup.xml << EOF
<?xml version="1.0" encoding="UTF-8" ?>
<symmetries HermannMauguinSymbol="$HMs">
    <title>$1</title>
        <lattice scale="1.889727" a="$a" b="$b" c="$c"
           bc="${alpha}" ac="${beta}" ab="${gamma}" 
           ncell="1 1 1" primcell="true"/>
        <WyckoffPositions>
`cat dummy`
        </WyckoffPositions>
</symmetries>
EOF

rm -f coords dummy
exit
