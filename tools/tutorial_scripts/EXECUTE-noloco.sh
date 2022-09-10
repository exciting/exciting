#!/bin/bash

vdwdf_path="$EXCITINGROOT/src/src_vdwdf"
cd output
xsf_files=`ls RHO3D_*.xsf`

for file in $xsf_files; do

cat > input.in << EOF
# Integration parameters
  $vdwdf_path            ! path to vdW-DF kernel 
  vdW-DF                 ! vdW-DF type (vdW-DF, vdW-DF2, VV09)
  4 4 4                  ! number of unit cells in x/y/z-direction
  1.d-16                 ! relative integration accuracy
  1.d-04                 ! absolute integration accuracy in Hartree
  10000                  ! number of Monte-Carlo sampling points
  1000000000             ! maximum number of function evaluations

# Density data file      
  bohr                   ! bohr/angs: charge density units [e/(bohr^3) or e/(angs^3)]
  ${file}                ! density
EOF

name=${file/RHO3D/noloco}
$vdwdf_path/noloco < input.in > ${name/.xsf/.out}

done

rm -f input.in
exit
