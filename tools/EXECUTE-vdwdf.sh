#!/bin/bash

# Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
# This file is distributed under the terms of the GNU Lesser General Public
# License. See the file COPYING for license details.

#====================================================
#
# This script shows how to calculate the vdW-DF total 
# energy using the EXCITING code
#
#====================================================

#-------------------------------
# 3D electron density grid (system dependent)
# the grid should be large enough to describe correctly
# the total density distribution and reliably calculate density
# gradients
#-------------------------------

nrho_x="40"
nrho_y="40"
nrho_z="60"

#-------------------------------
# NOLOCO integration parameters (system dependent)
# require EcNL convergence tests
#-------------------------------

# the supercell (integration volume) size
nx="1"
ny="1"
nz="1"

# the absolute integration accuracy (in Hartree)
eabs="0.001"

#====================================================
#   Input files and commands
#====================================================

exciting="$EXCITINGROOT/bin/excitingser"
xml2xsf="xsltproc $EXCITINGROOT/xml/visualizationtemplates/plot3d2xsf.xsl"
noloco="$EXCITINGROOT/src/vdwdf/noloco"

if (($#<1)); then
    echo "USAGE: $0 <input.xml>"
    exit
fi

# create a temporary work directory
if [ -d vdwdf ]; then
    rm -rf vdwdf
fi
mkdir vdwdf
cp $1 vdwdf/

if [ ! -e STATE.OUT ]; then
    echo "ERROR: The restart file STATE.OUT is not found!"
    echo "Please perform first a self-consitent cycle using PBE-GGA!"
    exit 1
fi
cp STATE.OUT vdwdf/
cd vdwdf

#====================================================
#   Post-GGA calculations
#====================================================
#
# Evdw-df = Epbe - Epbe_xc + (Erevpbe_x + Elda_c) + E_nl
#
sed -e 's/fromscratch/fromfile/g' \
    -e 's/maxscl="200"/maxscl="1"/g' -i $1
    
# To evaluate first 4 terms and the total electron density
sed -e '/xctype/d' \
    -e '/<libxc/,/\/>/d' \
    -e '/<properties>/,/<\/properties>/d' \
    -i $1

awk -v nx=$nrho_x -v ny=$nrho_y -v nz=$nrho_z \
  '{if($0!~"</groundstate>"){print $0}\
      else{printf "\
    <libxc\n\
      exchange=\"XC_GGA_X_PBE_R\"\n\
      correlation=\"XC_LDA_C_PW\"\n\
    />\n\
  </groundstate>\n\
  <properties>\n\
    <chargedensityplot>\n\
      <plot3d>\n\
        <box grid=\"%d %d %d\">\n\
          <origin coord=\"0 0 0\" />\n\
          <point coord=\"1 0 0\"/>\n\
          <point coord=\"0 1 0\"/>\n\
          <point coord=\"0 0 1\"/>\n\
        </box>\n\
      </plot3d>\n\
    </chargedensityplot>\n\
  </properties>\n", nx, ny, nz}}' $1 > .dummy
  
mv .dummy $1

# Run the Exciting
$exciting

cp INFO.OUT INFO_VDWDF.OUT

# for vdW-DF calculations one needs to save the 3d charge density
$xml2xsf RHO3d.xml > RHO3d.xsf

#====================================================
#   Calculation of Ec_NL
#====================================================
cat > noloco.in << EOF
# Integration parameters
  $EXCITINGROOT/src/vdwdf/
  vdW-DF                   ! vdW-DF type (vdW-DF, vdW-DF2, VV09)
  $nx $ny $nz              ! number of unit cells in x/y/z-direction
  1.d-16                   ! relative integration accuracy
  $eabs                    ! absolute integration accuracy in Hartree
  10000                    ! number of Monte-Carlo sampling points
  1000000000               ! maximum number of function evaluations

# Density data file      
  bohr                     ! bohr/angs: charge density units [e/(bohr^3) or e/(angs^3)]
  RHO3d.xsf
EOF

$noloco < noloco.in > EcNL.OUT

#====================================================
#   The final result
#====================================================
ex=`tail -60 INFO.OUT | grep "exchange" | awk '{print $3}'`
ec=`tail -60 INFO.OUT | grep "correlation" | awk '{print $3}'`
ecnl=`grep "Ec_NL=" EcNL.OUT | awk '{print $2}'`

cat >> INFO_VDWDF.OUT << EOF

+------------------------------+
|           NOLOCO             |
+------------------------------+

Energies :
   Ex_revPBE         : $ex
   Ec_LDA            : $ec
   Ec_NL             : $ecnl
   Exc_vdW-DF        : `echo "$ex+$ec+$ecnl" | bc -l`

Timings (CPU) :
`grep "Execution time:" EcNL.OUT`

EOF

cp INFO_VDWDF.OUT ..
exit
