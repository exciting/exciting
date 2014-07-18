#!/bin/bash
#
# Graphite:
#   Energy vs Interlayer Distance Test
#

#====================================================
#   Input files and commands
#====================================================

# Main structure and calculational parameters
model="input_graphite.xml"

# local (species and xsl paths)
sed -e 's=\$EXCITINGROOT='"$EXCITINGROOT"'=g' $model > draft.xml

# aliases for commands
exciting="$EXCITINGROOT/bin/excitingser"
xml2xsf="xsltproc $EXCITINGROOT/xml/visualizationtemplates/plot3d2xsf.xsl"

#====================================================
#   Run calculations
#====================================================

# directory to store results
results="output"
if [ -e ${results} ]; then 
    rm -rf ${results} 
fi
mkdir ${results}

#----------------------------------------------------

# Distance in Angstrom
distance="3.2 3.4 3.6 3.8 4.0 5.0 6.0 8.0"

# main loop
for d in $distance; do

# crystal c-length
c=$(echo "2.*$d" | bc -l | xargs printf "%5.3f")

# n_z - number of grid points along c-axis ( supposed, for a=2.461, n=14)
nz=$(echo "24*$c/2.4" | bc -l | xargs printf "%1.0f")

# Step 1: Generate the exciting input file
sed -e 's/XXXXX/'"$c"'/g' -e 's/YYYYY/'"$nz"'/g' draft.xml > input.xml

# Step 2: Run the Exciting
$exciting

# Step 3: Save results
cp INFO.OUT ${results}/INFO_d${d}_PBE.OUT

# for vdW-DF calculations one needs to save the 3d charge density
$xml2xsf RHO3D.xml > ${results}/RHO3D_d$d.xsf

# Step 4: post GGA run ( EvdW-df = Epbe - Epbe_xc + (Erevpbe_x + Eldac) + E_nl )
sed -e 's/fromscratch/fromfile/g' \
    -e '/<properties>/,/<\/properties>/d' \
    -e 's/maxscl="200"/maxscl="1"/g' -i input.xml

# revPBE run
sed -e 's/GGA_PBE/GGA_PBE_R/g' -i input.xml
$exciting
cp INFO.OUT ${results}/INFO_d${d}_revPBE.OUT

# LDA
sed -e 's/GGA_PBE_R/LDA_PW/g' -i input.xml
$exciting
cp INFO.OUT ${results}/INFO_d${d}_LDA.OUT

done

rm -rf draft.xml input.xml geometry.xml info.xml RHO3D.xml atoms.xml
rm -rf *.OUT

exit
