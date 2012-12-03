#!/bin/tcsh 
#
# AUTHOR
# juergen.spitaler@unileoben.ac.at
#
# DATE: Mon Aug 30 15:32:36 CEST 2010
#
# SYNTAX:
# E_formation.csh scf_file(s)
#
# EXPLANATION:
# Calculates the formation energy of a NiTiHf structure,
# using as ref. energies 
#   E0 = energy per Al atom in pure fcc Al, in eV
#   E1 = energy per Ti atom in pure fcc Ti, in eV
#

# Check if help is needed

if (  "$1" == "-h" || "$1" == "-H" ) then
    goto help
endif

foreach i ($argv)
  set n_Al = `awk 'BEGIN{nat=0};/species chemicalSymbol="Al"/{getline; while($0 !~ /species>/) {getline; nat++}}END{print nat}' $i` 
  set n_Ti = `awk 'BEGIN{nat=0}/species chemicalSymbol="Ti"/{getline;while($0 !~ /species>/) {getline; nat++}}END{print nat}' $i` 
    echo "n_Al, n_Ti =", $n_Al, $n_Ti
   ex_lene $i | awk -v fname=$i -v n_Al=$n_Al -v n_Ti=$n_Ti '\
    BEGIN{Ha=2*13.6058; E0= -242.82363185*Ha;   E1=-853.80938594*Ha;};\
      {ene=$1; fname = $3};\
    END{printf "n_Ti = %2d   n_Al = %2d   Ene = %15.6f eV   E_formation = %12.6f  ... %s\n", n_Ti, n_Al, ene*Ha, (ene*Ha - n_Al*E0 - n_Ti*E1)/(n_Ti+n_Al), fname}' 
end    
exit

help:
awk '/# *SYNTAX:/,$0 !~ / *#.*[a-zA-Z]+/;\
     /# *EXPLANATION:/,$0 !~ / *#.*[a-zA-Z]+/' < $0
