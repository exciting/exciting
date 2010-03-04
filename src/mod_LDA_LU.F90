
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_LDA_LU
!-------------------------!
!     LDA+U variables     !
!-------------------------!
! type of LDA+U to use (0: none)
      Integer :: ldapu
! maximum angular momentum
      Integer, Parameter :: lmaxlu = 3
      Integer, Parameter :: lmmaxlu = (lmaxlu+1) ** 2
! angular momentum for each species
      Integer :: llu (_MAXSPECIES_)
! U and J values for each species
      Real (8) :: ujlu (2, _MAXSPECIES_)
! LDA+U density matrix
      Complex (8), Allocatable :: dmatlu (:, :, :, :, :)
! LDA+U potential matrix in (l,m) basis
      Complex (8), Allocatable :: vmatlu (:, :, :, :, :)
! LDA+U energy for each atom
      Real (8), Allocatable :: engyalu (:)
! interpolation constant alpha for each atom (PRB 67, 153106 (2003))
      Real (8), Allocatable :: alphalu (:)
! energy from the LDA+U correction
      Real (8) :: engylu
End Module
!
