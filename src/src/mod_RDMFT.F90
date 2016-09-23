
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_RDMFT
!-------------------------------------------------------------!
!     reduced density matrix functional (RDMFT) variables     !
!-------------------------------------------------------------!
! non-local matrix elements for varying occupation numbers
      Real (8), Allocatable :: vnlrdm (:, :, :, :)
! Coulomb potential matrix elements
      Complex (8), Allocatable :: vclmat (:, :, :)
! derivative of kinetic energy w.r.t. natural orbital coefficients
      Complex (8), Allocatable :: dkdc (:, :, :)
! step size for occupation numbers
!replaced by inputstructurereal(8)::taurdmn
! step size for natural orbital coefficients
!replaced by inputstructurereal(8)::taurdmc
! xc functional
!replaced by inputstructureinteger::rdmxctype
! maximum number of self-consistent loops
!replaced by inputstructureinteger::rdmmaxscl
! maximum number of iterations for occupation number optimisation
!replaced by inputstructureinteger::maxitn
! maximum number of iteration for natural orbital optimisation
!replaced by inputstructureinteger::maxitc
! exponent for the functional
!replaced by inputstructurereal(8)::rdmalpha
! temperature
!replaced by inputstructurereal(8)::rdmtemp
! entropy
      Real (8) :: rdmentrpy
End Module
!
