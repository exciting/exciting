
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_OEP_HF
!----------------------------------------!
!     OEP and Hartree-Fock variables     !
!----------------------------------------!
! maximum number of core states over all species
      Integer :: ncrmax
! maximum number of OEP iterations
!replaced by inputstructureinteger::maxitoep
! initial value and scaling factors for OEP step size
!replaced by inputstructurereal(8)::tauoep(3)
! magnitude of the OEP residual
      Real (8) :: resoep
! kinetic matrix elements
      Complex (8), Allocatable :: kinmatc (:, :, :)
! complex versions of the exchange potential and field
      Complex (8), Allocatable :: zvxmt (:, :, :)
      Complex (8), Allocatable :: zvxir (:)
      Complex (8), Allocatable :: zbxmt (:, :, :, :)
      Complex (8), Allocatable :: zbxir (:, :)
End Module
!
