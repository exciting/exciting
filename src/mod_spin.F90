
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_spin
!--------------------------------!
!     spin related variables     !
!--------------------------------!
! spinpol is .true. for spin-polarised calculations
!replaced by inputstructurelogical::spinpol
! spinorb is .true. for spin-orbit coupling
!replaced by inputstructurelogical::spinorb
! fixspin type: 0 = none, 1 = global, 2 = local, 3 = global + local
!replaced by inputstructureinteger::fixspin
! dimension of magnetisation and magnetic vector fields (1 or 3)
      Integer :: ndmag
! ncmag is .true. if the magnetisation is non-collinear, i.e. when ndmag = 3
      Logical :: ncmag
! fixed total spin magnetic moment
!replaced by inputstructurereal(8)::momfix(3)
! fixed spin moment global effective field in Cartesian coordinates
      Real (8) :: bfsmc (3)
! muffin-tin fixed spin moments
!replaced by inputstructurereal(8)::mommtfix(3, _MAXATOMS_, _MAXSPECIES_)
! muffin-tin fixed spin moment effective fields in Cartesian coordinates
      Real (8) :: bfsmcmt (3, _MAXATOMS_, _MAXSPECIES_)
! fixed spin moment field mixing parameter
!replaced by inputstructurereal(8)::taufsm
! second-variational spinor dimension (1 or 2)
      Integer :: nspinor
! external magnetic field in each muffin-tin in lattice coordinates
      Real (8) :: bflmt (3, _MAXATOMS_, _MAXSPECIES_)
! external magnetic field in each muffin-tin in Cartesian coordinates
!replaced by inputstructurereal(8)::bfcmt(3, _MAXATOMS_, _MAXSPECIES_)
! global external magnetic field in lattice coordinates
      Real (8) :: bfieldl (3)
! global external magnetic field in Cartesian coordinates
!replaced by inputstructurereal(8)::bfieldc(3)
! external magnetic fields are multiplied by reducebf after each iteration
!replaced by inputstructurereal(8)::reducebf
! spinsprl if .true. if a spin-spiral is to be calculated
!replaced by inputstructurelogical::spinsprl
! number of spin-dependent first-variational functions per state
      Integer :: nspnfv
! spin-spiral q-vector in lattice coordinates
!replaced by inputstructurereal(8)::vqlss(3)
! spin-spiral q-vector in Cartesian coordinates
      Real (8) :: vqcss (3)
End Module
!
