
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_APW_LO
!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! maximum allowable APW order
      Integer, Parameter :: maxapword = 4
! APW order
      Integer :: apword (0:_MAXLAPW_, _MAXSPECIES_)
! maximum of apword over all angular momenta and species
      Integer :: apwordmax
! APW initial linearisation energies
      Real (8) :: apwe0 (maxapword, 0:_MAXLAPW_, _MAXSPECIES_)
! APW linearisation energies
      Real (8), Allocatable :: apwe (:, :, :)
! APW derivative order
      Integer :: apwdm (maxapword, 0:_MAXLAPW_, _MAXSPECIES_)
! apwve is .true. if the linearisation energies are allowed to vary
      Logical :: apwve (maxapword, 0:_MAXLAPW_, _MAXSPECIES_)
! APW radial functions
      Real (8), Allocatable :: apwfr (:, :, :, :, :)
! derivate of radial functions at the muffin-tin surface
      Real (8), Allocatable :: apwdfr (:, :, :)
! maximum number of local-orbitals
      Integer, Parameter :: maxlorb = 100
! maximum allowable local-orbital order
      Integer, Parameter :: maxlorbord = 4
! number of local-orbitals
      Integer :: nlorb (_MAXSPECIES_)
! maximum nlorb over all species
      Integer :: nlomax
! total number of local-orbitals
      Integer :: nlotot
! local-orbital order
      Integer :: lorbord (maxlorb, _MAXSPECIES_)
! local-orbital angular momentum
      Integer :: lorbl (maxlorb, _MAXSPECIES_)
! maximum lorbl over all species
      Integer :: lolmax
! (lolmax+1)^2
      Integer :: lolmmax
! local-orbital initial energies
      Real (8) :: lorbe0 (maxlorbord, maxlorb, _MAXSPECIES_)
! local-orbital energies
      Real (8), Allocatable :: lorbe (:, :, :)
! local-orbital derivative order
      Integer :: lorbdm (maxlorbord, maxlorb, _MAXSPECIES_)
! lorbve is .true. if the linearisation energies are allowed to vary
      Logical :: lorbve (maxlorbord, maxlorb, _MAXSPECIES_)
! local-orbital radial functions
      Real (8), Allocatable :: lofr (:, :, :, :)
! energy step size for locating the band energy
!replaced by inputstructurereal(8)::deband
! minimum of the default linearisation energy over all APW and local-orbitals
! functions
      real(8) :: mine0
End Module
!
