
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_charge_and_moment
!-------------------------------------!
!     charge and moment variables     !
!-------------------------------------!
! tolerance for error in total charge
!replaced by inputstructurereal(8)::epschg
! total nuclear charge
      Real (8) :: chgzn
! total core charge
      Real (8) :: chgcr
! core leakage charge
      Real (8) :: chgcrlk
! total valence charge
      Real (8) :: chgval
! excess charge
!replaced by inputstructurereal(8)::chgexs
! total charge
      Real (8) :: chgtot
! calculated total charge
      Real (8) :: chgcalc
! interstitial region charge
      Real (8) :: chgir
! muffin-tin charges
      Real (8), Allocatable :: chgmt (:)
! total muffin-tin charge
      Real (8) :: chgmttot
! partial charges
      real(8), allocatable :: chgpart(:,:,:)
! charge distance
      real(8) :: chgdst
! effective Wigner radius
      Real (8) :: rwigner
! total moment
      Real (8) :: momtot (3)
! interstitial region moment
      Real (8) :: momir (3)
! muffin-tin moments
      Real (8), Allocatable :: mommt (:, :)
! total muffin-tin moment
      Real (8) :: mommttot (3)
End Module
!
