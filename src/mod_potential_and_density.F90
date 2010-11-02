
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_potential_and_density
      Use modinput
!-----------------------------------------!
!     potential and density variables     !
!-----------------------------------------!
! exchange-correlation functional type
 	integer::xctype(3)
! exchange-correlation functional description
      Character (512) :: xcdescr
! exchange-correlation functional spin treatment
      Integer :: xcspin
! exchange-correlation functional density gradient treatment
      Integer :: xcgrad
!
      Integer :: maxncv
! muffin-tin charge density
      Real (8), Allocatable :: rhomt (:, :, :)
! interstitial real-space charge density
      Real (8), Allocatable :: rhoir (:)
! muffin-tin charge density (reference)
      Real (8), Allocatable :: rhomtref (:, :, :)
! interstitial real-space charge density (reference)
      Real (8), Allocatable :: rhoirref (:)
! muffin-tin magnetisation vector field
      Real (8), Allocatable :: magmt (:, :, :, :)
! interstitial magnetisation vector field
      Real (8), Allocatable :: magir (:, :)
! muffin-tin Coulomb potential
      Real (8), Allocatable :: vclmt (:, :, :)
! interstitial real-space Coulomb potential
      Real (8), Allocatable :: vclir (:)
! order of polynomial for pseudocharge density
!replaced by inputstructureinteger::npsden
! muffin-tin exchange-correlation potential
      Real (8), Allocatable :: vxcmt (:, :, :)
! interstitial real-space exchange-correlation potential
      Real (8), Allocatable :: vxcir (:)
! muffin-tin exchange-correlation magnetic field
      Real (8), Allocatable :: bxcmt (:, :, :, :)
! interstitial exchange-correlation magnetic field
      Real (8), Allocatable :: bxcir (:, :)
! nosource is .true. if the field is to be made source-free
!replaced by inputstructurelogical::nosource
! muffin-tin effective potential
      Real (8), Allocatable :: veffmt (:, :, :)
! interstitial effective potential
      Real (8), Allocatable :: veffir (:)
! G-space interstitial effective potential
      Complex (8), Allocatable :: veffig (:)
! muffin-tin exchange energy density
      Real (8), Allocatable :: exmt (:, :, :)
! interstitial real-space exchange energy density
      Real (8), Allocatable :: exir (:)
! muffin-tin correlation energy density
      Real (8), Allocatable :: ecmt (:, :, :)
! interstitial real-space correlation energy density
      Real (8), Allocatable :: ecir (:)
! type of mixing to use for the potential
!replaced by inputstructureinteger::mixtype
! adaptive mixing parameters
!replaced by inputstructurereal(8)::beta0
!replaced by inputstructurereal(8)::betainc
!replaced by inputstructurereal(8)::betadec
End Module
!
