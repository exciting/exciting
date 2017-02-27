
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_Gkvector
!----------------------------------!
!     G+k-vector set variables     !
!----------------------------------!
! smallest muffin-tin radius times gkmax
!replaced by inputstructurereal(8)::rgkmax
! maximum |G+k| cut-off for APW functions
      Real (8) :: gkmax
! number of G+k-vectors for augmented plane waves
      integer, pointer :: ngk_ptr(:,:)
      Integer, Allocatable, target :: ngk (:, :)
! maximum number of G+k-vectors over all k-points
      integer, pointer :: ngkmax_ptr
      Integer, target :: ngkmax
! index from G+k-vectors to G-vectors
      Integer, Allocatable, target :: igkig (:, :, :)
! G+k-vectors in lattice coordinates
      Real (8), pointer :: vgkl_ptr(:, :, :, :)
      Real (8), Allocatable, target :: vgkl (:, :, :, :)
! G+k-vectors in Cartesian coordinates
      Real (8), Allocatable, target :: vgkc (:, :, :, :)
! length of G+k-vectors
      Real (8), Allocatable, target :: gkc (:, :, :)
! (theta, phi) coordinates of G+k-vectors
      Real (8), Allocatable, target :: tpgkc (:, :, :, :)
! structure factor for the G+k-vectors
      Complex (8), Allocatable, target :: sfacgk (:, :, :, :)
! dimensions of the FFT grid for APW functions
      Integer :: ngkfft(3)
! number of FFT grid points
      Integer :: ngktotfft 
! mapping to the FFT grid
      Integer, Allocatable :: igkfft (:,:)
      
End Module
!
