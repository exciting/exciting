
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_Gvector
!--------------------------------!
!     G-vector set variables     !
!--------------------------------!
! G-vector cut-off for interstitial potential and density
!replaced by inputstructurereal(8)::gmaxvr
! G-vector grid sizes
      Integer :: ngrid (3)
! total number of G-vectors
      Integer :: ngrtot
! integer grid intervals for each direction
      Integer :: intgv (3, 2)
! number of G-vectors with G < gmaxvr
      Integer :: ngvec
! G-vector integer coordinates
      Integer, Allocatable :: ivg (:, :)
! map from integer grid to G-vector array
      Integer, Allocatable :: ivgig (:, :, :)
! map from G-vector array to FFT array
      Integer, Allocatable :: igfft (:)
! G-vectors in Cartesian coordinates
      Real (8), Allocatable :: vgc (:, :)
! length of G-vectors
      Real (8), Allocatable :: gc (:)
! spherical harmonics of the G-vectors
      Complex (8), Allocatable :: ylmg (:, :)
! structure factor for the G-vectors
      Complex (8), Allocatable :: sfacg (:, :)
! G-space characteristic function: 0 inside the muffin-tins and 1 outside
      Complex (8), Allocatable :: cfunig (:)
! real-space characteristic function: 0 inside the muffin-tins and 1 outside
      Real (8), Allocatable :: cfunir (:)
! damping coefficient for characteristic function
!replaced by inputstructurereal(8)::cfdamp
End Module
!
