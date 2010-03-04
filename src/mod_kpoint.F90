
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_kpoint
!-------------------------------!
!     k-point set variables     !
!-------------------------------!
! autokpt is .true. if the k-point set is determined automatically
!replaced by inputstructurelogical::autokpt
! radius of sphere used to determine k-point density when autokpt is .true.
!replaced by inputstructurereal(8)::radkpt
! k-point grid sizes
!replaced by inputstructureinteger::ngridk(3)
! total number of k-points
      Integer :: nkpt
! k-point offset
!replaced by inputstructurereal(8)::vkloff(3)
! reducek is .true. if k-points are to be reduced (with crystal symmetries)
!replaced by inputstructurelogical::reducek
! locations of k-points on integer grid
      Integer, Allocatable :: ivk (:, :)
! k-points in lattice coordinates
      Real (8), Allocatable :: vkl (:, :)
! k-points in Cartesian coordinates
      Real (8), Allocatable :: vkc (:, :)
! k-point weights
      Real (8), Allocatable :: wkpt (:)
! map from non-reduced grid to reduced set
      Integer, Allocatable :: ikmap (:, :, :)
! total number of non-reduced k-points
      Integer :: nkptnr
! locations of non-reduced k-points on integer grid
      Integer, Allocatable :: ivknr (:, :)
! non-reduced k-points in lattice coordinates
      Real (8), Allocatable :: vklnr (:, :)
! non-reduced k-points in Cartesian coordinates
      Real (8), Allocatable :: vkcnr (:, :)
! non-reduced k-point weights
      Real (8), Allocatable :: wkptnr (:)
! map from non-reduced grid to non-reduced set
      Integer, Allocatable :: ikmapnr (:, :, :)
! k-point at which to determine effective mass tensor
!replaced by inputstructurereal(8)::vklem(3)
! displacement size for computing the effective mass tensor
!replaced by inputstructurereal(8)::deltaem
! number of displacements in each direction
!replaced by inputstructureinteger::ndspem
End Module
!
