
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_qpoint
!-------------------------------!
!     q-point set variables     !
!-------------------------------!
! q-point grid sizes
      Integer :: ngridq (3)
! total number of q-points
      Integer :: nqpt
! reduceq is .true. if q-points are to be reduced (with crystal symmetries)
!replaced by inputstructurelogical::reduceq
! locations of q-points on integer grid
      Integer, Allocatable :: ivq (:, :)
! map from non-reduced grid to reduced set
      Integer, Allocatable :: iqmap (:, :, :)
! q-points in lattice coordinates
      Real (8), Allocatable :: vql (:, :)
! q-points in Cartesian coordinates
      Real (8), Allocatable :: vqc (:, :)
! q-point weights
      Real (8), Allocatable :: wqpt (:)
! weights associated with the integral of 1/q^2
      Real (8), Allocatable :: wiq2 (:)
End Module
!
