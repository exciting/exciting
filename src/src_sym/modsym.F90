
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modsym
  implicit none
  ! true if symmetry group is abelian
  logical :: abelsg
  ! symmetry group multiplication table
  integer, allocatable :: sgmut(:,:)
  ! number of generators
  integer :: ngenr
  ! number of elements in orbits of generators
  integer, allocatable :: negenr(:)
  ! generators
  integer, allocatable :: genr(:)
  ! orbits of generators
  integer, allocatable :: orbgenr(:,:)
end module
