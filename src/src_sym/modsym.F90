


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modsym
  implicit none

  !---------------------------!
  !     general variables     !
  !---------------------------!
  ! maximum allowed number of symmetry operations (private to this module)
  integer, private, parameter :: maxsymcrs=192
  ! true if symmetry group is abelian
  logical :: abelsg
  ! symmetry group multiplication table
  integer, allocatable :: sgmut(:, :)

  !--------------------!
  !     generators     !
  !--------------------!
  ! number of generators
  integer :: ngenr
  ! number of elements in orbits of generators
  integer, allocatable :: negenr(:)
  ! generators
  integer, allocatable :: genr(:)
  ! orbits of generators
  integer, allocatable :: orbgenr(:, :)

  !-------------------!
  !     subgroups     !
  !-------------------!
  ! number of subgroups of space group
  integer :: nsubsymc
  ! space group subgroups
  integer :: subsymc(maxsymcrs, maxsymcrs)

  !---------------------------!
  !     conjugacy classes     !
  !---------------------------!
  ! number of classes of conjugated elements of spacegroup
  integer :: nsymccocl
  ! classes of conjugated elements of spacegroup
  integer :: symccocl(maxsymcrs, maxsymcrs)
  ! conjugacy relation between crystal symmetries
  logical :: tsymcocl(maxsymcrs, maxsymcrs)

end module
