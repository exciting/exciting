!
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module modsym
      Implicit None
!
  !---------------------------!
  !     general variables     !
  !---------------------------!
  ! maximum allowed number of symmetry operations (private to this module)
      Integer, Private, Parameter :: maxsymcrs = 192
  ! true if symmetry group is abelian
      Logical :: abelsg
  ! true if symmetry group contains spatial inversion symmetry
      Logical :: spainvsym
  ! symmetry group multiplication table
      Integer, Allocatable :: sgmut (:, :)
!
  !--------------------!
  !     generators     !
  !--------------------!
  ! number of generators
      Integer :: ngenr
  ! number of elements in orbits of generators
      Integer, Allocatable :: negenr (:)
  ! generators
      Integer, Allocatable :: genr (:)
  ! orbits of generators
      Integer, Allocatable :: orbgenr (:, :)
!
  !-------------------!
  !     subgroups     !
  !-------------------!
  ! number of subgroups of space group
      Integer :: nsubsymc
  ! space group subgroups
      Integer :: subsymc (maxsymcrs, maxsymcrs)
!
  !---------------------------!
  !     conjugacy classes     !
  !---------------------------!
  ! number of classes of conjugated elements of spacegroup
      Integer :: nsymccocl
  ! classes of conjugated elements of spacegroup
      Integer :: symccocl (maxsymcrs, maxsymcrs)
  ! conjugacy relation between crystal symmetries
      Logical :: tsymcocl (maxsymcrs, maxsymcrs)
!
End Module
