


! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module mod_expiGqr

type point_type
 real(8) :: coord(3)
end type

type point_type_array
  type(point_type), pointer :: point
end type

type indexmap2d_type
  integer, pointer :: idx1of2(:), idx2of2(:)
end type

type expiGqr_chunk_type
  ! here is the actual DATA for the state-combinations (both state indices)
  complex(8), allocatable :: ssdata(:, :)
  ! here is the actual DATA for the state-combinations (combined state-state index)
  complex(8), allocatable :: cdata(:)
end type

type expiGqr_type
  ! inputs (general)
  type(point_type) :: qpoint
  type(point_type_array), pointer :: gqlist(:), klist(:)
  type(indexmap2d_type) :: stateIdx
  ! inputs (precision/algorithm)
  integer :: lmaxrayleigh, lmaxwave
  ! map other variables, handle via pointers

  ! DATA (temporary)
  real(8), allocatable :: rintaa(:, :, :, :, :, :, :), rintloa(:, :, :, :, :, :), rintlolo(:, :, :, :, :)
  complex(8), allocatable :: gntaa(:, :, :, :, :), gntalo(:, :, :, :, :), gntloa(:, :, :, :, :), gntlolo(:, :, :, :, :)
  ! DATA (status)
  logical :: radial_integrals_calculated
  ! DATA
  type(expiGqr_chunk_type), pointer :: gqarray(:)
end type

!expiGqr%gqarray(igq)%ssdata(:,:)
!expiGqr%gqarray(igq)%stateIdx%idx1of2(j)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------


!!! #include "something.F90"


function expiGqr_init()


  write(*, *) "Hello World!"
end function





end module
