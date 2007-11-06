!BOP
!
! !MODULE: tetra_internal
!
! !INTERFACE:
      module tetra_internal
       
      integer(4) :: nirkp                 ! Number of irreducible k-points
      
      integer(4) :: ntet                  ! Total number of tetrahedra
      
      integer(4) :: ncore                 ! Maximum number of core states

      integer(4) :: nband                 ! Maximum number of bands.
       
      integer(4), pointer :: tetcorn(:,:) ! Id. numbers of the corner of
!                                           the tetrahedra
 
      integer(4), pointer :: tetweig(:)   ! Weight of each tetrahedron
      
      integer(4), pointer :: qweig(:,:)    ! Weight of each q-point
 
      integer(4), pointer :: tetln(:)     ! linked tetrahedron for q.
       
      real(8) :: vt                       ! Relative volume of the
!                                           tetrahedra      
      
      real(8)   , pointer :: eband(:,:)   ! Band energies

      real(8)   , pointer :: ecore(:)     ! core energies
       
      real(8)  :: omgga                  ! the frequency to be included

      integer(4)  :: sgnfrq              ! a sign to tell which weight to 
!                                          be calculated

!<sag>
      integer :: restype                 ! if "1" the resonant frequency term 
                                         !   is calculated
                                         ! if "2" the anti-resonant frequency
                                         !   term is calculated
                                         ! if "0" both parts are calculated
                                         !   together
      data restype / 0 /                 ! initialize with "0" according to
                                         !   original version
!</sag>

      end module tetra_internal
!EOC      

!<sag>
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module control
  implicit none
  ! level of debug output
  integer, save :: tetradbglv
  ! high default value to mimic original version that outputs a lot of
  ! information
  data tetradbglv / 1000 /
  ! handling of pointers
  integer, save :: pointerhandling
  ! default is 0 (original version); 1...explicit target assignment in
  ! routine "tetcw"
  data pointerhandling / 0 /

end module control

! set debug level
subroutine tetrasetdbglv(lev)
  use control
  implicit none
  integer, intent(in) :: lev
  tetradbglv=lev
end subroutine tetrasetdbglv

! set pointer handling
subroutine setpointerhandling(val)
  use control
  implicit none
  integer, intent(in) :: val
  pointerhandling=val
end subroutine setpointerhandling

!</sag>
