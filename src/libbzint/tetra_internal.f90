
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

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


      end module tetra_internal
!EOC      

! <sag>
module control
  implicit none
  ! interface parameter
  character(32), save :: tetraifc
  ! default to WIEN2k style (corresponding to original version)
  data tetraifc / 'wien2k' /
  ! level of debug output
  integer, save :: tetradbglv
  ! high default value to mimic original version
  data tetradbglv / 1000 /
  ! handling of pointers (problems with Portland compiler)
  integer, save :: pointerhandling
  ! default is 0 (original version); 1...explicit target assignment in
  ! routine "tetcw" (other routine to be followed)
  data pointerhandling / 0 /
  ! resonance type: 1...resonant term; 2...anti-resonant term; 0...both terms
  integer, save :: restype
  ! initialize with "0" according to original version
  data restype / 0 /
  ! switch wether k+q or k-q should be calculated
  logical, save :: kplusq
  ! initialize with ".false." according to original version
  data kplusq / .false. /
end module control

! set interface parameter
subroutine tetrasetifc(val)
  use control
  implicit none
  ! arguments
  character(*), intent(in) :: val
  ! local variables
  integer :: pass
  tetraifc=trim(adjustl(val))
  select case(trim(tetraifc))
     case('wien2k')
        pass=0
     case('exciting')
        pass=0
     case default
        pass=1
        write(*,*)
        write(*,'(a)') 'Error(libbzint): proper interface parameters are: &
             &"wien2k" and "exciting"'
        write(*,*)
        stop
  end select
end subroutine tetrasetifc

! set debug level
subroutine tetrasetdbglv(val)
  use control
  implicit none
  ! arguments
  integer, intent(in) :: val
  tetradbglv=val
end subroutine tetrasetdbglv

! set pointer handling
subroutine tetrasetpointerhandling(val)
  use control
  implicit none
  ! arguments
  integer, intent(in) :: val
  ! local variables
  integer :: pass
  select case(val)
     case(0)
        pass=0
     case(1)
        pass=0
     case(2)
        pass=0
     case default
        pass=1
        write(*,*)
        write(*,'(a)') 'Error(libbzint): proper pointer-handling parameters &
             &are: "0" and "1"'
        write(*,*)
        stop
  end select
  pointerhandling=val
end subroutine tetrasetpointerhandling

! set type of response function
subroutine tetrasetresptype(val)
  use control
  implicit none
  ! arguments
  integer, intent(in) :: val
  ! local variables
  integer :: pass
  select case(val)
     case(0)
        pass=0
     case(1)
        pass=0
     case(2)
        pass=0
     case default
        pass=1
        write(*,*)
        write(*,'(a)') 'Error(libbzint): proper response-type parameters are: &
             &"0", "1" and "2"'
        write(*,*)
        stop
  end select
  restype=val
end subroutine tetrasetresptype

! set treatment of q-shifted k-mesh
subroutine tetrasetkplusq(val)
  use control
  implicit none
  ! arguments
  logical, intent(in) :: val
  kplusq=val
end subroutine tetrasetkplusq

! report values of control module
subroutine tetrareportsettings
  use control
  implicit none
  write(*,'(a)')    'Info(libbzint): settings are:'
  write(*,'(a)')    '  interface        (default="wien2k") : '//trim(tetraifc)
  write(*,'(a,i6)') '  debug level      (default=1000)     :', tetradbglv
  write(*,'(a,i6)') '  pointer handling (default=0)        :', pointerhandling
  write(*,'(a,i6)') '  resonance type   (default=0)        :', restype
  write(*,'(a,l6)') '  kplusq           (default=F)        :', kplusq
  write(*,*)
end subroutine tetrareportsettings
!</sag>
