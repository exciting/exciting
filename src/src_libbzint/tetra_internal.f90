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
