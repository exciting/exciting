!BOP
!
! !MODULE: tetra_internal
!
! !INTERFACE:
      module tetra_internal

      integer :: fout = 10
      logical :: ldbg_bzint = .false.
      integer :: iop_integ = 0         ! the option to do the tetrahedron integration 
                                       !  0 -- analytical integration 
                                       !  1 -- numerical integration 
      integer :: n_gauq = 10 
      real(8), allocatable:: x_gauq(:),w_gauq(:) 
      real(8) :: eta_refreq=0.005      ! this is the value of the small imaginary part that 
                                       ! is needed for real frequency calculations 

      real(8),parameter:: weighttol = 1.0e+4
      real(8),parameter:: weightwarn = 1.0e+3
      real(8):: ztol_vol = 1.e-10
      real(8):: tol_taylor = 10.0
      real(8):: ztol_sorteq = 1.0*1.e-2

      integer :: nirkp 
      integer :: ntet                  ! Total number of tetrahedra
      integer :: ncore                 ! Maximum number of core states
      integer :: nband                 ! Maximum number of bands.
      integer :: mndg              
      integer,allocatable:: redtet(:)
      integer,pointer:: tetcorn(:,:) ! Id. numbers of the corner of the tetrahedra
      integer,pointer:: tetweig(:)   ! Weight of each tetrahedron
      integer,pointer:: qweig(:,:)    ! Weight of each q-point
      integer,pointer:: tetln(:)     ! linked tetrahedron for q.
      integer,pointer:: klinkq(:)    ! linked kpoint for q.
      real(8) :: vt                       ! Relative volume of the tetrahedra      
      real(8),pointer :: eband(:,:)   ! Band energies
      real(8),pointer :: ecore(:)     ! core energies
      real(8):: omgga                  ! the frequency to be included
      integer:: sgnfrq              ! a sign to tell which weight to be calculated
      real(8):: vol_small_tetra
      end module tetra_internal
!EOC      
