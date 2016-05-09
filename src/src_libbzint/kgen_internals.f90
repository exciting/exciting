!BOP
!
! !MODULE: kgen_internals

      module kgen_internals

! !PUBLIC TYPES:      

      integer :: nirkp                  ! Number of irreducible k-points
      integer, dimension(3) :: div      ! Number of subdivisions of the BZ in each direction
      integer, allocatable  :: ikpid(:) ! Order number of the irreducible k-points
      integer, allocatable  :: redkp(:) ! Identification number of the irreducible k-point associated to the general
                                        ! k-point i.
                        
      integer, dimension(3) :: shift     ! Shift of the sublattice from the origin
      integer, allocatable :: iio(:,:,:) ! The symmetry operations matrices. Dimension >= nsymt
      real(8) :: vt                      ! Volume of the tetrahedra
      real(8), dimension(3,3) :: gbas    ! Basis vectors of the reciprocal lattice
      
      end module kgen_internals
!EOP      
