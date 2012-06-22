
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !MODULE: kgen_internals

      module kgen_internals

! !PUBLIC TYPES:      

      integer(4) :: nirkp            ! Number of irreducible k-points

      integer(4), dimension(3) :: div      ! Number of subdivisions of 
!                                            the BZ in each direction

      integer(4), allocatable  :: ikpid(:) ! Order number of the 
!                                            irreducible k-points
                                            
      integer(4), allocatable  :: redkp(:) ! Identification number of 
!                                            the irreducible k-point
!                                            associated to the general
!                                            k-point i.
      integer(4), allocatable :: redtet(:)
                                            
      integer(4) :: mndg                                      
                                            
      integer(4), dimension(3) :: shift    ! Shift of the sublattice from
!                                            the origin
                                            
      integer(4), allocatable :: iio(:,:,:) ! The symmetry operations
!                                            matrices. 
!                                            Dimension >= nsymt

      real(8) :: vt                    ! Volume of the tetrahedra
      
      real(8), dimension(3,3) :: gbas  ! Basis vectors of the
!                                        reciprocal lattice
           
      
      end module kgen_internals
!EOP      
