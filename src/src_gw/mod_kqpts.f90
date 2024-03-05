
!------------------------------------------
! Original definitions of k-/q-point sets
! Supposed to be completely replaced by 
! mod_kpointset.f90
!------------------------------------------

module mod_kqpts
    use asserts, only: assert
    ! We need to change internally to `terminate_when_false` to avoid a circular dependency
    use modmpi, only: terminate_when_false => terminate_if_false
    use precision, only: i32

    implicit none
    
    !-------------------------------!
    ! tetrahedron method variables  !
    !-------------------------------!
    integer(4), allocatable :: idikp(:)
    integer(4), allocatable :: kqid(:,:)
    integer(4) :: dvq

    integer(4) :: ntetnr                    ! Total number of tetrahedra
    integer(4), allocatable :: wtetnr(:)    ! weight of each tetrahedron  for integration
    integer(4), allocatable :: tnodesnr(:,:)! index of the k-points corresponding to the nodes of each tetrahedra for integration

    integer(4), allocatable :: linkq(:,:)

    !--------------------------!
    !  Non-reduced G+k arrays  !
    !--------------------------!
    ! number of G+k-vectors for augmented plane waves
    Integer, Allocatable :: ngknr(:,:)
    ! index from G+k-vectors to G-vectors
    Integer, Allocatable :: igkignr(:,:,:)
    ! G+k-vectors in lattice coordinates
    Real (8), Allocatable :: vgklnr(:,:,:,:)
    ! G+k-vectors in Cartesian coordinates
    Real (8), Allocatable :: vgkcnr(:,:,:,:)
    ! length of G+k-vectors
    Real (8), Allocatable :: gkcnr(:,:,:)
    ! (theta, phi) coordinates of G+k-vectors
    Real (8), Allocatable :: tpgkcnr(:,:,:,:)
    ! structure factor for the G+k-vectors
    Complex (8), Allocatable :: sfacgknr(:,:,:,:)

    !--------------------------------!
    !     Small group of q-vectors   !
    !--------------------------------! 
    ! non-reduced number of q-points
    integer :: nqptnr
    ! number of the symmetry operations in the small group of q
    integer, allocatable :: nsymq(:) 
    ! q-dependent k-point weight
    real(8), allocatable :: wkpq(:,:)
    ! number of k-points in IBZ(q)
    integer, allocatable :: nkptq(:)
    ! index of the symmetry operation which rotates the k-point into equivalent one
    integer, allocatable :: iksymq(:,:)
    ! map the k-point index to the corresponding irreducible one
    integer, allocatable :: indkpq(:,:)
    ! map the irreducible k-point index to the corresponding from the non-reduced set
    integer, allocatable :: idikpq(:,:)
    ! rotation matrix for ylm's      
    complex(8), allocatable :: djmm(:,:)
      
    integer, allocatable :: nsymkstar(:,:), isymkstar(:,:,:)
    
    
    ! number of G-vectors for the bare coulomb matrix
    integer, allocatable :: ngbarc(:)
      
    ! map from G+q-vectors to G-vectors for the (increased) coulomb Gmax cutoff
    integer, allocatable :: igqigb(:,:)
      
    ! map from G-vectors to G+q-vectors for the (increased) coulomb Gmax cutoff
    integer, allocatable :: igigqb(:,:)
    
    
    
    
    ! reduced set of eigenvectors of barcoul matrix after barcevtol
    complex(8), allocatable :: vbas(:,:)
      
    ! transform matrix that diagonalized original bare Coulomb matrix 
    complex(8), allocatable :: barcvm(:,:)  
    
    !> Interface to be used for the ranges of k/q-point indexes defined in the input file
    type, private :: ranges_of_indexes
      integer(i32) :: first
      integer(i32) :: last
    contains
      procedure :: parse_first_last
      procedure :: sanity_checks
    end type 

contains 

    subroutine sanity_checks( this, first, last, maximum )
      class(ranges_of_indexes), intent(in) :: this
      integer, intent(in) :: first
      integer, intent(in) :: last
      integer, intent(in) :: maximum
    
      call terminate_when_false( first>0, 'first k/q-point must be positive' )
      ! Case qf <= 0 is interpreted as qf = n_qpt
      if( last <=0 ) then
        call terminate_when_false( first<=maximum, 'first k/q-point must <= the number of k/q-points' )
      else
        call terminate_when_false( first<=last, 'first k/q-point must be <= last k/q-point' )
        call terminate_when_false( last<=maximum, 'last k/q-point must be <= the number of k/q-points' )
      end if

    end subroutine
    

    subroutine parse_first_last( this, first, last, maximum )
      class(ranges_of_indexes), intent(inout)  :: this
      integer, intent(in)   :: first, last, maximum
    
      call this%sanity_checks( first, last, maximum )
      this%first = first
      this%last = merge( tsource=last, fsource=maximum, mask=last>=0 ) 
    
    end
    
    
end module
