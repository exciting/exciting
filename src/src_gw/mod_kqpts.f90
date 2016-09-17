
!------------------------------------------
! Original definitions of k-/q-point sets
! Supposed to be completely replaced by 
! mod_kpointset.f90
!------------------------------------------

module mod_kqpts

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

end module
