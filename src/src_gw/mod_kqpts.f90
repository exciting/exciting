
!------------------------------------------
! Original definitions of k-/q-point sets
! Supposed to be completely replaced by 
! mod_kpointset.f90
!------------------------------------------

module mod_kqpts
#include "asserts.fpp"
    use modinput, only: qpoints_type_array, kpoints_type_array
    ! We need to change internally to `terminate_when_false` to avoid a circular dependency
    use modmpi, only: terminate_when_false => terminate_if_false
    use precision, only: i32, dp

    implicit none
    
    !-------------------------------!
    ! tetrahedron method variables  !
    !-------------------------------!
    integer(i32), allocatable :: idikp(:)
    integer(i32), allocatable :: kqid(:,:)
    integer(i32) :: dvq

    integer(i32) :: ntetnr                    ! Total number of tetrahedra
    integer(i32), allocatable :: wtetnr(:)    ! weight of each tetrahedron  for integration
    integer(i32), allocatable :: tnodesnr(:,:)! index of the k-points corresponding to the nodes of each tetrahedra for integration

    integer(i32), allocatable :: linkq(:,:)

    !--------------------------!
    !  Non-reduced G+k arrays  !
    !--------------------------!
    ! number of G+k-vectors for augmented plane waves
    integer(i32), Allocatable :: ngknr(:,:)
    ! index from G+k-vectors to G-vectors
    integer(i32), Allocatable :: igkignr(:,:,:)
    ! G+k-vectors in lattice coordinates
    Real (dp), Allocatable :: vgklnr(:,:,:,:)
    ! G+k-vectors in Cartesian coordinates
    Real (dp), Allocatable :: vgkcnr(:,:,:,:)
    ! length of G+k-vectors
    Real (dp), Allocatable :: gkcnr(:,:,:)
    ! (theta, phi) coordinates of G+k-vectors
    Real (dp), Allocatable :: tpgkcnr(:,:,:,:)
    ! structure factor for the G+k-vectors
    Complex (dp), Allocatable :: sfacgknr(:,:,:,:)

    !--------------------------------!
    !     Small group of q-vectors   !
    !--------------------------------! 
    ! non-reduced number of q-points
    integer(i32) :: nqptnr
    ! number of the symmetry operations in the small group of q
    integer(i32), allocatable :: nsymq(:) 
    ! q-dependent k-point weight
    real(dp), allocatable :: wkpq(:,:)
    ! number of k-points in IBZ(q)
    integer(i32), allocatable :: nkptq(:)
    ! index of the symmetry operation which rotates the k-point into equivalent one
    integer(i32), allocatable :: iksymq(:,:)
    ! map the k-point index to the corresponding irreducible one
    integer(i32), allocatable :: indkpq(:,:)
    ! map the irreducible k-point index to the corresponding from the non-reduced set
    integer(i32), allocatable :: idikpq(:,:)
    ! rotation matrix for ylm's      
    complex(dp), allocatable :: djmm(:,:)
      
    integer(i32), allocatable :: nsymkstar(:,:), isymkstar(:,:,:)
    
    
    ! number of G-vectors for the bare coulomb matrix
    integer(i32), allocatable :: ngbarc(:)
      
    ! map from G+q-vectors to G-vectors for the (increased) coulomb Gmax cutoff
    integer(i32), allocatable :: igqigb(:,:)
      
    ! map from G-vectors to G+q-vectors for the (increased) coulomb Gmax cutoff
    integer(i32), allocatable :: igigqb(:,:)
    
    
    
    
    ! reduced set of eigenvectors of barcoul matrix after barcevtol
    complex(dp), allocatable :: vbas(:,:)
      
    ! transform matrix that diagonalized original bare Coulomb matrix 
    complex(dp), allocatable :: barcvm(:,:)  
    
    !> Interface to be used for the ranges of k/q-point indexes defined in the input file
    type, private :: ranges_of_indexes
      integer(i32) :: first
      integer(i32) :: last
    contains
      procedure :: parse_first_last
      procedure :: sanity_checks
    end type 

    !> Interface to the list of ranges of k/q-point indexes given in the input file
    type, public :: kpoints_sets
      type(ranges_of_indexes), allocatable :: sets(:)
      integer(i32), allocatable :: list_of_indexes(:)
    contains 
      generic, public :: parse_input => parse_input_qpoints, parse_input_kpoints
      procedure, private :: parse_input_qpoints
      procedure, private :: parse_input_kpoints
      procedure :: obtain_list_of_indexes
    end type

    public :: has_full_k_point_coverage

contains 

!> Return whether the selected k-/q-point indexes cover the full set exactly once.
pure logical function has_full_k_point_coverage(kpoint_indexes, n_kpoints) result(flag)
  !> Selected k-/q-point indexes.
  integer(i32), intent(in) :: kpoint_indexes(:)
  !> Total number of k-/q-points in the full set.
  integer(i32), intent(in) :: n_kpoints

  integer(i32) :: ik

  flag = size(kpoint_indexes) == n_kpoints .and. &
    all( [(any(kpoint_indexes == ik), ik = 1, n_kpoints)] )
end function has_full_k_point_coverage

!> Check that an input k-/q-point range is valid.
subroutine sanity_checks( this, first, last, maximum )
  !> Parsed range object.
  class(ranges_of_indexes), intent(in) :: this
  !> First requested k-/q-point index.
  integer, intent(in) :: first
  !> Last requested k-/q-point index, or non-positive to select through `maximum`.
  integer, intent(in) :: last
  !> Maximum available k-/q-point index.
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


!> Store a validated first/last k-/q-point range.
subroutine parse_first_last( this, first, last, maximum )
  !> Range object to update.
  class(ranges_of_indexes), intent(inout)  :: this
  !> First requested k-/q-point index.
  integer, intent(in) :: first
  !> Last requested k-/q-point index, or non-positive to select through `maximum`.
  integer, intent(in) :: last
  !> Maximum available k-/q-point index.
  integer, intent(in) :: maximum

  call this%sanity_checks( first, last, maximum )
  this%first = first
  this%last = merge( tsource=last, fsource=maximum, mask=last>=0 )

end subroutine


!> Parse q-point ranges from XML input.
subroutine parse_input_qpoints( this, qpoints_array, n_qpoints_max )
  !> q-point set object to update.
  class(kpoints_sets), intent(inout) :: this
  !> Parsed XML q-point ranges.
  type(qpoints_type_array), pointer, intent(in) :: qpoints_array(:)
  !> Maximum available q-point index.
  integer(i32), intent(in) :: n_qpoints_max

  integer(i32) :: n, i

  n = size( qpoints_array )
  allocate( this%sets(n) )
  do i = 1, n
    call this%sets(i)%parse_first_last( qpoints_array(i)%qpoints%qi, qpoints_array(i)%qpoints%qf, n_qpoints_max )
  end do

end subroutine


!> Parse k-point ranges from XML input.
subroutine parse_input_kpoints( this, kpoints_array, n_kpoints_max )
  !> k-point set object to update.
  class(kpoints_sets), intent(inout) :: this
  !> Parsed XML k-point ranges.
  type(kpoints_type_array), pointer, intent(in) :: kpoints_array(:)
  !> Maximum available k-point index.
  integer(i32), intent(in) :: n_kpoints_max

  integer(i32) :: n, i

  n = size( kpoints_array )
  allocate( this%sets(n) )
  do i = 1, n
    call this%sets(i)%parse_first_last( kpoints_array(i)%kpoints%ki, kpoints_array(i)%kpoints%kf, n_kpoints_max )
  end do

end subroutine


!> Build the explicit list of selected k-/q-point indexes.
subroutine obtain_list_of_indexes( this )
  !> k-/q-point set object to update.
  class(kpoints_sets), intent(inout) :: this

  integer(i32) :: i, j

  CALL_ASSERT( size( this%sets ) >= 1, 'sets must contain at least one element' )
  this%list_of_indexes = [ (i, i=this%sets(1)%first,this%sets(1)%last) ]
  do j = 2, size( this%sets )
    this%list_of_indexes = [ this%list_of_indexes, ( i, i=this%sets(j)%first, this%sets(j)%last ) ]
  end do
end subroutine


end module
