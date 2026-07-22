!> Module for constrained DFT (CDFT) calculations
!> Created Oct 2024 (Ronaldo)
module cdft
#include "asserts.fpp"
  use constants, only: real_one, real_zero, zzero
  use math_utils, only: all_close
  use mod_mpi_env, only: mpiinfo
  use modinput, only: groundstate_type
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
  use xlapack, only: hermitian_matrix_multiply, matrix_multiply

  implicit none

  private

  public :: deallocate_cdft_global_arrays, &
            determine_cdft_occupations, &
            initialize_cdft_global_arrays, &
            set_overlap_times_psi_gs, &
            set_status_to_finished_CDFT, &
            set_status_to_running_CDFT, &
            update_occupations_with_the_maximum_overlap_method, &
            cdft_gs_run_request, &
            gs_run_before_CDFT, scf, single_shot, skip

  !> Type that encapsulates the occupations in a CDFT calculation that have changed w.r.t. a DFT calculation
  type, public :: Occupations
    private
    !> k-point indexes
    integer(i32), allocatable :: i_kpoint(:) 
    !> KS state indexes
    integer(i32), allocatable :: i_state(:)
    !> occupation factors
    real(dp), allocatable :: occ_factor(:)
  contains
    private
    procedure :: allocate_arrays => Occupations_allocate_arrays
    procedure, public :: get_attributes => Occupations_get_attributes
    procedure, public :: get_from_file => Occupations_get_from_file
    procedure :: sanity_checks => Occupations_sanity_checks
    procedure, public :: set_attributes => Occupations_set_attributes
  end type

  !> Type that encapsulates the exciton coefficients
  !> as a structure of arrays
  type, public :: ExcitonCoefficients
    private
    !> k-point indexes
    integer(i32), allocatable :: i_kpoint(:) 
    !> valence-state indexes
    integer(i32), allocatable :: i_vb(:)
    !> Index of the conduction state
    integer(i32), allocatable :: i_cb(:)
    !> Coefficients \(A^{\lambda}_{\mathbf{k}vc}\)
    complex(dp), allocatable  :: coeffs(:) 
  contains
    private
    procedure, public :: set_attributes => ExcitonCoefficients_set_attributes
    procedure, public :: get_attributes => ExcitonCoefficients_get_attributes
    procedure, public :: get_from_file => ExcitonCoefficients_get_from_file
    procedure :: sanity_check => ExcitonCoefficients_sanity_check
    procedure :: get_norm => ExcitonCoefficients_get_normalization
  end type

  !> Enum with the requested GS run in CDFT calculation
  enum, bind(C)
    enumerator :: gs_run_before_CDFT
    enumerator :: single_shot, scf, skip
  end enum

  !> Enum with the status of a CDFT calculation
  enum, bind(C)
    enumerator :: status_CDFT_calculation
    enumerator :: running_GS, running_CDFT, finished_CDFT
  end enum

  !> Type to encapsulate the input elements required for a constrained DFT calculation
  type, public :: cdft_input_keys
    private
    !> If `.true.`, a constrained DFT calculation should be performed
    logical :: on = .false.
    !> If `.true.`, employ the maximum overlap method in the CDFT calculation
    logical :: maximum_overlap_method
    !> If `.true.`, read `STATE.OUT` for the initial guess regarding the electron density and the KS potential
    logical :: start_density_potential_from_file = .false.
    !> If `.true.`, employ the exciton coefficients to occupy the KS states
    logical :: use_exciton_coefficients
    !> If `.true.`, read the occupation numbers or the exciton coefficients from an external file
    logical :: use_external_file
    !> Name of the file containing the exciton coefficients or the occupation numbers
    character(len=:), allocatable :: file_name
    !> Exciton coefficients
    type(ExcitonCoefficients) :: exc
    !> Occupations to be updated (that differ from GS)
    type(Occupations) :: occ
  contains
    private
    procedure, public :: read_input_keys => cdft_input_keys_initialize
    procedure, public :: is_on => cdft_is_on
    procedure, public :: is_maximum_overlap_method_required => cdft_maximum_overlap_method
    procedure, public :: read_density_potential_from_file => cdft_start_density_potential_from_file
    procedure :: cdft_input_keys_mock_with_ExcitonCoefficients, cdft_input_keys_mock_with_Occupations
    generic, public :: mock => cdft_input_keys_mock_with_ExcitonCoefficients, cdft_input_keys_mock_with_Occupations
    procedure :: sanity_check => cdft_sanity_checks
  end type

  integer(kind(status_CDFT_calculation)), save :: status = running_GS
  complex(dp), public, protected, allocatable :: prod(:, :, :)
  complex(dp), public, protected, allocatable :: evecfv_gs(:, :, :)
  
  !> When running CDFT calculation, this file extension should be used in a (previous) 
  !> groundstate calculation to differentiate it from CDFT
  character(len=*), public, parameter :: file_extension_GS = ".OUT"
  !> File extension for the CDFT calculation
  character(len=*), public, parameter :: file_extension_CDFT = "_CDFT.OUT"

contains

!> Set the private variable `status` to the value `running_CDFT`
subroutine set_status_to_running_CDFT()
  status = running_CDFT
end subroutine

!> Set the private variable `status` to the value `finished_CDFT`
subroutine set_status_to_finished_CDFT()
  status = finished_CDFT
end subroutine

!> Returns `.true.`, when a CDFT calculation is running
pure logical function cdft_is_on( cdft_input ) result( is_on )
  class(cdft_input_keys), intent(in) :: cdft_input
  is_on = ( ( cdft_input%on ) .and. ( status == running_CDFT ) )
end function


!> Returns `.true.`, when the maximum overlap method is employed
pure logical function cdft_maximum_overlap_method(this) result(check)
  class(cdft_input_keys), intent(in) :: this

  check = this%is_on()
  if( check ) check = this%maximum_overlap_method
end function

!> Returns `.true.`, if a CDFT calculation is running and `start_density_potential_from_file` is `.true.`
pure logical function cdft_start_density_potential_from_file( cdft_input ) result( check )
  class(cdft_input_keys), intent(in) :: cdft_input
  check = ( cdft_input%is_on() ) .and. ( cdft_input%start_density_potential_from_file ) 
end function

!> Mock cdft_input_keys using ExcitonCoefficients for testing purposes
subroutine cdft_input_keys_mock_with_ExcitonCoefficients( this, exc )
  class(cdft_input_keys), intent(inout) :: this
  !> Exciton coefficients used to build `this`
  type(ExcitonCoefficients), intent(in) :: exc

  this%use_exciton_coefficients = .true.
  this%exc = exc
end subroutine

!> Mock cdft_input_keys using for testing purposes
subroutine cdft_input_keys_mock_with_Occupations( this, occ )
  class(cdft_input_keys), intent(inout) :: this
  !> Occupations used to build `this`
  type(Occupations), intent(in) :: occ

  this%use_exciton_coefficients = .false.
  this%occ = occ
end subroutine

!> Initialize the components of the type [[cdft_input_keys]]
subroutine cdft_input_keys_initialize( this, input_gs )
  class(cdft_input_keys), intent(inout) :: this
  !> Elements and attributes of groundstate defined in the input file
  type(groundstate_type), intent(in) :: input_gs

  integer(i32) :: i

  this%on = associated( input_gs%constrainedDFT )
  if( this%on ) then
    this%maximum_overlap_method = input_gs%constrainedDFT%MaximumOverlapMethod
    this%start_density_potential_from_file = input_gs%constrainedDFT%startDensityAndPotentialFromFile
    if( associated( input_gs%constrainedDFT%occupationChanges ) ) then
      this%use_exciton_coefficients = .false.
      this%use_external_file = .false.
      call this%occ%allocate_arrays( size( input_gs%constrainedDFT%occupationChanges%newOccupationarray ) )
      do i = 1, size( this%occ%i_kpoint )
        this%occ%i_kpoint(i) = input_gs%constrainedDFT%occupationChanges%newOccupationarray(i)%newOccupation%kPointIndex
        this%occ%i_state(i) = input_gs%constrainedDFT%occupationChanges%newOccupationarray(i)%newOccupation%stateIndex
        this%occ%occ_factor(i) = input_gs%constrainedDFT%occupationChanges%newOccupationarray(i)%newOccupation%occupation
      end do
    else
      this%use_exciton_coefficients = input_gs%constrainedDFT%useExcitonCoefficients
      this%use_external_file = input_gs%constrainedDFT%useExternalFile
      this%file_name = trim( input_gs%constrainedDFT%fileName )
      if( this%use_external_file ) then
        if( this%use_exciton_coefficients ) then
          call this%exc%get_from_file( this%file_name )
        else
          call this%occ%get_from_file( this%file_name )
        end if
      end if
    end if
    call this%sanity_check( input_gs )
  end if
end subroutine

!> Sanity checks for CDFT calculations
subroutine cdft_sanity_checks( this, input_gs )
  class(cdft_input_keys), intent(in) :: this
  !> Elements and attributes of groundstate defined in the input file
  type(groundstate_type), intent(in) :: input_gs

  character(len=*), parameter :: partially_compatible_solver = "Davidson"
  character(len=*), parameter :: procedure_name = "cdft_sanity_checks"
  character(len=*), parameter :: warning_header = "Warning(" // procedure_name // "): "

  if ( trim( input_gs%solver%type ) == partially_compatible_solver ) then
    if ( this%is_maximum_overlap_method_required() ) then
      call terminate_if_false( input_gs%solver%constructHS, &
        "Constrained DFT can only be used with " // partially_compatible_solver // " solver with constructHS = true" )
      if ( this%is_maximum_overlap_method_required() ) call warning( warning_header // "Constrained DFT with &
        maximum-overlap method is not yet fully supported when using the " // partially_compatible_solver // " solver. &
        Please treat the results with caution." )
    end if
  end if
  call terminate_if_false( associated(input_gs%constrainedDFT%occupationChanges) .or. this%use_external_file, &
    "Constrained DFT must have the element occupationChanges or use an external file")
end subroutine

!> Allocate arrays of [[Occupations]]
subroutine Occupations_allocate_arrays( this, n )
  class(Occupations), intent(inout) :: this
  !> Size of arrays
  integer(i32), intent(in) :: n

  CALL_ASSERT( n>0, "n must be positive" )
  allocate( this%i_kpoint(n), this%i_state(n), this%occ_factor(n) )
end subroutine

!> Get each attribute of [[Occupations]]
pure subroutine Occupations_get_attributes( occ, kpoint_indexes, state_indexes, occ_factors )
  class(Occupations), intent(in) :: occ
  !> k-point indexes
  integer(i32), allocatable, intent(out) :: kpoint_indexes(:) 
  !> KS state indexes
  integer(i32), allocatable, intent(out) :: state_indexes(:)
  !> occupation factors
  real(dp), allocatable, intent(out) :: occ_factors(:)

  kpoint_indexes = occ%i_kpoint
  state_indexes = occ%i_state
  occ_factors = occ%occ_factor
end subroutine

!> Read a file containing the occupations
subroutine Occupations_get_from_file( this, file_name )
  class(Occupations), intent(inout) :: this
  !> Name of the file to read
  character(len=*), intent(in) :: file_name

  integer(i32) :: i, unit, n_occupations

  open( newunit = unit, file = trim(file_name), action = 'read' )
  read( unit = unit, fmt = * ) n_occupations ! number of (nonzero) exciton coefficients
  call this%allocate_arrays( n_occupations )
  do i = 1, n_occupations
    read( unit=unit, fmt = * ) this%i_state(i), this%i_kpoint(i), this%occ_factor(i)
  end do
  close( unit = unit )
end subroutine

!> Sanity checks
subroutine Occupations_sanity_checks( this, occupation_factors_GS )
  class(Occupations), intent(in) :: this
  !> Occupation factors in ground state
  real(dp), contiguous, intent(in) :: occupation_factors_GS(:, :)

  integer(i32) :: n_kpt, n_states

  n_states = size( occupation_factors_GS, 1 )
  n_kpt = size( occupation_factors_GS, 2 )
  call terminate_if_false( all( this%i_kpoint <= n_kpt ), "Error: trying to force CDFT occupation to k-point out of borders")
  call terminate_if_false( all( this%i_state <= n_states ), "Error: trying to force CDFT occupation to KS state with index out of borders")
  if( any(this%occ_factor < 0._dp ) ) call warning("Forcing negative occupation")
  if( any(this%occ_factor > maxval(occupation_factors_GS) ) ) call warning("Forcing occupations that may be too large")
end subroutine

!> Set each attribute of [[Occupations]]
subroutine Occupations_set_attributes( occ, kpoint_indexes, state_indexes, occ_factor )
  class(Occupations), intent(inout) :: occ
  !> k-point indexes
  integer(i32), contiguous, intent(in) :: kpoint_indexes(:)
  !> state indexes
  integer(i32), contiguous, intent(in) :: state_indexes(:)
  !> array with the occupation factors
  real(dp), contiguous, intent(in) :: occ_factor(:)

  integer(i32) :: n_coeffs

  n_coeffs = size( kpoint_indexes )

  CALL_ASSERT( n_coeffs == size( state_indexes ), "state_indexes must have n_coeffs elements")
  CALL_ASSERT( n_coeffs == size( occ_factor ), "occ_factor must have n_coeffs elements")
  
  occ%i_kpoint = kpoint_indexes
  occ%i_state = state_indexes
  occ%occ_factor = occ_factor
end subroutine

!> Set each attribute of [[ExcitonCoefficient]]
subroutine ExcitonCoefficients_set_attributes( exc_coeff, kpoint_indexes, valence_indexes, conduction_indexes, coeffs )
  class(ExcitonCoefficients), intent(inout) :: exc_coeff
  !> k-point indexes
  integer(i32), contiguous, intent(in) :: kpoint_indexes(:)
  !> valence-band indexes
  integer(i32), contiguous, intent(in) :: valence_indexes(:)
  !> conduction-band indexes
  integer(i32), contiguous, intent(in) :: conduction_indexes(:)
  !> array with the exciton coefficients \(A^{\lambda}_{\mathbf{k}vc}\)
  complex(dp), contiguous, intent(in) :: coeffs(:)

  integer(i32) :: n_coeffs

  n_coeffs = size( kpoint_indexes )

  CALL_ASSERT( n_coeffs == size( valence_indexes ), "valence_indexes must have n_coeffs elements")
  CALL_ASSERT( n_coeffs == size( conduction_indexes ), "conduction_indexes must have n_coeffs elements")
  CALL_ASSERT( n_coeffs == size( coeffs ), "coeffs must have n_coeffs elements")
  
  exc_coeff%i_kpoint = kpoint_indexes
  exc_coeff%i_vb = valence_indexes
  exc_coeff%i_cb = conduction_indexes
  exc_coeff%coeffs = coeffs
end subroutine

!> Get each attribute of [[ExcitonCoefficient]]
pure subroutine ExcitonCoefficients_get_attributes( exc_coeff, kpoint_indexes, valence_indexes, conduction_indexes, coeffs )
  class(ExcitonCoefficients), intent(in) :: exc_coeff
  !> k-point indexes
  integer(i32), allocatable, intent(out) :: kpoint_indexes(:)
  !> valence-band indexes
  integer(i32), allocatable, intent(out) :: valence_indexes(:)
  !> conduction-band indexes
  integer(i32), allocatable, intent(out) :: conduction_indexes(:)
  !> array with the exciton coefficients \(A^{\lambda}_{\mathbf{k}vc}\)
  complex(dp), allocatable, intent(out) :: coeffs(:)

  kpoint_indexes = exc_coeff%i_kpoint
  valence_indexes = exc_coeff%i_vb
  conduction_indexes = exc_coeff%i_cb
  coeffs = exc_coeff%coeffs
end subroutine

!> Read a file containing the exciton coefficients
subroutine ExcitonCoefficients_get_from_file( exc_coeffs, file_name )
  class(ExcitonCoefficients), intent(inout) :: exc_coeffs
  !> Name of the file to read
  character(len=*), intent(in) :: file_name

  integer(i32) :: i, unit, n_coeffs
  integer(i32), allocatable :: iv(:), ic(:), ik(:)
  real(dp) :: re, im
  complex(dp), allocatable :: z(:)

  open( newunit = unit, file = trim(file_name), action = 'read' )
  read( unit = unit, fmt = * ) n_coeffs ! number of (nonzero) exciton coefficients
  allocate( iv(n_coeffs), ic(n_coeffs), ik(n_coeffs), z(n_coeffs) )
  do i = 1, n_coeffs
    read( unit=unit, fmt = * ) iv(i), ic(i), ik(i), re, im
    z(i) = cmplx(re, im, kind=dp)
  end do
  call exc_coeffs%set_attributes( ik, iv, ic, z )
  close( unit = unit )
end subroutine


subroutine initialize_cdft_global_arrays( psi_gs, first_k )
  !> index of the first k-point
  integer(i32), intent(in) :: first_k
  !> Groundstate KS wavefunctions
  complex(dp), contiguous, intent(in) :: psi_gs(:, :, first_k:)

  call deallocate_cdft_global_arrays()
  if (size( psi_gs, 3) /= 0) then
    associate( n_basis => size( psi_gs, 1 ), n_states => size( psi_gs, 2 ), last_k => ubound( psi_gs, 3 ) )
      allocate( prod(n_basis, n_states, first_k:last_k), source=zzero )
      allocate( evecfv_gs(n_basis, n_states, first_k:last_k), source=psi_gs )
    end associate
  else
    associate( n_basis => size( psi_gs, 1 ), n_states => size( psi_gs, 2 ))
      allocate( prod(n_basis, n_states, 0) )
      allocate( evecfv_gs(n_basis, n_states, 0) )
    end associate
  end if
end subroutine


subroutine set_overlap_times_psi_gs( ik, S )
  integer(i32), intent(in) :: ik
  complex(dp), contiguous, intent(in) :: S(:, :)

  integer(i32) :: n

  n = size( S , 1 )
  CALL_ASSERT( size( evecfv_gs, 1 ) >= n, "S is not compatible with evecfv_gs" )
  call hermitian_matrix_multiply(S, evecfv_gs(1:n, :, ik), prod(1:n, :, ik))
end subroutine


!> Deallocate the global arrays defined in this module
subroutine deallocate_cdft_global_arrays()
  if( allocated( prod ) ) deallocate( prod )
  if( allocated( evecfv_gs) ) deallocate( evecfv_gs )
end subroutine


!> Sanity checks for the exciton coefficients
subroutine ExcitonCoefficients_sanity_check( exc, occupations_GS )
  !> Array with the exciton coefficients
  class(ExcitonCoefficients), intent(in) :: exc
  !> Occupation factors of each KS state obtained in a groundstate calculation
  real(dp), contiguous, intent(in) :: occupations_GS(:, :)

  integer(i32) :: i, n_coeffs, n_kpoints, n_states
  real(dp) :: exc_normalization
  real(dp), parameter :: tol = 1e-6_dp

  n_coeffs = size( exc%coeffs )
  n_states = size( occupations_GS, 1 )
  n_kpoints = size( occupations_GS, 2 )
  ! Hard conditions
  call terminate_if_false( all( exc%i_kpoint <= n_kpoints ), &
    "One exciton coefficient has k-point index which is larger than n-kpoints: " // to_char(n_kpoints) )
  call terminate_if_false( all( exc%i_vb <= n_states ), &
    "One exciton coefficient has valence-band index which is larger than n_states: " // to_char(n_states) )
  call terminate_if_false( all( exc%i_cb <= n_states ), &
    "One exciton coefficient has conduction-band index which is larger than n_states: " // to_char(n_states) )
  ! Soft conditions
  do i = 1, n_coeffs
    if( occupations_GS(exc%i_vb(i), exc%i_kpoint(i)) <= tol ) &
      call warning("Exciton coefficient i: " // to_char(i) // " takes state " // to_char(exc%i_vb(i)) // " as a valence state, but it is unoccupied" )
    if( occupations_GS(exc%i_cb(i), exc%i_kpoint(i)) > tol ) &
      call warning("Exciton coefficient i: " // to_char(i) // " takes state " // to_char(exc%i_cb(i)) // " as a conduction state, but it is occupied" )
  end do
  exc_normalization = exc%get_norm()
  if( .not. all_close(exc_normalization, real_one, tol ) ) call warning("Exciton coefficients are not normalized to 1: " // to_char(exc_normalization) )
end subroutine


! (private function)
pure real(dp) function ExcitonCoefficients_get_normalization( exc ) result( x )
  class(ExcitonCoefficients), intent(in) :: exc

  x = sqrt( real( dot_product( exc%coeffs, exc%coeffs ), kind=dp ) )
end function

!> Determine the occupations to be constrained during the SCF cycle
subroutine determine_cdft_occupations( cdft_inp, kpt_weights, occupation_factors )
  !> CDFT input keys 
  type(cdft_input_keys), intent(in) :: cdft_inp
  !> k-point weights
  real(dp), contiguous, intent(in) :: kpt_weights(:)
  !> occupation factors
  real(dp), contiguous, intent(inout) :: occupation_factors(:, :)

  if( cdft_inp%use_exciton_coefficients ) then
    call cdft_inp%exc%sanity_check( occupation_factors )
    call occupy_cdft_from_exc_coeff( cdft_inp%exc, kpt_weights, occupation_factors )
  else
    call cdft_inp%occ%sanity_checks( occupation_factors )
    call occupy_cdft_from_forced_occs( cdft_inp%occ, occupation_factors )
  end if
end subroutine

!> Determine the occupation factors in a CDFT calculation using exciton coefficients
pure subroutine occupy_cdft_from_exc_coeff( exc, wkpt, occupation_factors )
  !> Type with exciton coefficients
  type(ExcitonCoefficients), intent(in)  :: exc
  !> The weight of each k-point
  real(dp), contiguous, intent(in) :: wkpt(:)
  !> The occupation factors of each KS state and each k-point
  real(dp), contiguous, intent(inout) :: occupation_factors(:, :)

  integer(i32)          :: i, ik, number_elements

  number_elements = size( exc%coeffs )
  do i = 1, number_elements
    ik = exc%i_kpoint(i)
    call change_occupations( exc%i_vb(i), exc%i_cb(i), ( abs( exc%coeffs(i) ) ** 2 )/wkpt(ik), occupation_factors(:, ik) )
  end do
end subroutine

!> Determine the occupation factors in a CDFT calculation using constrained occupation factors
pure subroutine occupy_cdft_from_forced_occs( occ, occupation_factors )
  !> Type containing occupations to constrain
  type(Occupations), intent(in)  :: occ
  !> Occupation factors of each KS state and each k-point
  real(dp), contiguous, intent(inout) :: occupation_factors(:, :)

  integer(i32)          :: i, ik, is, number_elements

  number_elements = size( occ%i_kpoint )
  do i = 1, number_elements
    ik = occ%i_kpoint(i)
    is = occ%i_state(i)
    occupation_factors(is, ik) = occ%occ_factor(i)
  end do
end subroutine

!> Update the occupation numbers following the maximum overlap method
subroutine update_occupations_with_the_maximum_overlap_method( psi, occupation_factors )
  !> Wavefunction coefficients (in terms of the LAPW+LO basis)
  complex(dp), intent(in) :: psi(:, :, :)
  !> Occupation factors
  real(dp), intent(inout) :: occupation_factors(:, :)

  integer(i32) :: ik, first_k, i, idx_max, n_states
  real(dp), allocatable :: occ_save(:)
  complex(dp), allocatable :: projection(:, :)
  logical, allocatable :: search(:)

  first_k = lbound( prod, 3 ) 
  n_states = size( psi, 2 )
  CALL_ASSERT( size(occupation_factors, 1) == n_states, "occupation_factors must have n_states elements along 1st dim")
  CALL_ASSERT( size(psi, 3) == size(occupation_factors, 2), "occupation_factors and psi have incompatible size")
  allocate( projection(n_states, n_states), occ_save(n_states), search(n_states) )
  do ik = 1, size( psi, 3 )
    call matrix_multiply( prod(:, :, ik+first_k-1), psi(:, :, ik), projection, 'C', 'N' )
    occ_save = occupation_factors(:, ik)
    search = .true.
    do i = 1, n_states
      idx_max = maxloc( abs(projection(:, i)), mask=search, dim=1 )
      search(idx_max) = .false.
      occupation_factors(i, ik) = occ_save(idx_max)
    end do 
  end do 
end subroutine

! (private subroutine)
pure subroutine change_occupations( index_vb, index_cb, delta, occ )
  integer(i32), intent(in) :: index_vb
  integer(i32), intent(in) :: index_cb
  real(dp), intent(in) :: delta
  real(dp), intent(inout) :: occ(:)

  occ(index_vb) = occ(index_vb) - delta
  occ(index_cb) = occ(index_cb) + delta
end subroutine

!> Given a string, get the corresponding [[gs_run_before_CDFT]]
function cdft_gs_run_request( string ) result(r)
  !> String containing the pre-cDFT GS run request
  character(len=*), intent(in) :: string
  integer(kind( gs_run_before_CDFT )) :: r

  select case ( trim( string ) )
    case ("singleShot")
      r = single_shot
    case ("scf")
      r = scf
    case ("skip")
      r = skip
    case default
      CALL_ASSERT( .false., "Unrecognized cdft_gs_run_request" )
  end select
end function

end module
