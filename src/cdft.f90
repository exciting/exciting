!> Module for constrained DFT (CDFT) calculations
!> Created Oct 2024 (Ronaldo)
module cdft
  use asserts, only: assert
  use constants, only: real_one, real_zero, zzero
  use math_utils, only: all_close
  use mod_mpi_env, only: mpiinfo
  use modinput, only: groundstate_type
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32
  use to_char_conversion, only: to_char

  implicit none

  private

  public :: occupy_cdft, &
            set_status_to_finished_CDFT, &
            set_status_to_running_CDFT

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
    procedure :: set_attributes => ExcitonCoefficients_set_attributes
    procedure :: get_attributes => ExcitonCoefficients_get_attributes
    procedure :: get_from_file => ExcitonCoefficients_get_from_file
    procedure :: sanity_check => ExcitonCoefficients_sanity_check
    procedure, private :: get_norm => ExcitonCoefficients_get_normalization
  end type

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
    !> If `.true.`, read `STATE.OUT` for the initial guess regarding the electron density and the KS potential
    logical :: start_density_potential_from_file = .false.
    !> If `.true.`, employ the exciton coefficients to occupy the KS states
    logical :: use_exciton_coefficients
    !> If `.true.`, read the occupation numbers or the exciton coefficients from an external file
    logical :: use_external_file
    !> Name of the file containing the exciton coefficients or the occupation numbers
    character(len=:), allocatable, public :: file_name
  contains
    procedure :: read_input_keys => cdft_input_keys_initialize
    procedure :: is_on => cdft_is_on
    procedure :: read_density_potential_from_file => cdft_start_density_potential_from_file
    procedure, private :: sanity_check => cdft_sanity_checks
  end type

  integer(kind(status_CDFT_calculation)), save :: status = running_GS
  
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


!> Returns `.true.`, if a CDFT calculation is running and `start_density_potential_from_file` is `.true.`
pure logical function cdft_start_density_potential_from_file( cdft_input ) result( check )
  class(cdft_input_keys), intent(in) :: cdft_input
  check = ( cdft_input%is_on() ) .and. ( cdft_input%start_density_potential_from_file ) 
end function


!> Initialize the components of the type [[cdft_input_keys]]
subroutine cdft_input_keys_initialize( this, input_gs )
  class(cdft_input_keys), intent(inout) :: this
  !> Elements and attributes of groundstate defined in the input file
  type(groundstate_type), intent(in) :: input_gs

  this%on = associated( input_gs%constrainedDFT )
  if( this%on ) then
    this%start_density_potential_from_file = input_gs%constrainedDFT%startDensityAndPotentialFromFile
    this%use_exciton_coefficients = input_gs%constrainedDFT%useExcitonCoefficients
    this%use_external_file = input_gs%constrainedDFT%useExternalFile
    this%file_name = trim( input_gs%constrainedDFT%fileName )
    call this%sanity_check( input_gs )
  end if
end subroutine


!> Sanity checks for CDFT calculations
subroutine cdft_sanity_checks( this, input_gs )
  class(cdft_input_keys), intent(in) :: this
  !> Elements and attributes of groundstate defined in the input file
  type(groundstate_type), intent(in) :: input_gs

  character(len=*), parameter :: compatible_solver = "Lapack"

  call terminate_if_false( trim( input_gs%solver%type ) == compatible_solver, &
    "Constrained DFT currently only implemented for solver " // compatible_solver )
  call terminate_if_false( this%use_exciton_coefficients, &
    "Constrained DFT currently only implemented for useExcitonCoefficients true" )
  call terminate_if_false( this%use_external_file, &
    "Constrained DFT currently only implemented for useExternalFile true" )
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

  call assert( n_coeffs == size( valence_indexes ), "valence_indexes must have n_coeffs elements")
  call assert( n_coeffs == size( conduction_indexes ), "conduction_indexes must have n_coeffs elements")
  call assert( n_coeffs == size( coeffs ), "coeffs must have n_coeffs elements")
  
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

  x = sqrt( dot_product( exc%coeffs, exc%coeffs ) )
end function


!> Evaluate the occupation factors in a CDFT calculation
subroutine occupy_cdft( exc, wkpt, occupation_factors )
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


! (private subroutine)
pure subroutine change_occupations( index_vb, index_cb, delta, occ )
  integer(i32), intent(in) :: index_vb
  integer(i32), intent(in) :: index_cb
  real(dp), intent(in) :: delta
  real(dp), intent(inout) :: occ(:)

  occ(index_vb) = occ(index_vb) - delta
  occ(index_cb) = occ(index_cb) + delta
end subroutine

end module
