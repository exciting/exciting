!> Module with propagators for real-time TDDFT calculations
module propagators
#include "asserts.fpp"
  use constants, only: zi
  use integration, only: RungeKutta4thOrder => ODESolver_RungeKutta4thOrder
  use matrix_exp, only: &
      & exp_hermitian => exp_hermitian_matrix_times_vectors, & 
      & exp_general => exp_general_matrix_times_vectors, &
      & exphouston_propagator => exphouston_hermitian_matrix_times_vectors
  use precision, only: dp, i32

  implicit none

  private

  public :: create_propagator, propagator_input_elements

  !> Enum with the solver type for the propagator
  enum, bind(C)
    enumerator :: propagator_methods
    enumerator :: SE, EMR, AETRS, CFM4, EH, EHM, RK4
  end enum

  !> A type to encapsulate data relative to the propagators that
  !> are given in the input file
  type :: propagator_input_elements
    private
    !> Propagator method used
    integer(kind( propagator_methods )) :: method
    !> Time step \( \Delta t \) employed in RT-TDDFT
    real(dp) :: dt_
    !> Order of the Taylor expansion
    integer(i32) :: order_taylor
    !> Tolerance required for diagonalization
    real(dp) :: tol_diagonalization
    !> Number of eigenvectors used in the Houston propagator expansion
    integer(i32) :: n_eigvecs_houston
  contains
    procedure, public :: initialize => initialize_propagator_input_elements
  end type

  !> Abstract type that should be extended by any concrete propagator
  type, public, abstract :: propagator
    private
    !> Time step \(\Delta t\)
    real(dp) :: dt
  contains
    procedure, public :: evolve => propagate_list_of_arrays
    procedure(propagate_single_array_), private, deferred :: propagate_single_array
    procedure(initialize_), public, deferred :: initialize
    procedure, public :: extrapolation_needed => propagator_requires_extrapolation
    procedure, public :: update_and_check => set_and_check_n_eigvecs_houston
    procedure, public :: time_step => propagator_time_step
  end type

  !> Abstract type for propagators that employ a Taylor expansion
  type, abstract, extends(propagator) :: propagator_with_Taylor_expansion
    private
    !> The order of the Taylor expansion
    integer(i32) :: order_Taylor
    !> Procedure pointer that, in the initialization, will point either to the
    !> exponential of either generic or Hermitian matrices
    procedure(exp_general), pointer, nopass :: exp_matrix => null()
  contains
    procedure, public :: initialize => initialize_Taylor
  end type

  !> Type to encapsulate the `SE` propagator 
  type, extends (propagator_with_Taylor_expansion) :: SE_propagator 
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_SE
  end type 

  !> Type to encapsulate the `EMR` propagator 
  type, extends (propagator_with_Taylor_expansion) :: EMR_propagator 
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_EMR
  end type 

  !> Type to encapsulate the `AETRS` propagator 
  type, extends (propagator_with_Taylor_expansion) :: AETRS_propagator 
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_AETRS
  end type 

  !> Type to encapsulate the `CFM4` propagator 
  type, extends (propagator_with_Taylor_expansion) :: CFM4_propagator 
    private
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_CFM4
  end type 

  !> Type to encapsulate the `RK4` propagator 
  type, extends (propagator) :: RK4_propagator 
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_RK4
    procedure, public :: initialize => initialize_RK4 
  end type 

  !> Abstract type for propagators that employ the Houston method: the
  !> exponential of a matrix is evaluated in an exact form through its 
  !> eigenvalues and eigenvectors
  type, abstract, extends(propagator) :: Houston_propagator
    private
    !> Tolerance used to diagonalize the matrix to be exponentiated
    real(dp) :: tol
    !> Number of eigenvectors used in the Houston propagator expansion
    integer(i32) :: n_eigvecs_houston
  contains
    procedure, public :: initialize => initialize_Houston
  end type

  !> Type to encapsulate the `EH` propagator 
  type, extends (Houston_propagator) :: EH_propagator 
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_EH
  end type 

  !> Type to encapsulate the `EHM` propagator 
  type, extends (Houston_propagator) :: EHM_propagator 
  contains
    procedure, private :: propagate_single_array => propagate_single_array_with_EHM
  end type 

  abstract interface
    subroutine propagate_single_array_( self, dim, H_dt, H_0, H_minus_dt, S, x )
      import :: propagator, dp, i32
      class(propagator) :: self      
      !> Actual dimensions of each matrix: `H_dt`, `H_0`, and `S`
      integer(i32), intent(in) :: dim
      !> Hamiltonian matrix \(H\) at time \(\Delta t\)
      complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
      !> Hamiltonian matrix \(H\) at time \(0\)
      complex(dp), contiguous, intent(in) :: H_0(:, :)
      !> Hamiltonian matrix \(H\) at time \(-\Delta t\)
      complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
      !> Overlap matrix \(S\)
      complex(dp), contiguous, intent(in) :: S(:, :)
      !> In: wavefunction \(\psi(0)\). Out: wavefunction \(\psi(\Delta t)\)
      complex(dp), contiguous, intent(inout) :: x(:, :)
    end subroutine
    subroutine initialize_( this, input_parameters, is_hermitian )
      import :: propagator, propagator_input_elements, dp
      class(propagator), intent(inout) :: this
      !> Type that encapsulates the input parameters relative to propagators
      type(propagator_input_elements), intent(in) :: input_parameters
      !> If `.true.`, the matrix to be exponentiated is hermitian
      logical, intent(in) :: is_hermitian
    end subroutine
  end interface 

  !> Trick to mimic a C++ constructor (applied to the type [[propagator_input_elements]])
  interface propagator_input_elements
    module procedure :: constructor_propagator_input_elements
  end interface

contains
  !> Get the enum corresponding to the string
  function propagator_type( string ) result( method )
    character(len=*), intent(in) :: string
    integer(kind(propagator_methods)) :: method

    select case( trim(string) )
      case('SE')
        method = SE
      case('EMR')
        method = EMR
      case('AETRS')
        method = AETRS
      case('CFM4')
        method = CFM4
      case('EH')
        method = EH
      case('EHM')
        method = EHM
      case('RK4')
        method = RK4
      case default
        CALL_ASSERT( .false., 'Unrecognized propagator method' )
    end select
  end function

  !> Return `.true.` if a propagator requires extrapolation of \(H(t)\)
  pure logical function propagator_requires_extrapolation( this ) result(check)
    class(propagator), intent(in) :: this
    
    select type( this )
      type is (SE_propagator)
        check = .false.
      type is (EH_propagator)
        check = .false.
      class default
        check = .true.
    end select
  end function

  !> Subroutine to initialize the components of [[propagator_input_elements]]
  subroutine initialize_propagator_input_elements( self, method, dt, order_Taylor, tol, n_eigvecs_houston )
    class(propagator_input_elements), intent(inout) :: self
    !> String containing the name of the propagator following [[propagator_type]]
    character(len=*), intent(in) :: method
    !> Time step
    real(dp), intent(in) :: dt
    !> The order used in a Taylor expansion
    integer(i32), intent(in) :: order_Taylor
    !> Tolerance used in the diagonalization
    real(dp), intent(in) :: tol
    !> Number of eigenvectors used in the Houston propagator expansion
    integer(i32) :: n_eigvecs_houston

    self%method = propagator_type( method )
    self%dt_ = dt
    self%order_Taylor = order_Taylor
    self%tol_diagonalization = tol
    self%n_eigvecs_houston = n_eigvecs_houston
  end subroutine

  !> Subroutine to set the value of n_eigvecs_houston to the default one 
  !> in case negative value is provided by user, and checks whether it 
  !> lies in range from  n_eigvecs_houston_min to n_eigvecs_houston_max. Returns 
  !> success = .True. if the check was succsessful.
  subroutine set_and_check_n_eigvecs_houston( self, n_eigvecs_houston_min, &
    n_eigvecs_houston_max, n_eigvecs_houston_default, success )
    class(propagator), intent(inout) :: self
    !> Minimal adequate value of n_eigvecs_houston
    integer(i32), intent(in) :: n_eigvecs_houston_min
    !> Maximal adequate value of n_eigvecs_houston
    integer(i32), intent(in) :: n_eigvecs_houston_max
    !> Default value of n_eigvecs_houston which should be used if a user provides 
    !> n_eigvecs_houston < 0
    integer(i32), intent(in) :: n_eigvecs_houston_default
    !> If .true., propagation parameters are consistent with each other
    logical, intent(out) :: success

    success = .true.
    select type( self )
    class is ( Houston_propagator )
      if ( self%n_eigvecs_houston < 0 ) self%n_eigvecs_houston = n_eigvecs_houston_default
      if ( self%n_eigvecs_houston < n_eigvecs_houston_min .or. &
        self%n_eigvecs_houston > n_eigvecs_houston_max ) success = .false.
    end select

  end subroutine

  !> Return the component `dt` of [[propagator]]
  pure real(dp) function propagator_time_step( self ) result( time_step )
    class(propagator), intent(in) :: self
    time_step = self%dt
  end function

  !> Return a variable of type [[propagator_input_elements]], using all its attributes.
  !> This is a trick to mimic a C++ constructor (applied here to the type [[propagator_input_elements]])
  function constructor_propagator_input_elements( method, dt, order_Taylor, tol, n_eigvecs_houston ) result(params)
    type(propagator_input_elements) :: params
    !> String containing the propagator method to be used
    character(len=*), intent(in) :: method
    !> Time step \( \Delta t \) employed in RT-TDDFT
    real(dp), intent(in) :: dt
    !> Order of the Taylor expansion
    integer(i32), intent(in) :: order_Taylor
    !> Tolerance required for diagonalization
    real(dp), intent(in) :: tol
    !> Number of eigenvectors used in the Houston propagator expansion
    integer(i32) :: n_eigvecs_houston

    call params%initialize( method, dt, order_Taylor, tol, n_eigvecs_houston )
  end function

  !> Initializes the propagator class ([[propagator]])
  !> with a concrete propagator, depending on the propagator method contained in 
  !> `input_parameters`
  subroutine create_propagator( prop, input_parameters, is_hermitian )
    class(propagator), allocatable, intent(out) :: prop
    !> Type that encapsulates the input parameters relative to propagators
    type(propagator_input_elements), intent(in) :: input_parameters
    !> If `.true.`, the matrix to be exponentiated is hermitian
    logical, intent(in) :: is_hermitian

    select case( input_parameters%method )
    case( SE )
      allocate( SE_propagator :: prop )
    case( EMR )
      allocate( EMR_propagator :: prop )
    case( AETRS )
      allocate( AETRS_propagator :: prop )
    case( CFM4 ) 
      allocate( CFM4_propagator :: prop )
    case( RK4 ) 
      allocate( RK4_propagator :: prop)
    case( EH ) 
      allocate( EH_propagator :: prop )
    case( EHM ) 
      allocate( EHM_propagator :: prop )
    case default
      CALL_ASSERT( .false., 'propagator method not recognized')
    end select
    call prop%initialize( input_parameters, is_hermitian )
  end subroutine

  !> Initialize a propagator of class [[propagator_with_Taylor_expansion]]. 
  !> The arguments are documented in [[initialize_]]
  subroutine initialize_Taylor( this, input_parameters, is_hermitian )
    class(propagator_with_Taylor_expansion), intent(inout) :: this
    type(propagator_input_elements), intent(in) :: input_parameters
    logical, intent(in) :: is_hermitian

    this%dt = input_parameters%dt_
    this%order_Taylor = input_parameters%order_Taylor
    if( is_hermitian ) then
      this%exp_matrix => exp_hermitian
    else
      this%exp_matrix => exp_general
    end if
  end subroutine

  !> Initialize a propagator of class [[RK4_propagator]].
  !> The arguments are documented in [[initialize_]]
  subroutine initialize_RK4( this, input_parameters, is_hermitian )
    class(RK4_propagator), intent(inout) :: this
    type(propagator_input_elements), intent(in) :: input_parameters
    logical, intent(in) :: is_hermitian

    this%dt = input_parameters%dt_
  end subroutine

  !> Initialize a propagator of class [[Houston_propagator]].
  !> The arguments are documented in [[initialize_]]
  subroutine initialize_Houston( this, input_parameters, is_hermitian )
    class(Houston_propagator), intent(inout) :: this
    type(propagator_input_elements), intent(in) :: input_parameters
    logical, intent(in) :: is_hermitian

    this%dt = input_parameters%dt_
    this%tol = input_parameters%tol_diagonalization
    this%n_eigvecs_houston = input_parameters%n_eigvecs_houston
  end subroutine

  !> Propagate a list of arrays. The 3rd dimension goes, usually, over the various k-points.
  subroutine propagate_list_of_arrays( self, list_of_H_dt, list_of_H_0, list_of_H_minus_dt, list_of_S, psi, dims )
    class(propagator) :: self
    !> The list of \( H \) matrices at time \(\Delta t\)
    complex(dp), contiguous, optional, intent(in) :: list_of_H_dt(:, :, :)
    !> The list of \( H \) matrices at time \(0\)
    complex(dp), contiguous, intent(in) :: list_of_H_0(:, :, :)
    !> The list of \( H \) matrices at time \(-\Delta t\)
    complex(dp), contiguous, optional, intent(in) :: list_of_H_minus_dt(:, :, :)
    !> The list of overlap matrices \(S\)
    complex(dp), contiguous, intent(in) :: list_of_S(:, :, :)
    !> The list of wavefunctions to be propagated
    complex(dp), contiguous, intent(inout) :: psi(:, :, :)
    !> Actual dimensions of each matrix in `list_of_H_dt`, `list_of_H_0`, and `list_of_S`
    integer(i32), intent(in), optional :: dims(:)

    integer(i32) :: i, m, case_H
    integer(i32), allocatable :: dims_(:)

    enum, bind(C)
      enumerator :: H0_only, H0_and_Hdt, H0_and_Hminusdt
    end enum
    
    m = size( list_of_H_0, 3 )
    case_H = H0_only
    if( present(list_of_H_dt) ) then
      case_H = H0_and_Hdt
      CALL_ASSERT( size(list_of_H_dt, 3) == m, 'matrix has 3rd dim different from m' )
      CALL_ASSERT( .not. present(list_of_H_minus_dt), "both H_dt and H_minus_dt cannot be passed at the same time")
    end if
    if( present(list_of_H_minus_dt) ) then
      case_H = H0_and_Hminusdt
      CALL_ASSERT( size(list_of_H_minus_dt, 3) == m, 'matrix has 3rd dim different from m' )
    end if
    CALL_ASSERT( size(list_of_S, 3) == m, 'matrix has 3rd dim different from m' )
    CALL_ASSERT( size(psi, 3) == m, 'matrix has 3rd dim different from m' )
    if( present(dims) ) then
      CALL_ASSERT( size(dims) == m, 'array must have m elements')
      dims_ = dims
    else
      dims_ = spread( size( list_of_H_0, 1 ), dim=1, ncopies=m )
    end if
    select case(case_H)
      case(H0_and_Hminusdt)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(m, self, psi, list_of_H_minus_dt, list_of_H_0, list_of_S, dims_)
      do i = 1, m
        call self%propagate_single_array( dims_(i), H_0=list_of_H_0(:, :, i), H_minus_dt=list_of_H_minus_dt(:, :, i), S=list_of_S(:, :, i), x=psi(:, :, i) )
      end do
!$OMP END PARALLEL DO
      case(H0_only)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(m, self, psi, list_of_H_0, list_of_S, dims_)
      do i = 1, m
        call self%propagate_single_array( dims_(i), H_0=list_of_H_0(:, :, i), S=list_of_S(:, :, i), x=psi(:, :, i) )
      end do
!$OMP END PARALLEL DO
      case(H0_and_Hdt)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(m, self, psi, list_of_H_dt, list_of_H_0, list_of_S, dims_)
      do i = 1, m
        call self%propagate_single_array( dims_(i), H_dt=list_of_H_dt(:, :, i), H_0=list_of_H_0(:, :, i), S=list_of_S(:, :, i), x=psi(:, :, i) )
      end do
!$OMP END PARALLEL DO
    end select
    
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[SE]]
  !> \[ \psi(\Delta t) = \mathrm{exp}(-\mathrm{i}\Delta t S^{-1}H(0)) \psi(0).\]
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_SE( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(SE_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)
    
    CALL_ASSERT( associated(self%exp_matrix), 'exp_matrix not associated')
    call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_0, S, x)
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[EMR]]
  !> \[ \psi(\Delta t) = \mathrm{exp}(-\mathrm{i}\Delta t S^{-1}H(\Delta t/2)) \psi(0).\]
  !> \(H(\Delta t/2)\) is estimated as
  !> \[ H(\Delta t/2) = \frac{1}{2}\left[ H(\Delta t) + H(0)\right].\]
  !> or as
  !> \[ H(\Delta t/2) = H(0) + \frac{1}{2}\left[ H(0) - H(-\Delta t)\right].\]
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_EMR( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(EMR_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)
    
    complex(dp), allocatable :: H_aux(:, :)

    CALL_ASSERT( present(H_dt) .neqv. present(H_minus_dt), "only one optional argument must be present" )
    CALL_ASSERT( associated(self%exp_matrix), 'exp_matrix not associated')
    if( present(H_dt) ) then
      H_aux = 0.5_dp*( H_dt + H_0 )
    else
      H_aux = 0.5_dp * (3._dp * H_0 - H_minus_dt)
    end if
    call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[AETRS]]
  !> \[ \psi(\Delta t) = \mathrm{exp}(-\mathrm{i}\Delta t S^{-1}H(\Delta t)/2)
  !> \mathrm{exp}(-\mathrm{i}\Delta t S^{-1}H(0)/2) \psi(0).\]
  !> If needed, \(H(\Delta t)\) is obtained as
  !> \[H(\Delta t) = H(0) + (H(0)-H(-\Delta t))\]
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_AETRS( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(AETRS_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)

    complex(dp), allocatable :: H_aux(:, :)

    CALL_ASSERT( present(H_dt) .neqv. present(H_minus_dt), "only one optional argument must be present" )
    CALL_ASSERT( associated(self%exp_matrix), 'exp_matrix not associated')
    allocate( H_aux, mold=H_0 )
    H_aux = 0.5_dp*H_0
    call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
    if( present(H_dt) ) then
      H_aux = 0.5_dp*H_dt
    else
      H_aux = H_0 - 0.5_dp*H_minus_dt
    end if
    call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[CFM4]]
  !> \[ \psi(\Delta t) = 
  !> \mathrm{exp}[-\mathrm{i}\Delta t S^{-1}(a_1H(\Delta t_1)+a_2H(\Delta t_2))]
  !> \mathrm{exp}[-\mathrm{i}\Delta t S^{-1}(a_2H(\Delta t_1)+a_1H(\Delta t_2))] 
  !>     \psi(0),\]
  !> where \(a_1 = 1/4 - \sqrt{3}/6\), and \(a_2 = 1/4 + \sqrt{3}/6\). The times
  !> \(\Delta t_1 = f_1 \Delta t \), \(\Delta t_2 = f_2 \Delta t \), where 
  !> \( f_1 = 1/2 - \sqrt{3}/6 \) and \( f_2 = 1/2 + \sqrt{3}/6\). Furthermore
  !> \(H(f\Delta t)\) is estimated as
  !> \[ H(f\Delta t) = fH(\Delta t) + (1-f)H(0)\]
  !> or as
  !> \[ H(f\Delta t) = H(0) + f[H(0)-H(-\Delta t)] = (1+f)H(0)-fH(-\Delta t)\]
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_CFM4( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(CFM4_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)

    ! Factors that multiply the hamiltonian in the following propagator:
    ! Commutator-Free Magnus expansion of 4th order
    real(dp), parameter       :: f1 =  0.5_dp - sqrt(3._dp) / 6._dp ! 1/2 - sqrt(3)/6
    real(dp), parameter       :: f2 =  0.5_dp + sqrt(3._dp) / 6._dp ! 1/2 + sqrt(3)/6
    real(dp), parameter       :: a1 =  0.25_dp - sqrt(3._dp) / 6._dp ! 1/4 - sqrt(3)/6
    real(dp), parameter       :: a2 =  0.25_dp + sqrt(3._dp) / 6._dp ! 1/4 + sqrt(3)/6
    real(dp), parameter       :: b_0 = a1*(1-f2) + a2*(1-f1)
    real(dp), parameter       :: b_dt = a1*f2 + a2*f1
    real(dp), parameter       :: c_0 = a1*(1-f1) + a2*(1-f2)
    real(dp), parameter       :: c_dt = a1*f1 + a2*f2
    real(dp), parameter       :: d_0 = a1*(1+f2) + a2*(1+f1)
    real(dp), parameter       :: d_minus_dt = - (a1*f2 + a2*f1)
    real(dp), parameter       :: e_0 = a1*(1+f1) + a2*(1+f2)
    real(dp), parameter       :: e_minus_dt = -(a1*f1 + a2*f2)
    
    complex(dp), allocatable :: H_aux(:, :)

    CALL_ASSERT( present(H_dt) .neqv. present(H_minus_dt), "only one optional argument must be present" )
    CALL_ASSERT( associated(self%exp_matrix), 'exp_matrix not associated')
    if( present(H_dt) ) then
      ! a1*( (1-f2)*H_0 + f2*H_dt ) + a2*( (1-f1)*H_0 + f1*H_dt )
      H_aux = b_0*H_0 + b_dt*H_dt
      call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
      ! a1*( (1-f1)*H_0 + f1*H_dt ) + a2*( (1-f2)*H_0 + f2*H_dt )
      H_aux = c_0*H_0 + c_dt*H_dt
      call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
    else
      ! a1*( (1+f2)*H_0 - f2*H_minus_dt ) + a2*( (1+f1)*H_0 - f1*H_minus_dt )
      H_aux = d_0*H_0 + d_minus_dt*H_minus_dt
      call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
      ! a1*( (1+f1)*H_0 - f1*H_minus_dt ) + a2*( (1+f2)*H_0 - f2*H_minus_dt )
      H_aux = e_0*H_0 + e_minus_dt*H_minus_dt
      call self%exp_matrix( self%order_Taylor, -zi*self%dt, H_aux, S, x)
    end if
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[RK4]], 
  !> see documentation of the subroutine `RungeKutta4thOrder` for the details.
  !> \(H(-\Delta t)\) is evaluated as
  !> \[H(-\Delta t) = H(0) - [H(\Delta t)-H(0)].\]
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_RK4( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(RK4_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)
    
    CALL_ASSERT( present(H_dt) .neqv. present(H_minus_dt), "only one optional argument must be present" )
    if( present(H_minus_dt) ) then
      call RungeKutta4thOrder( self%dt, zi, H_0, H_minus_dt, S, x )
    else
      call RungeKutta4thOrder( self%dt, zi, H_0, 2*H_0 - H_dt, S, x )
    end if
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[EH]],
  !> see documentation of the subroutine `exphouston_propagator` for the details.
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_EH( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(EH_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)

    call exphouston_propagator( -zi*self%dt, H_0(1:dim, 1:dim), S(1:dim, 1:dim), x(1:dim, :), self%tol, self%n_eigvecs_houston )
  end subroutine

  !> Propagate a single KS wavefunction according to the propagator [[EHM]]
  !> see documentation of the subroutine `exphouston_propagator` for the details.
  !> \(H(\Delta t/2)\) is estimated as
  !> \[ H(\Delta t/2) = \frac{1}{2}\left[ H(\Delta t) + H(0)\right]\]
  !> or as
  !> \[ H(\Delta t/2) = H(0) + \frac{1}{2}\left[ H(0) - H(-\Delta t) \right]\]
  !> The arguments are documented in [[propagate_single_array_]]
  subroutine propagate_single_array_with_EHM( self, dim, H_dt, H_0, H_minus_dt, S, x )
    class(EHM_propagator) :: self
    integer(i32), intent(in) :: dim
    complex(dp), contiguous, optional, intent(in) :: H_dt(:, :)
    complex(dp), contiguous, intent(in) :: H_0(:, :)
    complex(dp), contiguous, optional, intent(in) :: H_minus_dt(:, :)
    complex(dp), contiguous, intent(in) :: S(:, :)
    complex(dp), contiguous, intent(inout) :: x(:, :)

    complex(dp), allocatable :: H_aux(:, :)

    CALL_ASSERT( present(H_dt) .neqv. present(H_minus_dt), "only one optional argument must be present" )
    if( present(H_dt) ) then
      H_aux = 0.5_dp*( H_0(1:dim, 1:dim) + H_dt(1:dim, 1:dim) )
    else
      H_aux = 0.5_dp*( 3*H_0(1:dim, 1:dim) - H_minus_dt(1:dim, 1:dim) )
    end if
    call exphouston_propagator( -zi*self%dt, H_aux, S(1:dim, 1:dim), x(1:dim, :), self%tol, self%n_eigvecs_houston )
  end subroutine

end module