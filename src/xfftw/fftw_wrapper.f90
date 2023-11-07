!> FFTW wrappers for calculating fast fourier transforms (FFT).
!> For more information refer to the documentation of FFTW http://www.fftw.org/.
module xfftw
  use constants, only: zzero, zone, uninit_c_int
  use precision, only: dp 
  use asserts, only: assert
  use modmpi, only: 
  use grid_utils, only: mesh_1d
  use, intrinsic :: ISO_C_binding

  implicit none

#ifdef FFTW3_INTERFACE   
  include "fftw3.f03"
#else
  integer(c_int), parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1, FFTW_ESTIMATE=uninit_c_int
#endif


  private
  public :: fft_type, FFTW_FORWARD, FFTW_BACKWARD, abort_if_not_fftw3
  
  !> Effort to search for optimal FFT plans. When the FFTW interface to MKL DFTI 
  !> is used, this has no effect.
  integer(c_int), parameter :: default_effort = FFTW_ESTIMATE


  !> Type for setting up FFTW plans.
  type fft_type
    !> Dimensions of the FFT.
    integer(c_int), allocatable, private :: dims(:)
    !> Rank.
    integer, private :: rank = uninit_c_int
    !> How many FFTs are planned.
    integer, private :: howmany = uninit_c_int
    !> Sign of the FFT. For forward FFT insert [[FFTW_FORWARD]] (-1), 
    !> for backward FFT insert [[FFTW_BACKWARD]] (1). Default is -1.
    integer(c_int), private :: sign = uninit_c_int
    !> FFTW flag that determines the effort with
    !> which the plan was generated
    integer(c_int), private :: effort = uninit_c_int
    !> Actual FFTW plan.
    type(c_ptr), private :: plan = c_null_ptr
  contains
    generic :: initialize => initialize_noarray, &
                             initialize_inplace
    procedure :: initialize_noarray, &
                 initialize_inplace

    procedure :: delete
    procedure :: execute
  end type fft_type

  contains

  !> Abort exciting if xfftw is used but exciting is not linked FFTW. Call this routine before starting workflows that need fftw.
  subroutine abort_if_not_fftw3(mpi_env, message)
    use modmpi, only: mpiinfo, terminate_mpi_env, mpiglobal
    !> MPI environment to terminate.
    type(mpiinfo), intent(inout), optional :: mpi_env
    !> Message to print to the terminal
    character(*), intent(in), optional :: message

    character(:), allocatable :: message_local
    type(mpiinfo) :: mpi_env_local

    character(*), parameter :: default_message = "Error: Exciting must be linked to FFTW3 to run this routine."
    
    message_local = default_message
    if(present(message)) message_local = message

    mpi_env_local = mpiglobal
    if(present(mpi_env)) mpi_env_local = mpi_env

#ifndef FFTW3_INTERFACE 
    call terminate_mpi_env(mpi_env_local, message=message)
#endif 
  end subroutine

  subroutine initialize_noarray(this, dims, sign)
    class(fft_type), intent(out) :: this
    !> Number of grid points per dimension on which
    !> the function to fourier transform is defined.
    integer, intent(in) :: dims(:)
    !> Sign of the fourier exponential.
    !> <li> forward fft requires -1 </li>
    !> <li> backward fft requires 1 </li>
    !> Default is -1 (forward)
    integer, intent(in) :: sign
    
    integer(c_int), parameter :: effort = default_effort 
    complex(dp), allocatable :: in_out(:)

    this%dims = dims
    this%rank = size(dims)
    this%sign = sign
    this%effort = effort

#ifdef FFTW3_INTERFACE    
    allocate(in_out(product(this%dims)))
    this%plan = fftw_plan_dft( &
      rank=this%rank, &
      n=this%dims(this%rank:1:-1), & ! reverse dims because fftw excepts c_arrays 
      in=in_out, &
      out=in_out, &
      sign=this%sign, &
      flags=this%effort &
    )
    deallocate(in_out)
#endif

  end subroutine initialize_noarray


  !> Initialize fft_type for FFT of a function \( f \).
  subroutine initialize_inplace(this, dims, sign, in_out)
    class(fft_type), intent(out) :: this
    !> Number of grid points per dimension on which
    !> the function to fourier transform is defined.
    integer, intent(in) :: dims(:)
    !> Sign of the fourier exponential.
    !> <li> forward fft requires -1 </li>
    !> <li> backward fft requires 1 </li>
    !> Default is -1 (forward)
    integer, intent(in) :: sign
    !> Array that holds the function \( f \).
    complex(dp), intent(inout) :: in_out(:)

    integer(c_int), parameter :: effort = default_effort 

    call assert(size(in_out) == product(dims), 'size(in_out) /= product(dims).')

    this%dims = dims
    this%rank = size(dims)
    this%sign = sign
    this%effort = effort

#ifdef FFTW3_INTERFACE     
    this%plan = fftw_plan_dft( &
      rank=this%rank, &
      n=this%dims(this%rank:1:-1), & ! reverse dims because fftw excepts c_arrays 
      in=in_out, &
      out=in_out, &
      sign=this%sign, &
      flags=this%effort &
    )
#endif 
  end subroutine initialize_inplace

  !> Execute FFT in place.
  subroutine execute(this, in_out, rescale_forward)
    !> Plan for the fft saved in as fft_type.
    !> If no plan is given, a default plan is generated.
    !> See [init_fft_plancomplex_dp].
    class(fft_type), intent(in) :: this
    !> On input, array holds the function \( f \).
    !> On output, array holds its fourier transform \( \hat{f} \).
    complex(dp), intent(inout) :: in_out(:)
    !> If `.true.`, the result of a forward FFT is rescaled with `1 / product(this%dims).
    !> Default is `.true.`.
    logical, intent(in), optional :: rescale_forward

    logical :: rescale_forward_local
    complex(dp) :: scaling

    rescale_forward_local = .true.
    if(present(rescale_forward)) rescale_forward_local = rescale_forward

#ifdef USE_ASSERT
    call assert(c_associated(this%plan), 'execute: c_associated(this%plan) == .false.')
    call assert(size(in_out) == product(this%dims), 'size(in_out) /= product(this%dims).')
#endif 

#ifdef FFTW3_INTERFACE
    call fftw_execute_dft(this%plan, in_out, in_out)
#endif 

    if (this%sign == FFTW_FORWARD .and. rescale_forward_local) then
      in_out = in_out / size(in_out)
    end if
  end subroutine execute


  !> Delete fft_type plan object
  subroutine delete(this)
    !> fft_type object
    class(fft_type) :: this 

#ifdef FFTW3_INTERFACE
    if (c_associated(this%plan)) call dfftw_destroy_plan(this%plan)
#endif

    if(allocated(this%dims)) deallocate(this%dims)
    this%rank = uninit_c_int
    this%sign = uninit_c_int
    this%effort = uninit_c_int  
  end subroutine delete

end module xfftw
