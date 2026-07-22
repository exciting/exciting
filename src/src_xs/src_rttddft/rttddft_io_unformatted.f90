module rttddft_io_unformatted
#include "asserts.fpp"
  use exciting_mpi, only: mpiinfo
  use file_utils, only: add_default_extension, delete_file
  use MD, only: trajectory
  use MD_io, only: read_trajectory, write_trajectory
  use mod_corestate, only: rhocr
  use mod_misc, only: filext
  use mod_potential_and_density, only: vxcir, vxcmt, rhoir, rhomt, vclir, vclmt, veffir, veffmt, veffig, meffig
  use os_utils, only: path_exists
  use precision, only: dp, i32
  use rttddft_file_formats, only: file_handler, restart_format, binary, hdf5
  use rttddft_file_names, only: filename_phases, filename_pmat, filename_pmat_mt, &
    filename_rho_vks, filename_wavefunction, filename_wavefunction_second_variation, &
    kpt_latt_name, RTTDDFT_suffix, &
    suffix_wavefunction_gnd, suffix_wavefunction_t, suffix_wavefunction_t_minus_dt
  use rttddft_io_hdf5, only: dataset_exists, read_array_hdf5, read_three_arrays_hdf5, write_array_hdf5, write_three_arrays_hdf5
#ifdef MPI
  use rttddft_io_parallel, only: read_array, read_three_arrays, write_array, write_three_arrays
#else
  use rttddft_io_serial, only: read_array, read_three_arrays, write_array, write_three_arrays
#endif
  use to_char_conversion, only: to_char

  implicit none
  
  private

  ! Methods
  public :: delete_pmat_binary_file, delete_pmat_mt_binary_file, delete_wavefunction_binary_file, &
    file_pmat_exists, file_pmat_mt_exists, &
    get_filename_pmat, get_filename_pmat_mt, &
    read_phases, read_pmat, read_pmat_mt, read_state_Ehrenfest_MD, read_wavefunction, &
    write_phases, write_pmat, write_pmat_mt, write_state_Ehrenfest_MD, write_wavefunction

  ! Variables
  public :: groundstate, t, t_minus_dt

  interface read_wavefunction
    module procedure :: read_wavefunction_non_spin_spiral
    module procedure :: read_wavefunction_spin_polarized
  end interface

  interface write_wavefunction
    module procedure :: write_wavefunction_non_spin_polarized
    module procedure :: write_wavefunction_spin_polarized
  end interface

  enum, bind(C)
    enumerator :: wavefunction_case
    enumerator :: groundstate, t, t_minus_dt
  end enum

contains 
  !> Check if a binary file or HDF5 dataset exists
  logical function rttddft_file_exists( name, handler, mpi_env ) result( exists )
    !> File name or dataset name
    character(len=*), intent(in) :: name
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), optional, intent(in) :: mpi_env

    integer(kind(restart_format)) :: file_format
    
    CALL_ASSERT( present( handler ) .eqv. present( mpi_env ), "mpi_env must be passed when handler is present" )

    file_format = binary
    if( present( handler ) ) file_format = handler%file_format
    select case(file_format)
      case( binary )
        exists = path_exists( trim( name ) )
      case( hdf5 )
        call handler%assert_consistency( )
        exists = dataset_exists( handler%file_name, handler%path, trim( name ), mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
        exists = .false.
    end select
  end function

  logical function file_pmat_exists( handler, mpi_env ) result(ok)
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), optional, intent(in) :: mpi_env

    logical :: special

    special = .false.
    if( present( handler ) ) then
      special = ( handler%file_format == hdf5 )
    end if
    if( .not. special ) then
      ok = rttddft_file_exists( get_filename_pmat( ), handler, mpi_env )
    else
      ok = rttddft_file_exists( get_filename_pmat( )//"-1", handler, mpi_env ) &
        .and. rttddft_file_exists( get_filename_pmat( )//"-2", handler, mpi_env ) &
        .and. rttddft_file_exists( get_filename_pmat( )//"-3", handler, mpi_env )
    end if
  end function

  function get_filename_pmat() result(name)
    character(len=:), allocatable :: name
    name = add_default_extension( filename_pmat )
  end function

  !> Read the momentum matrix elements from file
  subroutine read_pmat( first_kpt, pmat_x, pmat_y, pmat_z, mpi_env, handler )
    !> Index of the first `k-point` to be considered
    integer,intent(in) :: first_kpt
    !> Momentum matrix elements - x component
    complex(dp), contiguous, intent(inout)   :: pmat_x(:, :, first_kpt:)
    !> Momentum matrix elements - y component
    complex(dp), contiguous, intent(inout)   :: pmat_y(:, :, first_kpt:)
    !> Momentum matrix elements - z component
    complex(dp), contiguous, intent(inout)   :: pmat_z(:, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    integer(kind(restart_format)) :: file_format
    
    file_format = binary
    if( present( handler ) ) file_format = handler%file_format
    select case(file_format)
      case( binary )
        call read_three_arrays( get_filename_pmat( ), first_kpt, pmat_x, pmat_y, pmat_z, mpi_env=mpi_env )
      case( hdf5 )
        call handler%assert_consistency( )
        call read_three_arrays_hdf5( handler%file_name, handler%path, get_filename_pmat( ), first_kpt, &
          pmat_x, pmat_y, pmat_z, mpi_env=mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Write the momentum matrix elements to file
  subroutine write_pmat( first_kpt, pmat_x, pmat_y, pmat_z, mpi_env, handler, n_kpt )
    !> Index of the first `k-point` to be considered in the sum
    integer, intent(in) :: first_kpt
    !> Momentum matrix elements - x component
    complex(dp), contiguous, intent(inout)   :: pmat_x(:, :, first_kpt:)
    !> Momentum matrix elements - y component
    complex(dp), contiguous, intent(inout)   :: pmat_y(:, :, first_kpt:)
    !> Momentum matrix elements - z component
    complex(dp), contiguous, intent(inout)   :: pmat_z(:, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> Number of k-points
    integer(i32), optional, intent(in) :: n_kpt
    
    integer(kind(restart_format)) :: file_format

    CALL_ASSERT( present(n_kpt) .eqv. present(handler), "n_kpt must be passed when handler is present")

    file_format = binary
    if( present( handler ) ) file_format = handler%file_format
    select case(file_format)
      case( binary )
        call write_three_arrays( get_filename_pmat( ), first_kpt, pmat_x, pmat_y, pmat_z, mpi_env=mpi_env )
      case( hdf5 )
        call handler%assert_consistency( )
        call write_three_arrays_hdf5( handler%file_name, handler%path, get_filename_pmat( ), &
          pmat_x, pmat_y, pmat_z, first_kpt, n_kpt, mpi_env=mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Delete the file containing `pmat` in binary format
  subroutine delete_pmat_binary_file( )
    integer(i32) :: i_error
    call delete_file( get_filename_pmat( ), i_error )
  end subroutine

  !> Check if file with `pmat_mt` exists
  logical function file_pmat_mt_exists( handler, mpi_env )
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), optional, intent(in) :: mpi_env

    file_pmat_mt_exists = rttddft_file_exists( get_filename_pmat_mt( ), handler, mpi_env )
  end function

  function get_filename_pmat_mt() result(name)
    character(len=:), allocatable :: name
    name = add_default_extension( filename_pmat_mt )
  end function

  !> Read the muffin-tin part of the momentum matrix (`pmat_mt`) from file
  subroutine read_pmat_mt( first_kpt, pmat_mt, mpi_env, handler )
    !> Index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> Muffin-tin part of the momentum matrix
    complex(dp), intent(out)  :: pmat_mt(:, :, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    integer(kind(restart_format)) :: file_format
    
    file_format = binary
    if( present( handler ) ) file_format = handler%file_format
    select case(file_format)
      case( binary )
        call read_array( get_filename_pmat_mt( ), first_kpt, pmat_mt, mpi_env )
      case( hdf5 )
        call handler%assert_consistency( )
        call read_array_hdf5( handler%file_name, handler%path, get_filename_pmat_mt( ), first_kpt, pmat_mt, mpi_env=mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Write the muffin-tin part of the momentum matrix (`pmat_mt`) to file
  subroutine write_pmat_mt( first_kpt, pmat_mt, mpi_env, handler, n_kpt )
    !> Index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> Muffin-tin part of the momentum matrix
    complex(dp), intent(in)   :: pmat_mt(:, :, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> Number of k-points
    integer(i32), optional, intent(in) :: n_kpt
    
    integer(kind(restart_format)) :: file_format

    CALL_ASSERT( present(n_kpt) .eqv. present(handler), "n_kpt must be passed when handler is present")

    file_format = binary
    if( present( handler ) ) file_format = handler%file_format
    select case(file_format)
      case( binary )
        call write_array( get_filename_pmat_mt( ), first_kpt, pmat_mt, mpi_env )
      case( hdf5 )
        call handler%assert_consistency( )
        call write_array_hdf5( handler%file_name, handler%path, get_filename_pmat_mt( ), pmat_mt, first_kpt, n_kpt, mpi_env=mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Delete the file containing `pmat` in binary format
  subroutine delete_pmat_mt_binary_file( )
    integer(i32) :: i_error
    call delete_file( get_filename_pmat_mt( ), i_error )
  end subroutine

  !> (Private) Return `suffix_wavefunction_t`, `suffix_wavefunction_t_minus_dt` or `suffix_wavefunction_gnd`
  pure function get_suffix_filename_wavefunction( psi_case ) result(suffix)
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)), intent(in) :: psi_case
    !> File name (to return)
    character(len=:), allocatable :: suffix

    select case( psi_case )
      case ( t )
        suffix = suffix_wavefunction_t
      case ( t_minus_dt)
        suffix = suffix_wavefunction_t_minus_dt
      case ( groundstate )
        suffix = suffix_wavefunction_gnd
      case default
        suffix = suffix_wavefunction_t
    end select
  end function

  !> (Private) Return `filename_wavefunction_previous`, `filename_wavefunction` or `filename_wavefunction_gnd`
  pure function get_filename_wavefunction( psi_case ) result(name)
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)), intent(in) :: psi_case
    !> File name (to return)
    character(len=:), allocatable :: name

    name = add_default_extension( filename_wavefunction // get_suffix_filename_wavefunction(psi_case) )
  end function

  !> (Private) Return `filename_wavefunction_previous`, `filename_wavefunction` or `filename_wavefunction_gnd`
  pure function get_filename_wavefunction_second_variation( psi_case ) result(name)
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)), intent(in) :: psi_case
    !> File name (to return)
    character(len=:), allocatable :: name

    name = add_default_extension( filename_wavefunction_second_variation // get_suffix_filename_wavefunction(psi_case) )
  end function

  !> Read wavefunction coefficients from file. Similar to [[getevecfv]], but 
  !> does not need to split files and can be used by multiple MPI procs simultaneously.
  subroutine read_wavefunction_non_spin_spiral( psi_case, first_kpt, kpt_latt, psi, psi_sv, mpi_env, handler )
    !> Enum telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)) :: psi_case
    !> First k-point treated by this (MPI)rank
    integer(i32), intent(in) :: first_kpt
    !> k-points in lattice coordinates
    real(dp), intent(in) :: kpt_latt(:, first_kpt:)
    !> Basis-expansion coefficients of the (not spin-spiral) KS-WFs
    complex(dp), contiguous, target, intent(out) :: psi(:, :, first_kpt:)
    !> Basis-expansion coefficients of the second-variational KS-WFs
    complex(dp), contiguous, optional, intent(out) :: psi_sv(:, :, first_kpt:)
    !> MPI environment (needed to read in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    integer(i32), parameter :: n_spin = 1
    complex(dp), contiguous, pointer :: ptr(:, :, :, :)

    ! Map wavefunction to spin-polarized wavefunction
    associate( m => size(psi, 1), n => size(psi, 2), last_kpt => ubound( psi, 3 ) )
      ptr(1:m, 1:n, 1:n_spin, first_kpt:last_kpt) => psi
    end associate
    call read_wavefunction_spin_polarized( psi_case, first_kpt, kpt_latt, ptr, psi_sv, &
      mpi_env, handler )
  end subroutine

  !> Same as [[read_wavefunction_non_spin_polarized]], but for the spin polarized case
  subroutine read_wavefunction_spin_polarized( psi_case, first_kpt, kpt_latt, psi, psi_sv, mpi_env, handler )
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)) :: psi_case
    !> First k-point treated by this (MPI)rank
    integer(i32), intent(in) :: first_kpt
    !> k-points in lattice coordinates
    real(dp), contiguous, intent(in) :: kpt_latt(:, first_kpt:)
    !> Basis-expansion coefficients of the (spin-polarized) KS-WFs
    complex(dp), contiguous, intent(out) :: psi(:, :, :, first_kpt:)
    !> Basis-expansion coefficients of the second-variational KS-WFs
    complex(dp), contiguous, optional, intent(out) :: psi_sv(:, :, first_kpt:)
    !> MPI environment (needed to read in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    integer(i32), parameter :: n_cartesian_coords = 3, n_spin_max = 2
    integer(kind(restart_format)) :: file_format

    file_format = binary
    if( present(handler) ) file_format = handler%file_format
    associate( n_spin => size(psi, 3) )
      CALL_ASSERT( n_spin <= n_spin_max, "psi has more spin polarizations than allowed")
      CALL_ASSERT( size(kpt_latt, 1) == n_cartesian_coords, to_char(n_cartesian_coords) // " cartesian components are expected" )
      CALL_ASSERT( size(kpt_latt, 2) == size(psi, 4), "kpt_latt and psi must be compatible.")
    end associate
    select case(file_format)
      case( binary )
        call read_array( get_filename_wavefunction(psi_case), first_kpt, psi, kpt_latt, mpi_env )
        if( present(psi_sv) ) call read_array( get_filename_wavefunction_second_variation(psi_case), first_kpt, psi_sv, kpt_latt, mpi_env )
      case( hdf5 )
        call handler%assert_consistency( )
        call read_array_hdf5( handler%file_name, handler%path, get_filename_wavefunction(psi_case), first_kpt, psi, kpt_latt, kpt_latt_name, mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine
  
  !> Write the wavefunction coefficients to file. Similar to [[putevecfv]], but 
  !> does not need to split files and can be used by multiple MPI procs simultaneously.
  subroutine write_wavefunction_non_spin_polarized( psi_case, first_kpt, kpt_latt, psi, mpi_env, handler, n_kpt )
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)) :: psi_case
    !> index of the first `k-point` to be considered in the sum
    integer(i32), intent(in) :: first_kpt
    !> k-points in lattice coordinates
    real(dp), contiguous, intent(inout) :: kpt_latt(:, first_kpt:)
    !> wavefunction coefficients
    complex(dp), contiguous, target, intent(inout) :: psi(:, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> Number of k-points
    integer(i32), optional, intent(in) :: n_kpt
    
    complex(dp), contiguous, pointer :: ptr(:, :, :, :)
    integer(i32), parameter :: n_spin = 1
    
    ! Map wavefunction to spin-polarized wavefunction
    associate( m => size(psi, 1), n => size(psi, 2), last_kpt => ubound( psi, 3 ) )
      ptr(1:m, 1:n, 1:n_spin, first_kpt:last_kpt) => psi
    end associate
    CALL_ASSERT( present(n_kpt) .eqv. present(handler), "n_kpt must be passed when handler is present")
    call write_wavefunction_spin_polarized( psi_case, first_kpt, kpt_latt, ptr, mpi_env, handler, n_kpt )
  end subroutine

  !> Same as [[write_wavefunction_non_spin_polarized]], but for the spin polarized case
  subroutine write_wavefunction_spin_polarized( psi_case, first_kpt, kpt_latt, psi, mpi_env, handler, n_kpt )
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)) :: psi_case
    !> index of the first `k-point` to be considered in the sum
    integer(i32), intent(in) :: first_kpt
    !> k-points in lattice coordinates
    real(dp), contiguous, intent(inout) :: kpt_latt(:, first_kpt:)
    !> wavefunction coefficients
    complex(dp), contiguous, intent(inout) :: psi(:, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> Number of k-points
    integer(i32), optional, intent(in) :: n_kpt

    integer(i32), parameter :: n_spin_max = 2, n_cartesian_coords = 3
    integer(kind(restart_format)) :: file_format

    CALL_ASSERT( present(n_kpt) .eqv. present(handler), "n_kpt must be passed when handler is present")
    file_format = binary
    if( present(handler) ) file_format = handler%file_format
    associate( n_spin => size(psi, 3) )
      CALL_ASSERT( n_spin <= n_spin_max, "psi has more spin polarizations than allowed")
      CALL_ASSERT( size(kpt_latt, 1) == n_cartesian_coords, "kpt_latt must have size 3 along 1st dim.")
      CALL_ASSERT( size(kpt_latt, 2) == size(psi, 4), "kpt_latt and psi must be compatible.")
      select case(file_format)
        case( binary )
          call write_array( get_filename_wavefunction(psi_case), first_kpt, psi, kpt_latt, mpi_env=mpi_env )
        case( hdf5 )
          call handler%assert_consistency( )
          call write_array_hdf5( handler%file_name, handler%path, get_filename_wavefunction(psi_case), &
            psi, first_kpt, n_kpt, kpt_latt, kpt_latt_name, mpi_env )
        case default
          CALL_ASSERT( .false., "Unrecognized format" )
      end select
    end associate
  end subroutine

  !> Delete the file in binary format which contains the KS wavefunctions
  subroutine delete_wavefunction_binary_file( psi_case )
    !> Enum containing telling if `psi` refers to \(t\), \(t-\Delta t\), or to groundstate
    integer(kind(wavefunction_case)) :: psi_case
    integer(i32) :: i_error
    call delete_file( get_filename_wavefunction( psi_case ), i_error )
  end subroutine

  !> Write current phases
  subroutine write_phases( phases_to_match, handler, mpi_env )
    !> Phases needed for matching in the MTP polarization calculation
    real(dp), intent(in) :: phases_to_match(:, :)
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: unit
    
    select case( handler%file_format )
      case( binary )
        if( mpi_env%is_root ) then
          open( newunit=unit, file=add_default_extension( filename_phases ), action="write", &
            form="unformatted", access="stream" )
          write( unit ) phases_to_match
          close( unit )
        end if
      case( hdf5 )
        call handler%assert_consistency( )
        call write_array_hdf5( handler%file_name, handler%path, add_default_extension( filename_phases ), &
          phases_to_match, mpi_env, serial_access=.true. )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Read current phases
  subroutine read_phases( phases_to_match, handler, mpi_env )
    !> Phases needed for matching in the MTP polarization calculation
    real(dp), intent(out) :: phases_to_match(:, :)
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: unit

    select case( handler%file_format )
      case( binary )
        open( newunit=unit, file=add_default_extension( filename_phases ), action="read", &
          form="unformatted", access="stream" )
        read( unit ) phases_to_match
        close( unit )
      case( hdf5 )
        call handler%assert_consistency( )
        call read_array_hdf5( handler%file_name, handler%path, add_default_extension( filename_phases ), &
          phases_to_match, mpi_env )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Write current state
  subroutine write_state_Ehrenfest_MD( nuclei_motion, handler, mpi_env )
    !> This argument packs nuclei positions and velocities
    class(trajectory), intent(inout) :: nuclei_motion
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    character(len=:), allocatable :: string

    string = trim( filext )
    filext = RTTDDFT_suffix // trim( filext )
    call write_rho_vks( handler, mpi_env )
    call write_trajectory( nuclei_motion, handler, mpi_env )
    filext = string
  end subroutine

  !> Read current state
  subroutine read_state_Ehrenfest_MD( nuclei_motion, handler, mpi_env )
    !> This argument packs nuclei positions and velocities
    class(trajectory), intent(inout) :: nuclei_motion
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    character(len=:), allocatable :: string

    string = trim( filext )
    filext = RTTDDFT_suffix // trim( filext )
    call read_rho_vks( handler, mpi_env )
    call read_trajectory( nuclei_motion, handler, mpi_env )
    filext = string
  end subroutine

  !> Write the core and valence densities and the KS potential
  subroutine write_rho_vks( handler, mpi_env )
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: unit

    select case( handler%file_format )
      case( binary )
        if( mpi_env%is_root ) then
          open( newunit=unit, file=add_default_extension( filename_rho_vks ), action="write", &
            form="unformatted", access="stream" )
          write( unit ) rhocr, rhomt, rhoir
          write( unit ) vclmt, vclir, vxcmt, vxcir, veffmt, veffir, veffig, meffig
          close( unit )
        end if
      case( hdf5 )
        call handler%assert_consistency( )
        call wrapper_write_real_array_hdf5( 'rhocr', rhocr )
        call wrapper_write_real_array_hdf5( 'rhomt', rhomt )
        call wrapper_write_real_array_hdf5( 'rhoir', rhoir )
        call wrapper_write_real_array_hdf5( 'vclmt', vclmt )
        call wrapper_write_real_array_hdf5( 'vclir', vclir )
        call wrapper_write_real_array_hdf5( 'vxcmt', vxcmt )
        call wrapper_write_real_array_hdf5( 'vxcir', vxcir )
        call wrapper_write_real_array_hdf5( 'veffmt', veffmt )
        call wrapper_write_real_array_hdf5( 'veffir', veffir )
        call wrapper_write_complex_array_hdf5( 'veffig', veffig )
        call wrapper_write_complex_array_hdf5( 'meffig', meffig )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
    contains 
      !> Wrapper to call [[write_array_hdf5]] for complex arrays
      subroutine wrapper_write_complex_array_hdf5( name, array )
        !> Dataset name
        character(len=*), intent(in) :: name
        !> Array to write
        complex(dp), contiguous, intent(in) :: array(:)
        call write_array_hdf5( handler%file_name, handler%path, filename_rho_vks // '-' // add_default_extension(name), &
          array, mpi_env, serial_access=.true. )
      end subroutine
      !> Same as [[wrapper_write_complex_array_hdf5]] for real arrays
      subroutine wrapper_write_real_array_hdf5( name, array )
        !> Dataset name
        character(len=*), intent(in) :: name
        !> Array to write
        real(dp), contiguous, intent(in) :: array(..)
        call write_array_hdf5( handler%file_name, handler%path, filename_rho_vks // '-' // add_default_extension(name), &
          array, mpi_env, serial_access=.true. )
      end subroutine
  end subroutine

  !> Read the core and valence densities and the KS potential
  subroutine read_rho_vks( handler, mpi_env )
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: unit

    select case( handler%file_format )
      case( binary )
        open( newunit=unit, file=add_default_extension( filename_rho_vks ), action="read", form="unformatted", access="stream" )
        read( unit ) rhocr, rhomt, rhoir
        read( unit ) vclmt, vclir, vxcmt, vxcir, veffmt, veffir, veffig, meffig
        close( unit )
      case( hdf5 )
        call handler%assert_consistency( )
        call wrapper_read_real_array_hdf5( 'rhocr', rhocr )
        call wrapper_read_real_array_hdf5( 'rhomt', rhomt )
        call wrapper_read_real_array_hdf5( 'rhoir', rhoir )
        call wrapper_read_real_array_hdf5( 'vclmt', vclmt )
        call wrapper_read_real_array_hdf5( 'vclir', vclir )
        call wrapper_read_real_array_hdf5( 'vxcmt', vxcmt )
        call wrapper_read_real_array_hdf5( 'vxcir', vxcir )
        call wrapper_read_real_array_hdf5( 'veffmt', veffmt )
        call wrapper_read_real_array_hdf5( 'veffir', veffir )
        call wrapper_read_complex_array_hdf5( 'veffig', veffig )
        call wrapper_read_complex_array_hdf5( 'meffig', meffig )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
    contains 
      !> Wrapper to call [[read_array_hdf5]] for complex arrays
      subroutine wrapper_read_complex_array_hdf5( name, array )
        !> Dataset name
        character(len=*), intent(in) :: name
        !> Array to write
        complex(dp), contiguous, intent(out) :: array(:)
        call read_array_hdf5( handler%file_name, handler%path, filename_rho_vks // "-" // add_default_extension(name), &
          array, mpi_env )
      end subroutine
      !> Same as [[wrapper_read_complex_array_hdf5]] for real arrays
      subroutine wrapper_read_real_array_hdf5( name, array )
        !> Dataset name
        character(len=*), intent(in) :: name
        !> Array to write
        real(dp), contiguous, intent(out) :: array(..)
        call read_array_hdf5( handler%file_name, handler%path, filename_rho_vks // "-" // add_default_extension(name), &
          array, mpi_env )
      end subroutine
  end subroutine
end module