module rttddft_io_test
  use constants, only: zi
  use exciting_mpi, only: mpiinfo
  use file_utils, only: delete_file
  use math_utils, only: all_close
  use mock_arrays, only: real_matrix_5x7, complex_matrix_5x7
  use modmpi, only: barrier
  use precision, only: dp, i32
  use rttddft_io, only: binary, close_files_vector_fields, delete_jpa_files, delete_pmat_binary_file, &
    delete_pmat_mt_binary_file, delete_wavefunction_binary_file, file_handler, file_pmat_exists, &
    file_pmat_mt_exists, hdf5, open_files_vector_fields, read_vector_field, read_pmat, read_pmat_mt, read_wavefunction, &
    t, t_minus_dt, write_vector_field, write_wavefunction, write_pmat, write_pmat_mt
  use rttddft_CurrentDensity, only: Current_Density_Field
  use rttddft_Polarization, only: Polarization
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: rttddft_io_test_driver

  real(dp), parameter :: tol = 1.0e-10_dp

#ifndef _HDF5_
  character(len=*), parameter :: no_hdf5_message = "Built without support to HDF5. Nothing to test here."
#endif
 
contains

  subroutine rttddft_io_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report

    character(len=*), parameter :: module_tested = 'rttddft_io'

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests
    call test_read_write_vector_field( mpiglobal, test_report )
    call test_read_write_pmat( mpiglobal, test_report )
    call test_read_write_pmat_mt( mpiglobal, test_report )
    call test_read_write_wavefunction( mpiglobal, test_report )

    ! report results
    call test_report%report( module_tested, kill_on_failure )

    ! Finalise test object
    call test_report%finalise()
  end subroutine

  !> Test [[rttddft_io:read_vector_field]] and [[rttddft_io:write_vector_field]]
  subroutine test_read_write_vector_field( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_id = "test_read_write_vector_field"
    character(len=*), parameter :: last = "last ", penult = "penultimate "
    character(len=*), parameter :: str_t = "time", str_j = "j_ind"
    character(len=*), parameter :: str_p = "p_vec", str_ai = "a_ind", str_at = "a_tot"
    character(len=*), parameter :: last_t = last // str_t, last_j = last // str_j
    character(len=*), parameter :: last_p = last // str_p, last_ai = last // str_ai, last_at = last // str_at
    character(len=*), parameter :: penult_j = penult // str_j, penult_p = penult // str_p
    character(len=*), parameter :: penult_ai = penult // str_ai, penult_at = penult // str_at
    real(dp) :: time
    real(dp), allocatable :: times(:), aux(:, :), aux2(:, :)
    type(Current_Density_Field) :: j_ind_t, j_ind_t_minus_dt
    type(Current_Density_Field), allocatable :: j_ind_list(:)
    type(Polarization) :: p_vec_t, p_vec_t_minus_dt
    type(Polarization), allocatable :: p_vec_list(:)
    type(Vector_Potential_Field) :: a_ind_t, a_tot_t, a_ind_t_minus_dt, a_tot_t_minus_dt
    type(Vector_Potential_Field), allocatable :: a_ind_list(:), a_tot_list(:)
    integer(i32) :: n, i, test_counter
    logical :: my_rank_writes

    my_rank_writes = mpiglobal%is_root
    n = 10
    times = [(i*1._dp, i = 1, n)]
    aux = reshape( real_matrix_5x7, [3, n] )
    aux2 = reshape( transpose(real_matrix_5x7), [3, n] )
    allocate( j_ind_list(n), p_vec_list(n), a_ind_list(n), a_tot_list(n) )
    do i = 1, n
      j_ind_list(i)%components = aux(:, i)
      p_vec_list(i)%components = aux2(:, i)
      a_ind_list(i)%components = aux(:, i) - aux2(:, i)**2
      a_tot_list(i)%components = aux2(:, i) + aux(:, i)**2
    end do
    if( my_rank_writes ) then
      call open_files_vector_fields( new=.true. )
      call write_vector_field( times, j_ind_list )
      call write_vector_field( times, p_vec_list )
      call write_vector_field( times, a_ind_list, a_tot_list )
      call close_files_vector_fields
    end if
    call barrier

    test_counter = 1
    call read_vector_field( time, j_ind_t )
    call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter ) )
    associate( j => j_ind_t%components, j_ref => j_ind_list(n)%components )
      call test_report%assert( all_close( j, j_ref, tol ), report_message( test_id, last_j, test_counter) )
      
      test_counter = 2
      call read_vector_field( time, j_ind_t, j_ind_t_minus_dt )
      call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
      call test_report%assert( all_close( j, j_ref, tol ), report_message( test_id, last_j, test_counter) )
    end associate
    associate( j => j_ind_t_minus_dt%components, j_ref => j_ind_list(n-1)%components )
      call test_report%assert( all_close( j, j_ref, tol ), report_message( test_id, penult_j, test_counter) )
    end associate
    
    test_counter = 3
    call read_vector_field( time, p_vec_t )
    call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
    associate( p => p_vec_t%components, p_ref => p_vec_list(n)%components )
      call test_report%assert( all_close( p, p_ref, tol ), report_message( test_id, last_p, test_counter) )

      test_counter = 4
      call read_vector_field( time, p_vec_t, p_vec_t_minus_dt )
      call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
      call test_report%assert( all_close( p, p_ref, tol ), report_message( test_id, last_p, test_counter) )
    end associate
    associate( p => p_vec_t_minus_dt%components, p_ref => p_vec_list(n-1)%components )
      call test_report%assert( all_close( p, p_ref, tol ), report_message( test_id, penult_p, test_counter) )
    end associate

    test_counter = 5
    call read_vector_field( time, a_ind_t, a_tot_t=a_tot_t )
    call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
    associate( ai => a_ind_t%components, ai_ref => a_ind_list(n)%components, at => a_tot_t%components, at_ref => a_tot_list(n)%components )
      call test_report%assert( all_close( ai , ai_ref, tol ), report_message( test_id, last_ai, test_counter) )
      call test_report%assert( all_close( at , at_ref, tol ), report_message( test_id, last_at, test_counter) )

      test_counter = 6
      call read_vector_field( time, a_ind_t, a_ind_t_minus_dt, a_tot_t, a_tot_t_minus_dt )
      call test_report%assert( time == times(n), report_message( test_id, last_t, test_counter) )
      call test_report%assert( all_close( ai , ai_ref, tol ), report_message( test_id, last_ai, test_counter) )
      call test_report%assert( all_close( at , at_ref, tol ), report_message( test_id, last_at, test_counter) )
    end associate
    associate( ai => a_ind_t_minus_dt%components, ai_ref => a_ind_list(n-1)%components, at => a_tot_t_minus_dt%components, at_ref => a_tot_list(n-1)%components )
      call test_report%assert( all_close( ai , ai_ref, tol ), report_message( test_id, penult_ai, test_counter) )
      call test_report%assert( all_close( at , at_ref, tol ), report_message( test_id, penult_at, test_counter) )
    end associate

    call barrier
    if( my_rank_writes ) call delete_jpa_files
  end subroutine

  !> Test [[rttddft_io:read_wavefunction]] and [[rttddft_io:write_wavefunction]]
  subroutine test_read_write_wavefunction( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_id = "test_read_write_wavefunction"
    character(len=*), parameter :: tested_array = "psi"
    integer(i32) :: i, m, n, n_spin
    integer(i32) :: n_kpt, n_kpt_per_proc, first_kpt, last_kpt, test_counter
#ifdef _HDF5_
    integer(i32) :: i_err
#endif
    real(dp), allocatable :: kpt_latt(:, :)
    complex(dp), allocatable :: psi(:, :, :), psi_ref(:, :, :)
    complex(dp), allocatable :: psi_spin(:, :, :, :), psi_spin_ref(:, :, :, :)
    type(file_handler) :: hdf5_handler
    
    n_kpt_per_proc = 2
    n_kpt = n_kpt_per_proc*mpiglobal%procs
    first_kpt = (mpiglobal%rank)*n_kpt_per_proc + 1
    last_kpt = first_kpt + n_kpt_per_proc - 1
    allocate( kpt_latt(3, n_kpt) )
    do i = 1, n_kpt
      kpt_latt(:, i) = [0._dp, 0._dp, real(i, dp)/n_kpt]
    end do

    ! Spin-unpolarized case, format: binary
    test_counter = 1
    m = 3; n = 2; n_spin = 1;
    allocate( psi(m, n, first_kpt:last_kpt), source=reshape( complex_matrix_5x7, [m, n, n_kpt_per_proc] ) )
    psi(1, 1, first_kpt) = psi(1, 1, first_kpt) + (mpiglobal%rank)*zi
    psi_ref = psi
    call write_wavefunction(t, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi, mpi_env=mpiglobal)
    call read_wavefunction(t, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi, mpi_env=mpiglobal)
    call test_report%assert( all_close( psi , psi_ref, tol ), report_message( test_id, tested_array // ' binary format ', test_counter) )
    call barrier()
    if( mpiglobal%is_root ) call delete_wavefunction_binary_file( t )
    call barrier()
    
    ! Spin-unpolarized case, format: HDF5
    test_counter = 2
    hdf5_handler%file_format = hdf5
    hdf5_handler%file_name = "rt.h5"
    hdf5_handler%path = "./"
#ifdef _HDF5_    
    call write_wavefunction( t_minus_dt, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi, mpi_env=mpiglobal, handler=hdf5_handler, n_kpt=n_kpt )
    call read_wavefunction( t_minus_dt, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi, mpi_env=mpiglobal, handler=hdf5_handler )
    call test_report%assert( all_close( psi , psi_ref, tol ), report_message( test_id, tested_array // ' HDF5 format ', test_counter) )
    call barrier
    if( mpiglobal%is_root ) call delete_file( hdf5_handler%file_name, i_err )
    call barrier
#else
    ! fake_message will never be printed - it serves as a hint to developers
    call test_report%assert( .true., no_hdf5_message )
#endif    

    ! Spin-polarized case, format: binary
    test_counter = 3
    m = 3; n = 2; n_spin = 2;
    allocate( psi_spin(m, n, n_spin, first_kpt:last_kpt), source=reshape( complex_matrix_5x7, [m, n, n_spin, n_kpt_per_proc] ) )
    psi_spin(1, 1, 1, first_kpt) = psi_spin(1, 1, 1, first_kpt) + (mpiglobal%rank)*zi
    psi_spin(1, 1, 2, first_kpt) = psi_spin(1, 1, 2, first_kpt)**2 - (mpiglobal%rank)*zi
    psi_spin_ref = psi_spin
    call write_wavefunction( t_minus_dt, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi_spin, mpiglobal )
    call read_wavefunction( t_minus_dt, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi_spin, mpi_env=mpiglobal )
    call test_report%assert( all_close( psi_spin , psi_spin_ref, tol ), report_message( test_id, tested_array // ' binary format - spin', test_counter) )
    call barrier
    if( mpiglobal%is_root ) call delete_wavefunction_binary_file( t_minus_dt )
    call barrier

    ! Spin-polarized case, format: HDF5 (not yet implemented)
    test_counter = 4
#ifdef _HDF5_     
    call write_wavefunction( t, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi_spin, mpiglobal, hdf5_handler, n_kpt )
    call read_wavefunction( t, first_kpt, kpt_latt(:, first_kpt:last_kpt), psi_spin, mpi_env=mpiglobal, handler=hdf5_handler )
    call test_report%assert( all_close( psi_spin , psi_spin_ref, tol ), report_message( test_id, tested_array // ' HDF5 format - spin', test_counter) )
    call barrier
    if( mpiglobal%is_root ) call delete_file( hdf5_handler%file_name, i_err )
    call barrier
#else
    call test_report%assert( .true., no_hdf5_message )
#endif    
  end subroutine

  !> Test [[rttddft_io:read_pmat]] and [[rttddft_io:write_pmat]]
  subroutine test_read_write_pmat( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_id = "test_read_write_pmat"
    character(len=*), parameter :: tested_array = "pmat"
#ifdef _HDF5_
    integer(i32) :: i_err
#endif
    integer(i32) :: n_kpt, n_kpt_per_proc, first_kpt, last_kpt, test_counter
    integer(i32), parameter :: n_cart = 3, m = 3, n = 2
    complex(dp), allocatable :: pmat(:, :, :, :), pmat_ref(:, :, :, :)
    type(file_handler) :: hdf5_handler
    
    n_kpt_per_proc = 2
    n_kpt = n_kpt_per_proc*mpiglobal%procs
    first_kpt = (mpiglobal%rank)*n_kpt_per_proc + 1
    last_kpt = first_kpt + n_kpt_per_proc - 1

    allocate( pmat(m, n, n_cart, first_kpt:last_kpt), &
              source=reshape( [complex_matrix_5x7, complex_matrix_5x7], &
                              [m, n, n_cart, n_kpt_per_proc] ) )
    pmat(1, 1, 1, first_kpt) = pmat(1, 1, 1, first_kpt) + (mpiglobal%rank)*zi
    pmat_ref = pmat
    
    ! Test read/write in binary format
    test_counter = 1
    call write_pmat( first_kpt, pmat(:, :, 1, :), pmat(:, :, 2, :), pmat(:, :, 3, :), mpiglobal )
    call test_report%assert( file_pmat_exists( ), report_message( test_id, tested_array // ' binary file not found', test_counter) )
    call read_pmat( first_kpt, pmat(:, :, 1, :), pmat(:, :, 2, :), pmat(:, :, 3, :), mpiglobal )
    call test_report%assert( all_close( pmat , pmat_ref, tol ), report_message( test_id, tested_array // ' binary format', test_counter) )
    call barrier
    if( mpiglobal%is_root ) call delete_pmat_binary_file()
    
    ! Test read/write in HDF5 format
    test_counter = 2
#ifdef _HDF5_  
    hdf5_handler%file_format = hdf5
    hdf5_handler%file_name = "rt.h5"
    hdf5_handler%path = "./"   
    call write_pmat( first_kpt, pmat_ref(:, :, 1, :), pmat_ref(:, :, 2, :), pmat_ref(:, :, 3, :), &
      mpiglobal, hdf5_handler, n_kpt )
    call test_report%assert( file_pmat_exists( hdf5_handler, mpiglobal ), &
      test_id // ' - ' // tested_array // ' HDF5 file not found, test - ' // to_char(test_counter) )
    call read_pmat( first_kpt, pmat(:, :, 1, :), pmat(:, :, 2, :), pmat(:, :, 3, :), &
      mpiglobal, hdf5_handler )
    call test_report%assert( all_close( pmat , pmat_ref, tol ), report_message( test_id, tested_array // ' HDF5 format', test_counter) )
    call barrier()
    if( mpiglobal%is_root ) call delete_file( hdf5_handler%file_name, i_err )
    call barrier()
#else
    ! Two .true. asserts are needed here
    call test_report%assert( .true., no_hdf5_message )
    call test_report%assert( .true., no_hdf5_message )
#endif  
  end subroutine

  !> Test [[rttddft_io:read_pmat_mt]] and [[rttddft_io:write_pmat_mt]]
  subroutine test_read_write_pmat_mt( mpiglobal, test_report )
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: test_id = "test_read_write_pmat_mt"
    character(len=*), parameter :: tested_array = "pmat_mt"
#ifdef _HDF5_
    integer(i32) :: i_err
#endif    
    integer(i32) :: n_kpt, n_kpt_per_proc, first_kpt, last_kpt, test_counter
    integer(i32), parameter :: n_cart = 3, n_atoms = 4, m = 2, n = 1
    complex(dp), allocatable :: pmat_mt(:, :, :, :, :), pmat_mt_ref(:, :, :, :, :)
    type(file_handler) :: hdf5_handler
    
    n_kpt_per_proc = 2
    n_kpt = n_kpt_per_proc*mpiglobal%procs
    first_kpt = (mpiglobal%rank)*n_kpt_per_proc + 1
    last_kpt = first_kpt + n_kpt_per_proc - 1

    allocate( pmat_mt(m, n, n_cart, n_atoms, first_kpt:last_kpt), &
              source=reshape( [complex_matrix_5x7, complex_matrix_5x7], &
                              [m, n, n_cart, n_atoms, n_kpt_per_proc] ) )
    pmat_mt(1, 1, 1, 1, first_kpt) = pmat_mt(1, 1, 1, 1, first_kpt)**2 + (mpiglobal%rank)*zi
    pmat_mt_ref = pmat_mt

    ! Test read/write in binary format
    test_counter = 1
    call write_pmat_mt( first_kpt, pmat_mt, mpiglobal )
    call test_report%assert( file_pmat_mt_exists( ), report_message( test_id, tested_array // ' binary file not found', test_counter) )
    call read_pmat_mt( first_kpt, pmat_mt, mpiglobal )
    call test_report%assert( all( abs(pmat_mt-pmat_mt_ref) <= tol ), report_message( test_id, tested_array // ' binary format', test_counter) )
    call barrier()
    if( mpiglobal%is_root ) call delete_pmat_mt_binary_file()

    ! Test read/write in HDF5 format
    test_counter = 2
#ifdef _HDF5_  
    hdf5_handler%file_format = hdf5
    hdf5_handler%file_name = "rt.h5"
    hdf5_handler%path = "./"   
    call write_pmat_mt( first_kpt, pmat_mt_ref, mpiglobal, hdf5_handler, n_kpt )
    call test_report%assert( file_pmat_mt_exists( hdf5_handler, mpiglobal ), &
      report_message( test_id, tested_array // ' HDF5 file not found', test_counter) )
    call read_pmat_mt( first_kpt, pmat_mt, mpiglobal, hdf5_handler )
    call test_report%assert( all( abs(pmat_mt-pmat_mt_ref) <= tol ), report_message( test_id, tested_array // ' HDF5 format', test_counter) )
    call barrier()
    if( mpiglobal%is_root ) call delete_file( hdf5_handler%file_name, i_err )
    call barrier()
#else
    call test_report%assert( .true., no_hdf5_message )
#endif  
  end subroutine

  !> (private) generates a report message, given a test case and test number
  function report_message( test_id, test_case, test_number ) result(message)
    !> Name of the subroutine calling `report_message`
    character(len=*), intent(in) :: test_id
    !> Test (sub)identifier, pointing out which case in [[test_id]] is tested
    character(len=*), intent(in) :: test_case
    !> Number for the test
    integer(i32), intent(in) :: test_number
    !> Message to return
    character(len=:), allocatable :: message
    message = test_id // ' - ' // test_case // ' does not match reference: test - ' // to_char(test_number)
  end function
end module