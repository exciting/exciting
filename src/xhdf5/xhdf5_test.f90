module xhdf5_test
  use precision, only: sp, dp
  use modmpi, only: mpiinfo, barrier, nofset, firstofset, lastofset
  use unit_test_framework, only : unit_test_type
  use xhdf5
  use os_utils, only: system_cmd
  use math_utils, only: mod1, all_close

  implicit none

  private
  public :: xhdf5_test_driver    
  

  contains


  !> Run tests for math tools
  subroutine xhdf5_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure  
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 47

#ifdef _HDF5_
    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    
    call test_xhdf5_file(test_report, mpiglobal)
    call barrier(mpiglobal)

    ! CHARACTER

    call test_xhdf5_write_and_read_character(test_report, mpiglobal)
    call barrier(mpiglobal)

    ! INTEGER

    call test_xhdf5_write_and_read_integer_rank_0(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_integer_rank_1(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_integer_rank_2(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_integer_rank_3(test_report, mpiglobal)
    call barrier(mpiglobal)

    ! REAL(SP)

    call test_xhdf5_write_and_read_real_sp_rank_0(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_sp_rank_1(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_sp_rank_2(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_sp_rank_3(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_sp_rank_4(test_report, mpiglobal)
    call barrier(mpiglobal)

    ! REAL(DP)

    call test_xhdf5_write_and_read_real_dp_rank_0(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_dp_rank_1(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_dp_rank_2(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_dp_rank_3(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_real_dp_rank_4(test_report, mpiglobal)
    call barrier(mpiglobal)

    ! COMPLEX(DP)

    call test_xhdf5_write_and_read_cmplx_dp_rank_1(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_cmplx_dp_rank_2(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_cmplx_dp_rank_3(test_report, mpiglobal)
    call barrier(mpiglobal)

    call test_xhdf5_write_and_read_cmplx_dp_rank_4(test_report, mpiglobal)
    call barrier(mpiglobal)
#else 
      print*, 'Built without HDF5. Nothing to do here.'
      call test_report%init(0, mpiglobal)
#endif

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('xhdf5', kill_on_failure)
    else
      call test_report%report('xhdf5')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine xhdf5_test_driver


  !> Test hdf5 file utilities
  subroutine test_xhdf5_file(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr

    real(dp), parameter :: real_vec(5) = [1._dp, -23.0123_dp, 0._dp, 123.3245_dp, 3142._dp]
    complex(dp), parameter :: cmplx_vec(5) = cmplx(real_vec, -real_vec, dp)
    integer, allocatable :: dset_shape(:), dset_shape_ref(:)

    ! Create a new HDF5 file
    call h5%initialize('test_file.h5', mpiglobal)
    call h5%finalize()

    ! Open an existing HDF5 file
    call h5%initialize('test_file.h5', mpiglobal)

    ! Create a new group
    call h5%initialize_group('./', 'new_group')
    call h5%initialize_group('./', 'new_group') ! Can initialize_group handle an existing group?

    ! Create a sub group
    call h5%initialize_group('new_group', 'sub_group')
    call h5%initialize_group('new_group', 'sub_group') ! Can initialize_group handle an existing group?

    ! Test link_exists
    call test_report%assert(h5%exists('new_group'), &
                           "Failed to create 'new_group'.")

    call h5%delete('new_group')
    call test_report%assert(.not. h5%exists('new_group'), &
                                 "Failed to delete 'new_group'.") 

    call h5%write('/', 'real_vec', real_vec, [1], [5])
    call h5%write('/', 'cmplx_vec', cmplx_vec, [1], [5])

    call h5%dataset_shape('/', 'real_vec', dset_shape)
    call test_report%assert(all(dset_shape == [5]), &
            'dataset shape for real vector not correctly detected.')

    deallocate(dset_shape)
    call h5%dataset_shape('/', 'cmplx_vec', dset_shape)
    call test_report%assert(all(dset_shape == [2, 5]), &
            'dataset shape for complex vector not correctly detected (without complex specification).')

    deallocate(dset_shape)
    call h5%dataset_shape('/', 'cmplx_vec', dset_shape, complex_dataset=.true.)
    call test_report%assert(all(dset_shape == [5]), &
            'dataset shape for complex vector not correctly detected (with complex specification).')

    deallocate(dset_shape)
    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_file


! CHARACTER


  !> Test read and write character.
  subroutine test_xhdf5_write_and_read_character(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr

    character(:), allocatable :: char_write, char_read
    logical :: bool_write, bool_read

    char_write = 'aloifhSADpi$%^flmdsa'

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')

    call h5%write('datasets', 'character', char_write)
    call h5%read('datasets', 'character', char_read)
    call test_report%assert(char_write .eq. char_read, &
            'character: Read array differs from written array.')

    bool_write = .true.
    call h5%write('datasets', 'bool_true', bool_write)
    call h5%read('datasets', 'bool_true', bool_read)
    call test_report%assert(bool_read .eqv. bool_write, &
            'bool: Read and write .true. failed.')

    bool_write = .false.
    call h5%write('datasets', 'bool_true', bool_write)
    call h5%read('datasets', 'bool_true', bool_read)
    call test_report%assert(bool_read .eqv. bool_write, &
            'bool: Read and write .false. failed.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_character


! INTEGER


  !> Test read and write for integer scalars.
  subroutine test_xhdf5_write_and_read_integer_rank_0(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr

    integer :: integer_write, integer_read

    integer_write = 15

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')

    call h5%write('datasets', 'integer_r1', integer_write)
    call h5%read('datasets', 'integer_r1', integer_read)
    call test_report%assert(integer_write == integer_read, &
            'integer rank 0: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_integer_rank_0


  !> Test read and write for integer rank 1 arrays.
  subroutine test_xhdf5_write_and_read_integer_rank_1(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, l

    integer :: integer_r1(8)
    integer, allocatable :: datachunk(:), datachunk_read(:)
    integer :: offset(1), dataset_shape(1), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(integer_r1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(integer_r1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(integer_r1), mpiglobal%procs)

    do l = 1, size(integer_r1)

      integer_r1(l) =  l

    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    datachunk = integer_r1(first : last)
    datachunk_read = spread(0, 1, size(datachunk))
    offset = first
    dataset_shape = shape(integer_r1)

    call h5%write('datasets', 'integer_r1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'integer_r1', datachunk_read, offset)
    call test_report%assert(all(datachunk == datachunk_read), &
            'integer rank 1: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_integer_rank_1


  !> Test read and write for integer rank 2 arrays.
  subroutine test_xhdf5_write_and_read_integer_rank_2(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    integer :: integer_r2(8, 3)
    integer, allocatable :: datachunk(:, :), datachunk_read(:, :)
    integer :: offset(2), dataset_shape(2), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(integer_r2, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(integer_r2, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(integer_r2, 2), mpiglobal%procs)

    do l = 1, size(integer_r2, 2)
      do k = 1, size(integer_r2, 1)

        integer_r2(k, l) =  (l - k)

      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over second rank
    datachunk = integer_r2(:, first : last)
    datachunk_read = datachunk
    datachunk_read = 0.
    offset = [1, first]
    dataset_shape = shape(integer_r2)

    call h5%write('datasets', 'integer_r2-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'integer_r2-2', datachunk_read, offset)
    call test_report%assert(all(datachunk == datachunk_read), &
            'integer rank 2, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(integer_r2, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(integer_r2, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(integer_r2, 1), mpiglobal%procs)

    datachunk = integer_r2(first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0.
    offset = [first, 1]
    dataset_shape = shape(integer_r2)

    call h5%write('datasets', 'integer_r2-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'integer_r2-1', datachunk_read, offset)
    call test_report%assert(all(datachunk == datachunk_read), &
            'integer rank 2, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_integer_rank_2


  !> Test read and write for integer rank 3 arrays.
  subroutine test_xhdf5_write_and_read_integer_rank_3(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(sp) :: integer_r3(8, 3, 5)
    real(sp), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
    integer :: offset(3), dataset_shape(3), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(integer_r3, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(integer_r3, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(integer_r3, 3), mpiglobal%procs)

    do l = 1, size(integer_r3, 3)
      do k = 1, size(integer_r3, 2)
        do j = 1, size(integer_r3, 1)

          integer_r3(j, k, l) =  (j - k) + l

        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over third rank
    datachunk = integer_r3(:, :, first : last)
    datachunk_read = datachunk
    datachunk_read = 0
    offset = [1, 1, first]
    dataset_shape = shape(integer_r3)

    call h5%write('datasets', 'integer_r3-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'integer_r3-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 3, dim 3 distributed: Read array differs from written array.')


    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(integer_r3, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(integer_r3, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(integer_r3, 2), mpiglobal%procs)

    datachunk = integer_r3(:, first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0
    offset = [1, first, 1]
    dataset_shape = shape(integer_r3)

    call h5%write('datasets', 'integer_r3-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'integer_r3-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 3, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(integer_r3, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(integer_r3, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(integer_r3, 1), mpiglobal%procs)        
            
    datachunk = integer_r3(first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = 0
    offset = [first, 1, 1]
    dataset_shape = shape(integer_r3)

    call h5%write('datasets', 'integer_r3-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'integer_r3-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 3, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_integer_rank_3


  ! REAL(SP)


  !> Test read and write for real(sp) scalars.
  subroutine test_xhdf5_write_and_read_real_sp_rank_0(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(sp) :: real_write, real_read

    real_write = 0.213_sp

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')

    call h5%write('datasets', 'real_r1', real_write)
    call h5%read('datasets', 'real_r1', real_read)
    call test_report%assert(all_close(real_write, real_read), &
            'real(sp) rank 0: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_sp_rank_0


  !> Test read and write for real(sp) rank 1 arrays.
  subroutine test_xhdf5_write_and_read_real_sp_rank_1(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(sp) :: real_r1(8)
    real(sp), allocatable :: datachunk(:), datachunk_read(:)
    integer :: offset(1), dataset_shape(1), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)

    do l = 1, size(real_r1)

      real_r1(l) =  real(l, sp)

    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    datachunk = real_r1(first : last)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = first
    dataset_shape = shape(real_r1)

    call h5%write('datasets', 'real_r1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 1: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_sp_rank_1


  !> Test read and write for real(sp) rank 2 arrays.
  subroutine test_xhdf5_write_and_read_real_sp_rank_2(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(sp) :: real_r2(8, 3)
    real(sp), allocatable :: datachunk(:, :), datachunk_read(:, :)
    integer :: offset(2), dataset_shape(2), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)

    do l = 1, size(real_r2, 2)
      do k = 1, size(real_r2, 1)

        real_r2(k, l) =  real((l - k), sp)

      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over second rank
    datachunk = real_r2(:, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = [1, first]
    dataset_shape = shape(real_r2)

    call h5%write('datasets', 'real_r2-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r2-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 2, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)

    datachunk = real_r2(first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = [first, 1]
    dataset_shape = shape(real_r2)

    call h5%write('datasets', 'real_r2-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r2-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 2, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_sp_rank_2


  !> Test read and write for real(sp) rank 3 arrays.
  subroutine test_xhdf5_write_and_read_real_sp_rank_3(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(sp) :: real_r3(8, 3, 5)
    real(sp), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
    integer :: offset(3), dataset_shape(3), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)

    do l = 1, size(real_r3, 3)
      do k = 1, size(real_r3, 2)
        do j = 1, size(real_r3, 1)

          real_r3(j, k, l) =  real((j - k), sp) / l

        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over third rank
    datachunk = real_r3(:, :, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = [1, 1, first]
    dataset_shape = shape(real_r3)

    call h5%write('datasets', 'real_r3-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r3-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 3, dim 3 distributed: Read array differs from written array.')


    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)

    datachunk = real_r3(:, first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = [1, first, 1]
    dataset_shape = shape(real_r3)

    call h5%write('datasets', 'real_r3-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r3-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 3, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)        
            
    datachunk = real_r3(first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = [first, 1, 1]
    dataset_shape = shape(real_r3)

    call h5%write('datasets', 'real_r3-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r3-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 3, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_sp_rank_3


  !> Test read and write for real(sp) rank 4 arrays.
  subroutine test_xhdf5_write_and_read_real_sp_rank_4(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l, m

    real(sp) :: real_r4(8, 3, 5, 4)
    real(sp), allocatable :: datachunk(:, :, :, :), datachunk_read(:, :, :, :)
    integer :: offset(4), dataset_shape(4), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)

    do m = 1, size(real_r4, 4)
      do l = 1, size(real_r4, 3)
        do k = 1, size(real_r4, 2)
          do j = 1, size(real_r4, 1)

            real_r4(j, k, l, m) = real((j - k) * (l - m), sp)
          end do
        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over fourth rank
    datachunk = real_r4(:, :, :, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._sp
    offset = [1, 1, 1, first]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-4', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-4', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 4, dim 4 distributed: Read array differs from written array.')


    ! Distribute writing over third rank
    chunk_size = nofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)

    datachunk = real_r4(:, :, first : last, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, sp)
    offset = [1, 1, first, 1]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 4, dim 3 distributed: Read array differs from written array.')

    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)

    datachunk = real_r4(:, first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, sp)
    offset = [1, first, 1, 1]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 4, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)

    datachunk = real_r4(first : last, :, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, sp)
    offset = [first, 1, 1, 1]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(sp) rank 4, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')  
  end subroutine test_xhdf5_write_and_read_real_sp_rank_4


  ! REAL(DP)


  !> Test read and write for real(dp) scalars.
  subroutine test_xhdf5_write_and_read_real_dp_rank_0(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(dp) :: real_write, real_read

    real_write = 0.213_dp

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')

    call h5%write('datasets', 'real_r1', real_write)
    call h5%read('datasets', 'real_r1', real_read)
    call test_report%assert(all_close(real_write, real_read), &
            'real(dp) rank 0: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_dp_rank_0


  !> Test read and write for real(dp) rank 1 arrays.
  subroutine test_xhdf5_write_and_read_real_dp_rank_1(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(dp) :: real_r1(8)
    real(dp), allocatable :: datachunk(:), datachunk_read(:)
    integer :: offset(1), dataset_shape(1), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)

    do l = 1, size(real_r1)

      real_r1(l) =  real(l, dp)

    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    datachunk = real_r1(first : last)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = first
    dataset_shape = shape(real_r1)

    call h5%write('datasets', 'real_r1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 1: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_dp_rank_1


  !> Test read and write for real(dp) rank 2 arrays.
  subroutine test_xhdf5_write_and_read_real_dp_rank_2(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(dp) :: real_r2(8, 3)
    real(dp), allocatable :: datachunk(:, :), datachunk_read(:, :)
    integer :: offset(2), dataset_shape(2), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)

    do l = 1, size(real_r2, 2)
      do k = 1, size(real_r2, 1)

        real_r2(k, l) =  real((l - k), dp)

      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over second rank
    datachunk = real_r2(:, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [1, first]
    dataset_shape = shape(real_r2)

    call h5%write('datasets', 'real_r2-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r2-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 2, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)

    datachunk = real_r2(first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [first, 1]
    dataset_shape = shape(real_r2)

    call h5%write('datasets', 'real_r2-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r2-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 2, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_dp_rank_2


  !> Test read and write for real(dp) rank 3 arrays.
  subroutine test_xhdf5_write_and_read_real_dp_rank_3(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    real(dp) :: real_r3(8, 3, 5)
    real(dp), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
    integer :: offset(3), dataset_shape(3), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)

    do l = 1, size(real_r3, 3)
      do k = 1, size(real_r3, 2)
        do j = 1, size(real_r3, 1)

          real_r3(j, k, l) =  real((j - k), dp) / l

        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over third rank
    datachunk = real_r3(:, :, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [1, 1, first]
    dataset_shape = shape(real_r3)

    call h5%write('datasets', 'real_r3-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r3-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 3, dim 3 distributed: Read array differs from written array.')


    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)

    datachunk = real_r3(:, first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [1, first, 1]
    dataset_shape = shape(real_r3)

    call h5%write('datasets', 'real_r3-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r3-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 3, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)        
            
    datachunk = real_r3(first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [first, 1, 1]
    dataset_shape = shape(real_r3)

    call h5%write('datasets', 'real_r3-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r3-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 3, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_real_dp_rank_3


  !> Test read and write for real(dp) rank 4 arrays.
  subroutine test_xhdf5_write_and_read_real_dp_rank_4(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l, m

    real(dp) :: real_r4(8, 3, 5, 4)
    real(dp), allocatable :: datachunk(:, :, :, :), datachunk_read(:, :, :, :)
    integer :: offset(4), dataset_shape(4), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)

    do m = 1, size(real_r4, 4)
      do l = 1, size(real_r4, 3)
        do k = 1, size(real_r4, 2)
          do j = 1, size(real_r4, 1)

            real_r4(j, k, l, m) = real((j - k) * (l - m), dp)
          end do
        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over fourth rank
    datachunk = real_r4(:, :, :, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [1, 1, 1, first]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-4', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-4', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 4, dim 4 distributed: Read array differs from written array.')


    ! Distribute writing over third rank
    chunk_size = nofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)

    datachunk = real_r4(:, :, first : last, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, 1, first, 1]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 4, dim 3 distributed: Read array differs from written array.')

    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)

    datachunk = real_r4(:, first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, first, 1, 1]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 4, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)

    datachunk = real_r4(first : last, :, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [first, 1, 1, 1]
    dataset_shape = shape(real_r4)

    call h5%write('datasets', 'real_r4-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'real_r4-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 4, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')  
  end subroutine test_xhdf5_write_and_read_real_dp_rank_4


  ! COMPLEX(DP)


  !> Test read and write for complex(dp) rank 1 arrays.
  subroutine test_xhdf5_write_and_read_cmplx_dp_rank_1(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    complex(dp) :: cmplx_r1(8)
    complex(dp), allocatable :: datachunk(:), datachunk_read(:)
    integer :: offset(1), dataset_shape(1), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(cmplx_r1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r1), mpiglobal%procs)

    do l = 1, size(cmplx_r1)

      cmplx_r1(l) =  cmplx(l, -l, dp)

    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')

    datachunk = cmplx_r1(first : last)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = first
    dataset_shape = shape(cmplx_r1)

    call h5%write('datasets', 'cmplx_r1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex rank 1: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_cmplx_dp_rank_1


  !> Test read and write for complex(dp) rank 2 arrays.
  subroutine test_xhdf5_write_and_read_cmplx_dp_rank_2(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    complex(dp) :: cmplx_r2(8, 3)
    complex(dp), allocatable :: datachunk(:, :), datachunk_read(:, :)
    integer :: offset(2), dataset_shape(2), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(cmplx_r2, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r2, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r2, 2), mpiglobal%procs)

    do l = 1, size(cmplx_r2, 2)
      do k = 1, size(cmplx_r2, 1)

        cmplx_r2(k, l) =  cmplx((l + k), (l - k), dp)

      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over second rank
    datachunk = cmplx_r2(:, first : last)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [1, first]
    dataset_shape = shape(cmplx_r2)

    call h5%write('datasets', 'cmplx_r2-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r2-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 2, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(cmplx_r2, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r2, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r2, 1), mpiglobal%procs)

    datachunk = cmplx_r2(first : last, :)
    datachunk_read = datachunk
    datachunk_read = 0._dp
    offset = [first, 1]
    dataset_shape = shape(cmplx_r2)

    call h5%write('datasets', 'cmplx_r2-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r2-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'real(dp) rank 2, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine test_xhdf5_write_and_read_cmplx_dp_rank_2


  !> Test read and write for complex(dp) rank 3 arrays.
  subroutine test_xhdf5_write_and_read_cmplx_dp_rank_3(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l

    complex(dp) :: cmplx_r3(8, 3, 5)
    complex(dp), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
    integer :: offset(3), dataset_shape(3), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(cmplx_r3, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r3, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r3, 3), mpiglobal%procs)

    do l = 1, size(cmplx_r3, 3)
      do k = 1, size(cmplx_r3, 2)
        do j = 1, size(cmplx_r3, 1)

          cmplx_r3(j, k, l) = cmplx((j - k) * l, (k - j) * l, dp)

        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over third rank
    datachunk = cmplx_r3(:, :, first : last)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, 1, first]
    dataset_shape = shape(cmplx_r3)

    call h5%write('datasets', 'cmplx_r3-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r3-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 3, dim 3 distributed: Read array differs from written array.')


    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(cmplx_r3, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r3, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r3, 2), mpiglobal%procs)

    datachunk = cmplx_r3(:, first : last, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, first, 1]
    dataset_shape = shape(cmplx_r3)

    call h5%write('datasets', 'cmplx_r3-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r3-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 3, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(cmplx_r3, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r3, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r3, 1), mpiglobal%procs)        
            
    datachunk = cmplx_r3(first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [first, 1, 1]
    dataset_shape = shape(cmplx_r3)

    call h5%write('datasets', 'cmplx_r3-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r3-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 3, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
  end subroutine

  
  !> Test read and write for complex(dp) rank 4 arrays.
  subroutine test_xhdf5_write_and_read_cmplx_dp_rank_4(test_report, mpiglobal)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> MPI communicator
    type(mpiinfo), intent(in) :: mpiglobal

    type(xhdf5_type) :: h5
    integer :: ierr, j, k ,l, m

    complex(dp) :: cmplx_r4(8, 3, 5, 4)
    complex(dp), allocatable :: datachunk(:, :, :, :), datachunk_read(:, :, :, :)
    integer :: offset(4), dataset_shape(4), chunk_size, first, last

    chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 4), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r4, 4), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r4, 4), mpiglobal%procs)

    do m = 1, size(cmplx_r4, 4)
      do l = 1, size(cmplx_r4, 3)
        do k = 1, size(cmplx_r4, 2)
          do j = 1, size(cmplx_r4, 1)

            cmplx_r4(j, k, l, m) = cmplx((j - k) * (l - m), (k - j) * (m - l), dp)

          end do
        end do
      end do 
    end do 

    call h5%initialize('test_file.h5', mpiglobal)
    call h5%initialize_group('./', 'datasets')


    ! Distribute writing over fourth rank
    datachunk = cmplx_r4(:, :, :, first : last)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, 1, 1, first]
    dataset_shape = shape(cmplx_r4)

    call h5%write('datasets', 'cmplx_r4-4', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r4-4', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 4, dim 4 distributed: Read array differs from written array.')


    ! Distribute writing over third rank
    chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 3), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r4, 3), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r4, 3), mpiglobal%procs)

    datachunk = cmplx_r4(:, :, first : last, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, 1, first, 1]
    dataset_shape = shape(cmplx_r4)

    call h5%write('datasets', 'cmplx_r4-3', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r4-3', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 4, dim 3 distributed: Read array differs from written array.')


    ! Distribute writing over second rank
    chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 2), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r4, 2), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r4, 2), mpiglobal%procs)

    datachunk = cmplx_r4(:, first : last, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [1, first, 1, 1]
    dataset_shape = shape(cmplx_r4)

    call h5%write('datasets', 'cmplx_r4-2', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r4-2', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 4, dim 2 distributed: Read array differs from written array.')


    ! Distribute writing over first rank
    chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 1), mpiglobal%procs)
    first = firstofset(mpiglobal%rank, size(cmplx_r4, 1), mpiglobal%procs)
    last = lastofset(mpiglobal%rank, size(cmplx_r4, 1), mpiglobal%procs)

    datachunk = cmplx_r4(first : last, :, :, :)
    datachunk_read = datachunk
    datachunk_read = cmplx(0.0, 0.0, dp)
    offset = [first, 1, 1, 1]
    dataset_shape = shape(cmplx_r4)

    call h5%write('datasets', 'cmplx_r4-1', datachunk, offset, dataset_shape)
    call h5%read('datasets', 'cmplx_r4-1', datachunk_read, offset)
    call test_report%assert(all_close(datachunk, datachunk_read), &
            'complex(dp) rank 4, dim 1 distributed: Read array differs from written array.')

    call h5%finalize()
    if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')  
  end subroutine test_xhdf5_write_and_read_cmplx_dp_rank_4
  
end module