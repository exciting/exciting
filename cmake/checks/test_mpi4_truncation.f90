program test_mpi4_truncation
    use mpi_f08
    implicit none
    integer(kind = MPI_COUNT_KIND), allocatable :: recvcounts(:)
    integer(kind = MPI_ADDRESS_KIND), allocatable :: displs(:)
    type(MPI_Datatype) :: zero_type
    integer :: ierr, comm_size, i
    
    ! Distinct variables to prevent aliasing aborts
    integer, target :: send_dummy
    integer, target :: recv_dummy

    call MPI_Init( ierr )
    call MPI_Comm_size( MPI_COMM_WORLD, comm_size, ierr )
    allocate( recvcounts(comm_size) )
    allocate( displs(comm_size) )

    ! BUG CHECK: Sendcount Truncation (64-bit scalar)
    call MPI_Type_contiguous( 0, MPI_INTEGER, zero_type, ierr )
    call MPI_Type_commit( zero_type, ierr )
    do i = 1, comm_size
        recvcounts(i) = 6861637225_MPI_COUNT_KIND
        displs(i) = 0_MPI_ADDRESS_KIND
    end do
    call MPI_Allgatherv( send_dummy, 6861637225_MPI_COUNT_KIND, zero_type, &
      recv_dummy, recvcounts, displs, zero_type, MPI_COMM_WORLD, ierr ) 
    call MPI_Type_free( zero_type, ierr )
    call MPI_Finalize( ierr )
end program test_mpi4_truncation