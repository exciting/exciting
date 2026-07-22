program test_mpi4_in_place
    use mpi_f08
    implicit none
    integer(kind = MPI_COUNT_KIND), allocatable :: recvcounts(:)
    integer(kind = MPI_ADDRESS_KIND), allocatable :: displs(:)
    integer, allocatable :: recvbuf(:)
    integer :: ierr, rank, comm_size, i

    call MPI_Init( ierr )
    call MPI_Comm_size( MPI_COMM_WORLD, comm_size, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )

    allocate( recvcounts(comm_size) )
    allocate( displs(comm_size) )
    allocate( recvbuf(comm_size) )

    ! ==============================================================================
    ! BUG CHECK: MPI_IN_PLACE Pointer Translation in MPI-4.0 Large Count Wrappers
    ! ==============================================================================
    ! Context:
    ! The Fortran wrapper may fail to correctly translate the Fortran MPI_IN_PLACE 
    ! sentinel to the underlying C MPI_IN_PLACE sentinel when delegating the call. 
    ! Consequently, the C layer treats MPI_IN_PLACE as a standard, valid memory 
    ! address and attempts to read from or write to it.
    !
    ! How this test works:
    ! We populate small arrays and call MPI_Allgatherv using MPI_IN_PLACE with the
    ! MPI_DATATYPE_NULL as data type needed to be sent. If the wrapper is broken, this call crashes.
    ! ==============================================================================
    do i = 1, comm_size
        recvcounts(i) = 1_MPI_COUNT_KIND
        displs(i) = int( i - 1, kind = MPI_ADDRESS_KIND )
    end do
    recvbuf = rank
    call MPI_Allgatherv( MPI_IN_PLACE, 0_MPI_COUNT_KIND, MPI_DATATYPE_NULL, &
      recvbuf, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr )
    call MPI_Finalize( ierr )
end program test_mpi4_in_place
