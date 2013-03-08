subroutine cpFileToNodes(filnam)
    Use modmain
    Use modmpi
    Use m_getunit
    implicit none
    ! interface
    Character (*), Intent (In) :: filnam
    integer::recl
    !local variables
#ifdef MPI
    character,allocatable :: buffer(:)
    integer::nstsv_,ik,un,finrank,finprocs, filesize,chunk,chunksize
    real(8)::vkl_(3)
    Call getunit (un)
    if(firstinnode)then
    Call mpi_comm_rank (firstinnode_comm, finrank, ierr)
    call mpi_comm_size(firstinnode_comm, finprocs, ierr)
     if(  rank.eq.0.and.finrank.ne.0 ) then
     write(*,*) "Sorry wrong assumption about sub comm numbering"
     write(*,*) "world rank:",rank,"firstinnode",finrank
     endif

    chunksize=1000
    allocate(buffer( chunksize))

    inquire (file=trim(filnam), size=filesize)
    write(*,*)"bevor bcast"
    call  MPI_BCAST(filesize, 1, MPI_INTEGER, 0, firstinnode_comm, ierr)
    write(*,*)"after.."
    if(rank.eq. 0)then
        Open (Unit=un, File=trim(filnam), Form='unformatted', &
        Action='read', Access='direct', Recl=chunksize/4)
    else
        Open (Unit=un, File=trim(filnam), Form='unformatted', &
        Action='write', Access='direct', Recl=chunksize/4)
    endif

    do ik=1,filesize/chunksize +1
        if(ik.le.filesize/chunksize) then
            chunk=chunksize
        else
            chunk=mod(filesize,chunksize)
        endif
        if(chunk.ne.0)then
            if (rank .eq. 0) then
                read (un,Rec=ik) buffer(1:chunk)
            endif
            call  MPI_BCAST(buffer, chunk, MPI_Character, 0, firstinnode_comm, ierr)
            if (rank .ne. 0) then
                Write (un,Rec=ik) buffer(1:chunk)
            endif
        endif
    end do
    Close (un)
    endif
#endif
end subroutine
