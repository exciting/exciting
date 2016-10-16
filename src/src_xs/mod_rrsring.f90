module mod_rrsring
  use modmpi
  use modscl

  implicit none

  private 

  integer(4) :: mytarget, mysource
  integer(4) :: i, id, ishift, ndj, nall, rest
  integer(4), allocatable :: myreadlist(:), myreceivelist(:,:), mysendlist(:,:)

  public :: nall, rest
  public :: ndj
  public :: mytarget, mysource
  public :: myreceivelist, mysendlist
  public :: setup_rrslists

  contains 

    subroutine setup_rrslists(nelements)
      integer(4), intent(in) :: nelements

      ! Receive cycles
      !   The nelemtes registers of some file are read by npoc processes
      !   in N or N+1 read cycles. Where nelements/nproc = N rest R, so
      !   that in N cycles every process reads a register, while in
      !   the N+1'th and last cycle only the first R processes read 
      !   from file.
      !   After the read is done, the data gets shifted nporc-1 times
      !   from neighbour to next neighbour in one direction (along 
      !   a process ring, in mathematically negative direction),
      !   since each process needs data form every 
      !   file register to build the local matrix.
      ndj  = ceiling(dble(nelements)/dble(nproc)) 
      nall = nelements/nproc
      rest = mod(nelements, nproc)

      ! Determine from which rank to receive data (always the process previous in ring)
      mysource = mypcol1d_r - 1 
      if(mysource < 0) mysource = nproc -1 

      ! Determine to which rank to send data (always next process in ring)
      mytarget = mypcol1d_r + 1
      if(mytarget > nproc-1) mytarget = 0 

      ! The registers to read by each process are distributed
      ! around the process ring one by one.
      if(allocated(myreadlist)) deallocate(myreadlist)
      allocate(myreadlist(ndj))
      do id = 1, ndj
        myreadlist(id) = mypcol1d_r + (id-1)*nproc + 1
        ! An entry of 0 means no read from file
        ! in cycle N+1.
        if(myreadlist(id) > nelements) then
          myreadlist(id) = 0
        end if
      end do

      !!*************************************************!!
      !! Make lists of which elements to receive         !!
      !! each cycle and which to send.                   !!
      !!*************************************************!!
      if(allocated(myreceivelist)) deallocate(myreceivelist)
      allocate(myreceivelist(ndj, nproc))
      if(allocated(mysendlist)) deallocate(mysendlist)
      allocate(mysendlist(ndj, nproc))
      ! If there are more elements then processes:
      ! All processes read new data and send it around
      do id = 1, nall
        myreceivelist(id, 1) = 0
        mysendlist(id, 1) = myreadlist(id)
        do ishift = 2, nproc
          myreceivelist(id, ishift) = mysendlist(id,ishift-1) - 1
          if(mod(myreceivelist(id, ishift), nproc) == 0) then
            myreceivelist(id, ishift) = myreceivelist(id, ishift) + nproc
          end if
          mysendlist(id, ishift) = myreceivelist(id, ishift)
        end do
        if(nproc /= 1) then 
          mysendlist(id, nproc) = -1
        end if
      end do
      ! If there are less elements then processes or the processes
      ! do not divide the elements without rest.
      ! Adjust for potential last (or only) partial data injection.
      if(rest > 0) then
        myreceivelist(ndj,:) = -1
        do i = 1, rest
          myreceivelist(ndj, i) = nall*nproc + (rest-i+1)
        end do
        myreceivelist(ndj, :) = cshift(myreceivelist(ndj,:), rest-mypcol1d_r-1)
        mysendlist(ndj,:) = myreceivelist(ndj, :)
        if(myreadlist(ndj) > 0) then
          myreceivelist(ndj, 1) = 0
        end if
        if(mypcol1d_r+1 < rest) then
          mysendlist(ndj, nproc) = -1
        end if
        if(mypcol1d_r == nproc-1) then 
          mysendlist(ndj, nproc) = -1
        end if
      end if

      deallocate(myreadlist)

    end subroutine setup_rrslists

end module mod_rrsring
