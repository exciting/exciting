! Copyright (C) 2006-2008 C. Ambrosch-Draxl. C. Meisenbichler S. Sagmeister
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
! !MODULE:  modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   In case of compiled without MPI support it defines the
!   mpi specific variables such that the code behaves exactly as
!   the unmodified scalar version
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (S. Sagmeister)
!   Added allgatherv interface, August 2010 (S. Sagmeister)
!   Added subroutines/functions to documentation scheme, 2016 (Aurich)
!
module modmpi
#ifdef MPI
  use mpi
#endif

  implicit none
 
  integer(4) :: rank
  integer(4) :: procs
  integer(4) :: ierr, procnamelen, nnodes
  integer(4) :: world_group, firstinnode_group, firstinnode_comm
  integer(4), allocatable :: nodechefs(:)
  logical :: splittfile, firstinnode
  
  integer(4), parameter :: tag = 25 

  private :: tag

  contains

    !BOP
    ! !ROUTINE:
    ! !INTERFACE:
    subroutine initmpi
    ! !DESCRIPTION:
    !   Initializes MPI and sets procs and rank numbers.
    !   Sets splitfile and initializes the first in node
    !   list. Or if -DMPI is not used sets procs=1 and rank=1
    !   and splitfile=.false.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme (Aurich)
    !EOP
    !BOC
#ifdef MPI
      ! Initialize MPI and get number of processes 
      ! and current rank
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, procs, ierr)
      call mpi_comm_rank(mpi_comm_world, rank, ierr)

      nnodes = 1
      splittfile = .true.
      firstinnode = .true.

      allocate(nodechefs(procs))
      call get_isfirstinnode(200*procs)
#endif
#ifndef MPI
      procs = 1
      rank = 0
      splittfile = .false.
#endif
    end subroutine initmpi
    !EOC

    !BOP
    ! !ROUTINE: get_isfirstinnode
    ! !INTERFACE:
    subroutine get_isfirstinnode(strsize)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: strsize ! Size of neighor string
    !
    ! !DESCRIPTION:
    ! The routine uses the processor name returned by
    ! {\tt mpi\_get\_processor\_name} to determine the first rank
    ! on each node (nodechef). For those ranks a mpi group
    ! and communicator is created. Also sets the firstinnode
    ! flag in the modmpi module.
    !
    ! !REVISION HISTORY:
    ! Added to documentation scheme. (Aurich)
    !EOP
    !BOC

      integer(4), intent(in)::strsize
#ifdef MPI
      integer(4) :: request, recvstatus(mpi_status_size), i
      character(len = strsize) :: neighbors, neighborssend
      character(200) :: procname
      logical :: lbuffer

      procname = ''
      neighbors = ''
      neighborssend = ''

      call mpi_get_processor_name(procname, procnamelen, ierr)

      ! The following code is illustrated by the example 
      ! with 4 threads and 2 physical processors A and B.
      !   rank 0 has an empty neighbors string (ns) "" and a procname (pn) of "A", 
      !   it sends ",A" to rank 1
      !   rank 1 now has ns=",A" and pn="A", it sends ",A,A" to rank 2
      !   rank 2 now has ns=",A,A" and pn="B", it sends ",A,A,B" to rank 3
      !   rank 3 now has ns=",A,A,B" and pn="B", it has no one to send to.
      !   Each rank checks if pn is in ns, if not it is the nodechef.
      if(rank .gt. 0) then
        call mpi_recv(neighbors, strsize, mpi_character,&
          & rank-1, tag, mpi_comm_world, recvstatus, ierr)
      end if
      if(rank .lt. procs-1) then
        write(neighborssend,*) adjustl(trim(neighbors))//","//adjustl(trim(procname))
        call mpi_send(neighborssend, strsize, mpi_character,&
          & rank+1, tag, mpi_comm_world, ierr)
      endif
      if(index(neighbors, adjustl(trim(procname)) ).gt.0) then
        firstinnode = .false.
      else
        firstinnode = .true.
      endif

      ! The following collects all the rank numbes that are nodechefs 
      ! in the first 1:nnodes elements of nodechefs and chreates
      ! a mpi group and communicator for internode communication.
      nnodes = 0
      do i = 0, procs-1
        lbuffer = firstinnode
        call mpi_bcast(lbuffer, 1, mpi_logical, i, mpi_comm_world, ierr)
        if(lbuffer)then
           nnodes = nnodes+1
           nodechefs(nnodes)= i
        endif
      end do
      call mpi_comm_group(mpi_comm_world, world_group, ierr)
      call mpi_group_incl(world_group, nnodes, nodechefs, firstinnode_group, ierr)
      call mpi_comm_create(mpi_comm_world, firstinnode_group, firstinnode_comm, ierr)
#endif
    end subroutine
    !EOC

    !BOP
    ! !ROUTINE: finitmpi
    ! !INTERFACE: 
    subroutine finitmpi
    ! !DESCRIPTION:
    !   If -DMPI calls {\tt mpi\_finalize}, otherwise does nothing.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
#ifdef MPI
      call mpi_finalize(ierr)
#endif
    end subroutine finitmpi
    !EOC


!-------------generalized partition!
    !BOP
    ! !ROUTINE: nofset
    ! !INTERFACE:
    function nofset(process, set)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: process ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! OUT:
    ! integer(4) :: nofset  ! Number of elements for that rank
    !
    ! !DESCRIPTION:
    !   Calculates number of items for proc.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      integer(4) :: nofset
      integer(4), intent(in) :: process, set

      if(process .lt. 0) then
        nofset = 0
      else
        nofset = set / procs
        if((mod(set, procs) .gt. process)) nofset = nofset + 1
      end if
    end function nofset
    !EOC

    !BOP
    ! !ROUTINE: firstofset
    ! !INTERFACE:
    function firstofset(process, set)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: process ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! OUT:
    ! integer(4) :: firstofset ! Index of the total set for the first index 
    !                          ! of the current subset
    !
    ! !DESCRIPTION:
    !   Computes the total element index for first element of the subset of
    !   the current rank.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      integer(4) :: firstofset
      integer(4), intent(in) :: process, set
      integer(4) :: i

      firstofset = 1
      do i = 0, min( process-1, set-1)
         firstofset = firstofset + nofset(i, set)
      end do
      if(set .le. process) firstofset = set
    end function firstofset
    !EOC

    !BOP
    ! !ROUTINE: lastofset
    ! !INTERFACE:
    function lastofset(process, set)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: process ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! OUT:
    ! integer(4) :: lastofset  ! Index of the total set for the first index 
    !                          ! of the current subset
    !
    ! !DESCRIPTION:
    !   Computes the total element index for last element of the subset of
    !   the current rank.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      integer(4) :: lastofset
      integer(4), intent(in) :: process, set
      integer(4)::i
      lastofset = 0
      do i = 0, process
         lastofset = lastofset + nofset(i, set)
      end do
    end function lastofset
    !EOC

    !BOP
    ! !ROUTINE: procofindex
    ! !INTERFACE:
    function procofindex(k, set)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: k    ! MPI rank
    ! integer(4) :: set  ! Total number of elements to distribute
    ! OUT:
    ! integer(4) :: procofindex  ! Rank that holds the subset where the 
    !                            ! global index k is included.
    !
    ! !DESCRIPTION:
    !   Computed the rank than holds the subset of elements in which 
    !   the element k of the total set is included.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      integer(4) :: procofindex
      integer(4), intent(in) :: k, set
      integer(4) :: iproc
      procofindex = 0
      do iproc = 0, procs - 1
         if(k .gt. lastofset(iproc, set)) procofindex = procofindex + 1
      end do
    end function procofindex
    !EOC

    !BOP
    ! !ROUTINE: lastproc
    ! !INTERFACE:
    function lastproc(row, set)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: row    
    ! integer(4) :: set    ! Total number of elements to distribute
    ! OUT:
    ! integer(4) :: lastproc 
    !                        
    !
    ! !DESCRIPTION:
    !   (??)
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4) :: lastproc
      integer(4), intent(in) :: row, set
      if(row .ne. nofset(0, set)) then
         lastproc = procs
      else
         lastproc = modulo(set, procs)
         if(lastproc .eq. 0) lastproc = procs
      end if
      lastproc = lastproc - 1
    end function lastproc
    !EOC

!------------------interface to mpi_barrier for xs-part
    !BOP
    ! !ROUTINE: barrier
    ! !INTERFACE:
    subroutine barrier
    ! !DESCRIPTION:
    !   If -DMPI calls {\tt mpi\_barrier}, else nothing.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      implicit none
      ! do nothing if only one process
#ifndef MPI
      if(procs .eq. 1) return
#endif
      ! call the mpi barrier
#ifdef MPI
      call mpi_barrier(mpi_comm_world, ierr)
#endif
    end subroutine barrier
    !EOC

    !BOP
    ! !ROUTINE: endloopbarrier
    ! !INTERFACE:
    subroutine endloopbarrier(set, mult)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: set
    ! integer(4) :: mult
    !
    ! !DESCRIPTION:
    !   (??)
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: set, mult
      integer(4) :: i
      do i = 1, (nofset(0, set)-nofset(rank, set)) * mult
         call barrier
      end do
    end subroutine endloopbarrier
    !EOC

!------------------wrappers for mpi communication

    !BOP
    ! !ROUTINE: endloopbarrier
    ! !INTERFACE:
    subroutine mpi_allgatherv_ifc(set, rlen, ibuf, rbuf, rlpbuf, zbuf)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: set    ! Number of elements in distributed set 
    ! integer(4) :: rlen   ! Number of data elements in such an element
    ! Inout:
    ! integer(4), optional :: ibuf(*) ! Buffers to send/recive
    ! real(8), optional :: rbuf(*)    !
    ! real(4), optional :: rlbuf(*)   !
    ! complex(8), optional :: zbuf(*) !
    !
    ! !DESCRIPTION:
    !   Wrapper routine for {\tt MPI\_ALLGATHER} for different 
    !   data types which is adapted for the k-point set 
    !   distribution scheme.
    !   Note: Needs -DMPI1 (??)
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: set, rlen
      integer(4), intent(inout), optional :: ibuf(*)
      real(8), intent(inout), optional :: rbuf(*)
      real(4), intent(inout), optional :: rlpbuf(*)
      complex(8), intent(inout), optional :: zbuf(*)
#ifndef MPI1
#define BUFFER mpi_in_place
#endif
#ifdef MPI1
      complex(8), allocatable :: bufz(:)
      real(8), allocatable :: bufr(:)
      real(4), allocatable :: bufrlp(:)
      integer(4), allocatable :: bufi(:)
#endif
      integer(4), allocatable :: buf_n(:), buf_dspls(:)
      integer(4) :: j
      logical :: ti, tr, trlp, tz
      
      ti = present(ibuf)
      tr = present(rbuf)
      trlp = present(rlpbuf)
      tz = present(zbuf)

      if(count((/ti, tr, trlp, tz/)).ne.1) then
        write(*,*)
        write(*,'("Error (mpi_allgatherv_ifc): Exactly one array must be defined.")')
        write(*,*)
        stop
      end if

#ifdef MPI
      allocate(buf_n(procs), buf_dspls(procs))

      ! number of elements in send buffer (flattened array)
      buf_n =(/(rlen*nofset(j, set), j = 0, procs-1)/)

      ! Complex doubles

      ! displacements within receive buffer(flattened array)
      buf_dspls =(/(rlen*(firstofset(j, set)-1), j = 0, procs-1)/)

      ! use recieve buffer as sendbuffer by specifying mpi_in_place

      ! Integers
      if(ti) then
#ifdef MPI1
#define BUFFER bufi
        allocate(bufi(buf_n(rank+1)))
        bufi(:)= ibuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
#endif
        call mpi_allgatherv(BUFFER, &
          buf_n(rank+1), &
          mpi_integer, &
          ibuf, &
          buf_n, &
          buf_dspls, &
          mpi_integer, &
          mpi_comm_world, &
          ierr)
#ifdef MPI1
        deallocate(bufi)
#undef BUFFER
#endif
      end if

      ! Doubles
      if(tr) then
#ifdef MPI1
#define BUFFER bufr
        allocate(bufr(buf_n(rank+1)))
        bufr(:)= rbuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
#endif
        call mpi_allgatherv(BUFFER, &
          buf_n(rank+1), &
          mpi_double_precision, &
          rbuf, &
          buf_n, &
          buf_dspls, &
          mpi_double_precision, &
          mpi_comm_world, &
          ierr)
#ifdef MPI1
        deallocate(bufr)
#undef BUFFER
#endif
      end if

      ! Floats
      if(trlp) then
#ifdef MPI1
#define BUFFER bufrlp
        allocate(bufrlp(buf_n(rank+1)))
        bufrlp(:)= rlpbuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
#endif
        call mpi_allgatherv(BUFFER, &
          buf_n(rank+1), &
          mpi_real4, &
          rlpbuf, &
          buf_n, &
          buf_dspls, &
          mpi_real4, &
          mpi_comm_world, &
          ierr)
#ifdef MPI1
        deallocate(bufrlp)
#undef BUFFER
#endif
      end if

      ! Complex doubles
      if(tz) then
#ifdef MPI1
#define BUFFER bufz
        allocate(bufz(buf_n(rank+1)))
        bufz(:)= zbuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
#endif
        call mpi_allgatherv(BUFFER, &
          buf_n(rank+1), &
          mpi_double_complex, &
          zbuf, &
          buf_n, &
          buf_dspls, &
          mpi_double_complex, &
          mpi_comm_world, &
          ierr)
#ifdef MPI1
        deallocate(bufz)
#undef BUFFER
#endif
      end if
      deallocate(buf_n, buf_dspls)
#endif
    end subroutine
    !EOC

! Service functions to partition k points
! still kept around but should be replaced by generig functions
! firstofset nofset lastofset ....
    function nofk(process)
      use modmain, only: nkpt
      integer(4) :: nofk
      integer(4), intent(in) :: process
      nofk = nkpt / procs
      if((mod(nkpt, procs) .gt. process)) nofk = nofk + 1
    end function nofk

    function firstk(process)
       integer(4) :: firstk
       integer(4), intent(in) :: process
       integer(4)::i
       firstk = 1
       do i = 0, process - 1
          firstk = firstk + nofk(i)
       end do
    end function firstk

    function lastk(process)
       integer(4) :: lastk, i
       integer(4), intent(in) :: process
       lastk = 0
       do i = 0, process
          lastk = lastk + nofk(i)
       end do
    end function lastk

    function procofk(k)
       integer(4) :: procofk
       integer(4), intent(in) :: k
       integer(4) :: iproc
       procofk = 0
       do iproc = 0, procs - 1
          if(k .gt. lastk(iproc)) procofk = procofk + 1
       end do
    end function procofk


end module modmpi
