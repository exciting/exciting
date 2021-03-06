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
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (S. Sagmeister)
!   Added allgatherv interface, August 2010 (S. Sagmeister)
!
!
!

Module modmpi

#ifdef MPI
use mpi
#endif
implicit none
 
      Integer :: rank
      Integer :: procs
      Integer :: ierr,procnamelen,nnodes
      integer:: world_group,firstinnode_group,firstinnode_comm
      integer, allocatable:: nodechefs(:)
      Logical :: splittfile,firstinnode

!
!!$  character(256)::filename
Contains
      Subroutine initMPI
#ifdef MPI
    !        mpi init

         Call mpi_init (ierr)
         Call mpi_comm_size (mpi_comm_world, procs, ierr)
         Call mpi_comm_rank (mpi_comm_world, rank, ierr)
         nnodes=1
         splittfile = .True.
         firstinnode=.True.
         allocate(nodechefs(procs))
         call get_isfirstinnode(200*procs)
#endif
#ifndef MPI
         procs = 1
         rank = 0
         splittfile = .False.
#endif
      End Subroutine initMPI
!
subroutine get_isfirstinnode(strsize)
         integer,intent(in)::strsize
#ifdef MPI

	     integer::tag,request,recvstatus (MPI_STATUS_SIZE),i
	     character(len=strsize)::neighbors,neighborssend
	     character(200)::procname
	     logical::lbuffer
         procname=""
         neighbors=""
         neighborssend=""
         call MPI_Get_processor_name(procname, procnamelen,ierr)

         tag=25

         if (rank.gt.0)  call MPI_recv(neighbors,strsize,MPI_CHARACTER,rank-1,tag,&
         mpi_comm_world,recvstatus,ierr)
         if (rank.lt.procs-1)then
            write(neighborssend,*) adjustl(trim(neighbors))//","//adjustl(trim(procname))
      		call MPI_send(neighborssend,strsize,MPI_CHARACTER,rank+1,tag,&
         	&mpi_comm_world,ierr)
		 endif
         if(index(neighbors,adjustl(trim(procname)) ).gt.0) then
         firstinnode=.false.
         else
         firstinnode=.true.
         endif
         nnodes=0
         do i=0,procs-1
         lbuffer=firstinnode
         call mpi_bcast(lbuffer,1,mpi_logical,i,mpi_comm_world,ierr)

         if(lbuffer)then
            nnodes=nnodes+1
            nodechefs(nnodes)=i
         endif
         end do
         call MPI_Comm_group(MPI_Comm_WORLD, world_group,ierr)
         call MPI_GROUP_INCL(world_group, nnodes, nodechefs, firstinnode_group,ierr)
         call MPI_Comm_create(MPI_Comm_World, firstinnode_group, firstinnode_comm,ierr)
#endif
end subroutine
      Subroutine finitMPI
#ifdef MPI
         Call MPI_Finalize (ierr)
#endif
      End Subroutine finitMPI
!
! service functions to partition k points
! still kept around but should be replaced by generig functions
! firstofset nofset lastofset ....
      Function nofk (process)
         Use modmain, Only: nkpt
         Integer :: nofk
         Integer, Intent (In) :: process
         nofk = nkpt / procs
         If ((Mod(nkpt, procs) .Gt. process)) nofk = nofk + 1
      End Function nofk
!
      Function firstk (process)
         Integer :: firstk
         Integer, Intent (In) :: process
         integer::i
         firstk = 1
         Do i = 0, process - 1
            firstk = firstk + nofk (i)
         End Do
      End Function firstk
!
      Function lastk (process)
         Integer :: lastk,i
         Integer, Intent (In) :: process
         lastk = 0
         Do i = 0, process
            lastk = lastk + nofk (i)
         End Do
      End Function lastk
!
      Function procofk (k)
         Integer :: procofk
         Integer, Intent (In) :: k
         Integer :: iproc
         procofk = 0
         Do iproc = 0, procs - 1
            If (k .Gt. lastk(iproc)) procofk = procofk + 1
         End Do
      End Function procofk
!
!-------------generalized partition!
      Function nofset (process, set)
      !calculates number of items for proc
         Integer :: nofset
         Integer, Intent (In) :: process, set
         if (process.lt.0) then
         nofset=0
         else
         nofset = set / procs
         If ((Mod(set, procs) .Gt. process)) nofset = nofset + 1
		 endif
      End Function nofset
!
      Function firstofset (process, set)
         Integer :: firstofset
         Integer, Intent (In) :: process, set
         integer::i
         firstofset = 1
         Do i = 0, min( process-1,set-1)
            firstofset = firstofset + nofset (i, set)
         End Do
         if(set.le.process) firstofset=set
      End Function firstofset
!
      Function lastofset (process, set)
         Integer :: lastofset
         Integer, Intent (In) :: process, set
         integer::i
         lastofset = 0
         Do i = 0, process
            lastofset = lastofset + nofset (i, set)
         End Do
      End Function lastofset
!
      Function procofindex (k, set)
         Integer :: procofindex
         Integer, Intent (In) :: k, set
         Integer :: iproc
         procofindex = 0
         Do iproc = 0, procs - 1
            If (k .Gt. lastofset(iproc, set)) procofindex = procofindex &
           & + 1
         End Do
      End Function procofindex
!
      Function lastproc (row, set)
         Implicit None
         Integer :: lastproc
         Integer, Intent (In) :: row, set
         If (row .Ne. nofset(0, set)) Then
            lastproc = procs
         Else
            lastproc = modulo (set, procs)
            If (lastproc .Eq. 0) lastproc = procs
         End If
         lastproc = lastproc - 1
      End Function lastproc
!
!
!
!------------------interface to MPI_barrier for xs-part
      Subroutine barrier
         Implicit None
  ! do nothing if only one process
#ifndef MPI
         If (procs .Eq. 1) Return
#endif
  ! call the MPI barrier
#ifdef MPI
         Call MPI_barrier (mpi_comm_world, ierr)
#endif
      End Subroutine barrier


      Subroutine endloopbarrier (set, mult)
         Implicit None
         Integer, Intent (In) :: set, mult
         Integer :: i
         Do i = 1, (nofset(0, set)-nofset(rank, set)) * mult
            Call barrier
         End Do
      End Subroutine endloopbarrier
!
!------------------wrappers for MPI communication

      subroutine mpi_allgatherv_ifc(set,rlen,ibuf,rbuf,rlpbuf,zbuf)
        implicit none
        integer, intent(in) :: set,rlen
        integer, intent(inout), optional :: ibuf(*)
        real(8), intent(inout), optional :: rbuf(*)
        real(4), intent(inout), optional :: rlpbuf(*)
        complex(8), intent(inout), optional :: zbuf(*)
        ! local variables
#ifndef MPI1
#define BUFFER mpi_in_place
#endif
#ifdef MPI1
        complex(8), allocatable :: bufz(:)
        real(8), allocatable :: bufr(:)
        real(4), allocatable :: bufrlp(:)
        integer, allocatable :: bufi(:)
#endif
        integer, allocatable :: buf_n(:),buf_dspls(:)
        integer :: j
        logical :: ti,tr,trlp,tz
        ti=present(ibuf); tr=present(rbuf); trlp=present(rlpbuf); tz=present(zbuf)
        if (count((/ti,tr,trlp,tz/)).ne.1) then
          write(*,*)
          write(*,'("Error(mpi_allgatherv_ifc): exactly one array must be defined.")')
          write(*,*)
          stop
        end if
#ifdef MPI
        allocate(buf_n(procs),buf_dspls(procs))
        ! number of elements in send buffer (flattened array)
        buf_n=(/(rlen*nofset(j,set),j=0,procs-1)/)
        ! displacements within receive buffer (flattened array)
        buf_dspls=(/(rlen*(firstofset(j,set)-1),j=0,procs-1)/)
        ! use recieve buffer as sendbuffer by specifying mpi_in_place
        if (ti) then
#ifdef MPI1
#define BUFFER bufi
          allocate(bufi(buf_n(rank+1)))
          bufi(:)=ibuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
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
        if (tr) then
#ifdef MPI1
#define BUFFER bufr
          allocate(bufr(buf_n(rank+1)))
          bufr(:)=rbuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
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
        if (trlp) then
#ifdef MPI1
#define BUFFER bufrlp
          allocate(bufrlp(buf_n(rank+1)))
          bufrlp(:)=rlpbuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
#endif
          call mpi_allgatherv(BUFFER, &
            buf_n(rank+1), &
            MPI_REAL4, &
            rlpbuf, &
            buf_n, &
            buf_dspls, &
            MPI_REAL4, &
            mpi_comm_world, &
            ierr)
#ifdef MPI1
          deallocate(bufrlp)
#undef BUFFER
#endif
        end if
        if (tz) then
#ifdef MPI1
#define BUFFER bufz
          allocate(bufz(buf_n(rank+1)))
          bufz(:)=zbuf(buf_dspls(rank+1)+1:buf_dspls(rank+1)+buf_n(rank+1))
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
        deallocate(buf_n,buf_dspls)
#endif
      end subroutine

End Module modmpi
