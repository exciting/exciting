! Copyright (C) 2006-2008 C. Ambrosch-Draxl. C. Meisenbichler S. Sagmeister
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!BOP
! !MODULE:  modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   In case of compiled without MPI support it defines the
!   mpi specific variables such that the code behaves exactly as 
!   the unmodified scalar version 
!#include "/usr/local/include/mpif.h"

! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (Sagmeister)
!EOP
module  modmpi
#ifdef MPI
#include "mpif.h"
!  use mpi
#endif
  integer :: rank
  integer :: procs
  integer :: ierr
  logical :: splittfile

!!$  character(256)::filename
contains
  subroutine initMPI
#ifdef MPI
    !        mpi init
    call mpi_init(ierr)
    call mpi_comm_size( mpi_comm_world, procs, ierr)
    call mpi_comm_rank( mpi_comm_world, rank, ierr)
    splittfile=.true.
#endif
#ifndef MPI
    procs=1
    rank=0
    splittfile=.false.
#endif
  end subroutine initMPI

  subroutine finitMPI
#ifdef MPI 
    call MPI_Finalize(ierr)
#endif
  end subroutine finitmpi

! service functions to partition k points
! still kept around but should be replaced by generig functions
! firstofset nofset lastofset ....
function nofk(process)
use modmain, only:nkpt
 integer::nofk
 integer, intent(in)::process
 nofk=nkpt/procs
 if ((mod(nkpt,procs).gt.process)) nofk=nofk+1
end function nofk

function firstk(process)
integer::firstk
integer, intent(in)::process
firstk=1
do i=0,process-1
firstk=firstk+nofk(i)
end do
end function firstk

function lastk(process)
integer::lastk
integer, intent(in)::process
lastk=0
do i=0,process
lastk=lastk+nofk(i)
end do
end function lastk

function procofk(k)
   integer:: procofk
   integer, intent(in)::k
   integer::iproc
   procofk=0
   do iproc=0,procs-1
      if (k.gt.lastk(iproc)) procofk=procofk+1
   end do
end function procofk

!-------------generalized partition!
function nofset(process,set)
 integer::nofset
 integer, intent(in)::process,set
 nofset=set/procs
 if ((mod(set,procs).gt.process)) nofset=nofset+1
end function nofset

function firstofset(process,set)
integer::firstofset
integer, intent(in)::process,set
firstofset=1
do i=0,process-1
firstofset=firstofset+nofset(i,set)
end do
end function firstofset

function lastofset(process,set)
integer::lastofset
integer, intent(in)::process,set
lastofset=0
do i=0,process
lastofset=lastofset+nofset(i,set)
end do
end function lastofset

function procofindex(k,set)
   integer:: procofindex
   integer, intent(in)::k,set
   integer::iproc
  procofindex=0
   do iproc=0,procs-1
      if (k.gt.lastofset(iproc,set)) procofindex=procofindex+1
   end do
end function procofindex

function lastproc(row,set)
  implicit none
  integer :: lastproc
  integer, intent(in) :: row,set
  integer :: iproc
  if (row.ne.nofset(0,set)) then
     lastproc=procs
  else
     lastproc=modulo(set,procs)
     if (lastproc.eq.0) lastproc=procs
  end if
  lastproc=lastproc-1
end function lastproc



!------------------interface to MPI_barrier for xs-part
subroutine barrier
  implicit none
  ! do nothing if only one process
  if (procs.eq.1) return
  ! call the MPI barrier
#ifdef MPI
  call MPI_barrier(mpi_comm_world,ierr)

write(300+rank,*) 'barrier, rank=',rank
call flushifc(300+rank)

#endif
end subroutine barrier

subroutine endloopbarrier(set,mult)
  implicit none
  integer, intent(in) :: set,mult
  integer :: i,im
  do i=1,(nofset(0,set)-nofset(rank,set))*mult
     call barrier
  end do
end subroutine endloopbarrier

!------------------wrappers for MPI communication
subroutine zalltoallv(zarr,rlen,set)
  implicit none
  ! arguments
  complex(8), intent(inout) :: zarr(*)
  integer, intent(in) :: rlen,set
#ifdef MPI
  ! local variables
  integer :: mpireccnts(procs),mpirecdispls(procs)
  integer :: mpisndcnts(procs),mpisnddispls(procs)
  integer :: proc
  mpisndcnts(:)=nofset(rank,set)*rlen
  mpisnddispls(:)=(firstofset(rank,set)-1)*rlen
  do proc=0,procs-1
     mpireccnts(proc+1)=nofset(proc,set)*rlen
     mpirecdispls(proc+1)=(firstofset(proc,set)-1)*rlen
  end do
  call MPI_Alltoallv(zarr,mpisndcnts,mpisnddispls,MPI_DOUBLE_COMPLEX, &
       zarr,mpireccnts,mpirecdispls,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
#endif
end subroutine zalltoallv

subroutine ralltoallv(rarr,rlen,set)
  implicit none
  ! arguments
  real(8), intent(inout) :: rarr(*)
  integer, intent(in) :: rlen,set
#ifdef MPI
  ! local variables
  integer :: mpireccnts(procs),mpirecdispls(procs)
  integer :: mpisndcnts(procs),mpisnddispls(procs)
  integer :: proc
  mpisndcnts(:)=nofset(rank,set)*rlen
  mpisnddispls(:)=(firstofset(rank,set)-1)*rlen
  do proc=0,procs-1
     mpireccnts(proc+1)=nofset(proc,set)*rlen
     mpirecdispls(proc+1)=(firstofset(proc,set)-1)*rlen
  end do
  call MPI_Alltoallv(rarr,mpisndcnts,mpisnddispls,MPI_DOUBLE_COMPLEX, &
       rarr,mpireccnts,mpirecdispls,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
#endif
end subroutine ralltoallv

end module modmpi
