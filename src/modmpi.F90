! Copyright (C) 2006 C. Ambrosch-Draxl. C. Meisenbichler
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
!EOP
module  modmpi
#ifdef MPI
  use mpi
#endif
  integer :: rank
  integer :: procs
  integer :: ierr
  logical :: splittfile

  character(256)::filename
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

end module modmpi
