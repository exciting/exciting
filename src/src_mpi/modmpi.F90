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
!#include "/usr/local/include/mpif.h"
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (Sagmeister)
!
!
!
Module modmpi
#ifdef MPI
#include "../build/mpiconf.inc"
#endif
      Integer :: rank
      Integer :: procs
      Integer :: ierr
      Logical :: splittfile
!
!!$  character(256)::filename
Contains
      Subroutine initMPI
#ifdef MPI
    !        mpi init
         Call mpi_init (ierr)
         Call mpi_comm_size (mpi_comm_world, procs, ierr)
         Call mpi_comm_rank (mpi_comm_world, rank, ierr)
         splittfile = .True.
#endif
#ifndef MPI
         procs = 1
         rank = 0
         splittfile = .False.
#endif
      End Subroutine initMPI
!
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
         firstk = 1
         Do i = 0, process - 1
            firstk = firstk + nofk (i)
         End Do
      End Function firstk
!
      Function lastk (process)
         Integer :: lastk
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
         Integer :: nofset
         Integer, Intent (In) :: process, set
         nofset = set / procs
         If ((Mod(set, procs) .Gt. process)) nofset = nofset + 1
      End Function nofset
!
      Function firstofset (process, set)
         Integer :: firstofset
         Integer, Intent (In) :: process, set
         firstofset = 1
         Do i = 0, process - 1
            firstofset = firstofset + nofset (i, set)
         End Do
      End Function firstofset
!
      Function lastofset (process, set)
         Integer :: lastofset
         Integer, Intent (In) :: process, set
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
         Integer :: iproc
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
         If (procs .Eq. 1) Return
  ! call the MPI barrier
#ifdef MPI
         Call MPI_barrier (mpi_comm_world, ierr)
#endif
      End Subroutine barrier
!
!
      Subroutine endloopbarrier (set, mult)
         Implicit None
         Integer, Intent (In) :: set, mult
         Integer :: i, im
         Do i = 1, (nofset(0, set)-nofset(rank, set)) * mult
            Call barrier
         End Do
      End Subroutine endloopbarrier
!
!------------------wrappers for MPI communication
!
!
      Subroutine zalltoallv (zarr, rlen, set)
         Implicit None
  ! arguments
         Complex (8), Intent (Inout) :: zarr (*)
         Integer, Intent (In) :: rlen, set
#ifdef MPI
  ! local variables
         Integer :: mpireccnts (procs), mpirecdispls (procs)
         Integer :: mpisndcnts (procs), mpisnddispls (procs)
         Integer :: proc
         mpisndcnts (:) = nofset (rank, set) * rlen
         mpisnddispls (:) = (firstofset(rank, set)-1) * rlen
         Do proc = 0, procs - 1
            mpireccnts (proc+1) = nofset (proc, set) * rlen
            mpirecdispls (proc+1) = (firstofset(proc, set)-1) * rlen
         End Do
         Call MPI_Alltoallv (zarr, mpisndcnts, mpisnddispls, &
        & MPI_DOUBLE_COMPLEX, zarr, mpireccnts, mpirecdispls, &
        & MPI_DOUBLE_COMPLEX, mpi_comm_world, ierr)
#endif
      End Subroutine zalltoallv
!
!
      Subroutine ralltoallv (rarr, rlen, set)
         Implicit None
  ! arguments
         Real (8), Intent (Inout) :: rarr (*)
         Integer, Intent (In) :: rlen, set
#ifdef MPI
  ! local variables
         Integer :: mpireccnts (procs), mpirecdispls (procs)
         Integer :: mpisndcnts (procs), mpisnddispls (procs)
         Integer :: proc
         mpisndcnts (:) = nofset (rank, set) * rlen
         mpisnddispls (:) = (firstofset(rank, set)-1) * rlen
         Do proc = 0, procs - 1
            mpireccnts (proc+1) = nofset (proc, set) * rlen
            mpirecdispls (proc+1) = (firstofset(proc, set)-1) * rlen
         End Do
         Call MPI_Alltoallv (rarr, mpisndcnts, mpisnddispls, &
        & MPI_DOUBLE_PRECISION, rarr, mpireccnts, mpirecdispls, &
        & MPI_DOUBLE_PRECISION, mpi_comm_world, ierr)
#endif
      End Subroutine ralltoallv
!
End Module modmpi
