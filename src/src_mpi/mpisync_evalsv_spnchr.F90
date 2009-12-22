!
!
! Copyright (C) 2006 C. Ambrosch-Draxl. C. Meisenbichler
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!BOP
! !ROUTINE: mpigather
! !INTERFACE:
!
!
Subroutine mpisync_evalsv_spnchr
#ifdef MPI
  ! !USES:
      Use modmpi
      Use modmain
  ! !DESCRIPTION:
  !  Routine redistributes global vectors after k-parallel computation
  !
  ! !REVISION HISTORY:
  !   Created October SEPT 2006 (MULEOBEN)
  !   by Cristian Meisenbichler
  !EOP
      Implicit None
  !local variables
      Integer :: mpireccnts (procs), mpirecdispls (procs), mpisndcnts &
     & (procs), mpisnddispls (procs)
      Integer :: i, kgatherdispls (procs), kgatherrecvcnts (procs)
      Do i = 0, procs - 1
         kgatherdispls (i+1) = firstk (i) - 1
         kgatherrecvcnts (i+1) = nofk (i)
      End Do
      mpireccnts = kgatherrecvcnts * nstsv
      mpirecdispls = kgatherdispls * nstsv
      mpisndcnts = nofk (rank) * nstsv
      mpisnddispls = (firstk(rank)-1) * nstsv
      Call MPI_Alltoallv (evalsv, mpisndcnts, mpisnddispls, &
     & MPI_DOUBLE_PRECISION, evalsv, mpireccnts, mpirecdispls, &
     & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
!  mpireccnts=kgatherrecvcnts*nspinor*nstsv
!  mpirecdispls=kgatherdispls*nspinor*nstsv
!  mpisndcnts=nofk(rank)*nspinor*nstsv
!  mpisnddispls=(firstk(rank)-1)*nspinor*nstsv
!  call MPI_Alltoallv(spnchr, mpisndcnts, mpisnddispls, MPI_DOUBLE_PRECISION,&
!       spnchr, mpireccnts,  mpirecdispls,  MPI_DOUBLE_PRECISION ,  !MPI_COMM_WORLD,ierr)
!
#endif
!
End Subroutine mpisync_evalsv_spnchr
