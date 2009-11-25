!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepvnl (vnlcv, vnlvv)
      Use modmain
      Use modmpi
      Implicit None
! arguments
      Complex (8), Intent (Out) :: vnlcv (ncrmax, natmtot, nstsv, nkpt)
      Complex (8), Intent (Out) :: vnlvv (nstsv, nstsv, nkpt)
! local variables
      Integer :: ik, i
      Integer :: mpireccnts (procs), mpirecdispls (procs), mpisndcnts &
     & (procs), mpisnddispls (procs), kgatherdispls (procs), &
     & kgatherrecvcnts (procs)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
#ifdef MPIEXX
      Do ik = firstk (rank), lastk (rank)
         Write (*, '("Info(oepvnl): ", I6, " of ", I6, " k-points on pr&
        &oc:", I6)') ik, nkpt, rank
#endif
#ifndef MPIEXX
         Do ik = 1, nkpt
            Write (*, '("Info(oepvnl): ", I6, " of ", I6, " k-points")') ik, nkpt
#endif
            Call oepvnlk (ik, vnlcv(:, :, :, ik), vnlvv(:, :, ik))
         End Do
!$OMP END DO
!$OMP END PARALLEL
!
#ifdef MPIEXX
         Do i = 0, procs - 1
            kgatherdispls (i+1) = firstk (i) - 1
            kgatherrecvcnts (i+1) = nofk (i)
         End Do
         mpireccnts = kgatherrecvcnts * ncrmax * natmtot * nstsv
         mpirecdispls = kgatherdispls * ncrmax * natmtot * nstsv
         mpisndcnts (:) = nofk (rank) * ncrmax * natmtot * nstsv
         mpisnddispls (:) = (firstk(rank)-1) * ncrmax * natmtot * nstsv
         Call MPI_Alltoallv (vnlcv, mpisndcnts, mpisnddispls, &
        & MPI_DOUBLE_COMPLEX, vnlcv, mpireccnts, mpirecdispls, &
        & MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
         mpireccnts = kgatherrecvcnts * nstsv * nstsv
         mpirecdispls = kgatherdispls * nstsv * nstsv
         mpisndcnts (:) = nofk (rank) * nstsv * nstsv
         mpisnddispls (:) = (firstk(rank)-1) * nstsv * nstsv
         Call MPI_Alltoallv (vnlvv, mpisndcnts, mpisnddispls, &
        & MPI_DOUBLE_COMPLEX, vnlvv, mpireccnts, mpirecdispls, &
        & MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
#endif
!
         Return
   End Subroutine
