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
      Integer :: ik
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
        call mpi_allgatherv_ifc(nkpt,ncrmax*natmtot*nstsv,zbuf=vnlcv)
        call mpi_allgatherv_ifc(nkpt,nstsv*nstsv,zbuf=vnlvv)
#endif
         Return
   End Subroutine
