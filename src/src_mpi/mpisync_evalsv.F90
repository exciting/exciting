
! Copyright (C) 2006-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
!BOP
! !ROUTINE: mpisync_evalsv
! !INTERFACE:
Subroutine mpisync_evalsv
#ifdef MPI
  ! !USES:
      Use modmpi
      Use modmain
  ! !DESCRIPTION:
  !  Routine redistributes global vectors after k-parallel computation
  !
  ! !REVISION HISTORY:
  !   Created October 2006 (Christian Meisenbichler)
  !EOP
  !BOC
      Implicit None
  !local variables
      Integer :: mpireccnts (procs), mpirecdispls (procs)
      Integer :: i, kgatherdispls (procs), kgatherrecvcnts (procs)
      real(8), allocatable :: buf(:)
      allocate(buf(nofk(rank)*nstsv))
      Do i = 0, procs - 1
         kgatherdispls (i+1) = firstk (i) - 1
         kgatherrecvcnts (i+1) = nofk (i)
      End Do
      mpireccnts = kgatherrecvcnts * nstsv
      mpirecdispls = kgatherdispls * nstsv
      buf=reshape(evalsv(:,firstk(rank):lastk(rank)),(/nstsv*nofk(rank)/))
      call MPI_Allgatherv(buf,mpireccnts(rank+1),MPI_DOUBLE_PRECISION, &
        evalsv,mpireccnts,mpirecdispls,MPI_DOUBLE_PRECISION, &
        MPI_COMM_WORLD,ierr)
      deallocate(buf)
#endif
End Subroutine
!EOC
