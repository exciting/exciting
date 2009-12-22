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
Subroutine mpisumrhoandmag
#ifdef MPI
      Use modinput
  ! !USES:
      Use modmpi
      Use modmain
  ! !DESCRIPTION:
  !  Subroutine Builds the density sums after k-parallel rhovalk
  !  and broadcasts it to all nodes
  !
  ! !REVISION HISTORY:
  !   Created October SEPT 2006 (MULEOBEN)
  !   by Cristian Meisenbichler
  !EOP
      Implicit None
      Real (8), Allocatable :: buffer (:)
      Real (8), Allocatable :: buffer2d (:, :)
      Real (8), Allocatable :: buffer3d (:, :, :)
      Real (8), Allocatable :: buffer4d (:, :, :, :)
      Allocate (buffer3d(lmmaxvr, nrmtmax, natmtot))
      buffer3d = 0.d0
      Call MPI_barrier (MPI_COMM_WORLD, ierr)
      Call MPI_allreduce (rhomt, buffer3d, lmmaxvr*nrmtmax*natmtot, &
     & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      rhomt = buffer3d
      Deallocate (buffer3d)
      Allocate (buffer(ngrtot))
      buffer = 0.d0
      Call MPI_barrier (MPI_COMM_WORLD, ierr)
      Call MPI_allreduce (rhoir, buffer, ngrtot, MPI_DOUBLE_PRECISION, &
     & MPI_SUM, MPI_COMM_WORLD, ierr)
      rhoir = buffer
      Deallocate (buffer)
      If (associated(input%groundstate%spin)) Then
         Allocate (buffer4d(lmmaxvr, nrmtmax, natmtot, ndmag))
         buffer4d = 0.d0
         Call MPI_barrier (MPI_COMM_WORLD, ierr)
         Call MPI_allreduce (magmt, buffer4d, &
        & lmmaxvr*nrmtmax*natmtot*ndmag, MPI_DOUBLE_PRECISION, MPI_SUM, &
        & MPI_COMM_WORLD, ierr)
         magmt = buffer4d
         Deallocate (buffer4d)
         Allocate (buffer2d(ngrtot, ndmag))
         buffer2d = 0.d0
         Call MPI_barrier (MPI_COMM_WORLD, ierr)
         Call MPI_allreduce (magir, buffer2d, ngrtot*ndmag, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         magir = buffer2d
         Deallocate (buffer2d)
      End If
      If (rank .Eq. 0) write (60,*) "MPIINFO parallel rho mpireduced"
#endif
End Subroutine mpisumrhoandmag
