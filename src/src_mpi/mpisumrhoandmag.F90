
! Copyright (C) 2006 C. Ambrosch-Draxl. C. Meisenbichler
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!BOP
! !ROUTINE: mpigather
! !INTERFACE:


subroutine mpisumrhoandmag
#ifdef MPI
use modinput
  ! !USES:
  use modmpi
  use modmain
  ! !DESCRIPTION:
  !  Subroutine Builds the density sums after k-parallel rhovalk 
  !  and broadcasts it to all nodes
  !
  ! !REVISION HISTORY:
  !   Created October SEPT 2006 (MULEOBEN)
  !   by Cristian Meisenbichler
  !EOP
  implicit none
  real(8), allocatable :: buffer(:)
  real(8), allocatable :: buffer2d(:, :)
  real(8), allocatable :: buffer3d(:, :, :)
  real(8), allocatable :: buffer4d(:, :, :, :)
  allocate(buffer3d(lmmaxvr, nrmtmax, natmtot))
  buffer3d=0.d0
  call MPI_barrier(MPI_COMM_WORLD, ierr)
  call MPI_allreduce(rhomt, buffer3d, lmmaxvr * nrmtmax * natmtot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  rhomt=buffer3d
  deallocate(buffer3d)
  allocate(buffer(ngrtot))
  buffer=0.d0
  call MPI_barrier(MPI_COMM_WORLD, ierr)
  call MPI_allreduce(rhoir, buffer, ngrtot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  rhoir=buffer
  deallocate(buffer)
  if(associated(input%groundstate%spin)) then
     allocate(buffer4d(lmmaxvr, nrmtmax, natmtot, ndmag))
     buffer4d=0.d0
     call MPI_barrier(MPI_COMM_WORLD, ierr)
     call MPI_allreduce(magmt, buffer4d, lmmaxvr * nrmtmax * natmtot * ndmag, MPI_DOUBLE_PRECISION, MPI_SUM,&
    &MPI_COMM_WORLD, ierr)
     magmt=buffer4d
     deallocate(buffer4d)
     allocate(buffer2d(ngrtot, ndmag))
     buffer2d=0.d0
     call MPI_barrier(MPI_COMM_WORLD, ierr)
     call MPI_allreduce(magir, buffer2d, ngrtot * ndmag, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     magir=buffer2d
     deallocate(buffer2d)
  endif
if(rank.eq.0) write(60, *) "MPIINFO parallel rho mpireduced"
#endif
end subroutine mpisumrhoandmag
