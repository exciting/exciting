
! Copyright (C) 2006 C. Ambrosch-Draxl. C. Meisenbichler
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details

!BOP
! !ROUTINE: mpisumrhoandmag
! !INTERFACE:
Subroutine mpisumrhoandmag
#ifdef MPI
      Use modinput
! !USES:
      Use modmpi
      Use modmain
! !DESCRIPTION:
!  This subroutine adds up the partial sums for the charge density, the magnetisation,
!  and the partial charges for all processors and broadcasts it to all processors.
!
! !REVISION HISTORY:
!   Created October September 2006 (Christian Meisenbichler)
!   Modifications, August 2010 (Stephan Sagmeister)
!EOP
!BOC
      Implicit None
#ifdef MPI1
      Real (8), Allocatable :: buffer (:)
      Real (8), Allocatable :: buffer2d (:, :)
      Real (8), Allocatable :: buffer3d (:, :, :)
      Real (8), Allocatable :: buffer4d (:, :, :, :)
#endif
!------------------------!
!   muffin-tin density   !
!------------------------!
      Call MPI_barrier (MPI_COMM_WORLD, ierr)
#ifndef MPI1
      Call MPI_allreduce (mpi_in_place, rhomt, lmmaxvr*nrmtmax*natmtot, &
     & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#ifdef MPI1
      Allocate (buffer3d(lmmaxvr, nrmtmax, natmtot))
      buffer3d = 0.d0
      Call MPI_allreduce (rhomt, buffer3d, lmmaxvr*nrmtmax*natmtot, &
     & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      rhomt = buffer3d
      Deallocate (buffer3d)
#endif
!--------------------------!
!   interstitial density   !
!--------------------------!
      Call MPI_barrier (MPI_COMM_WORLD, ierr)
#ifndef MPI1
      Call MPI_allreduce (mpi_in_place, rhoir, ngrtot, MPI_DOUBLE_PRECISION, &
     & MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#ifdef MPI1
      Allocate (buffer(ngrtot))
      buffer = 0.d0
      Call MPI_allreduce (rhoir, buffer, ngrtot, MPI_DOUBLE_PRECISION, &
     & MPI_SUM, MPI_COMM_WORLD, ierr)
      rhoir = buffer
      Deallocate (buffer)
#endif
      If (associated(input%groundstate%spin)) Then
!------------------------------!
!   muffin-tin magnetization   !
!------------------------------!
         Call MPI_barrier (MPI_COMM_WORLD, ierr)
#ifndef MPI1
         Call MPI_allreduce (mpi_in_place, magmt, &
        & lmmaxvr*nrmtmax*natmtot*ndmag, MPI_DOUBLE_PRECISION, MPI_SUM, &
        & MPI_COMM_WORLD, ierr)
#endif
#ifdef MPI1
         Allocate (buffer4d(lmmaxvr, nrmtmax, natmtot, ndmag))
         buffer4d = 0.d0
         Call MPI_allreduce (magmt, buffer4d, &
        & lmmaxvr*nrmtmax*natmtot*ndmag, MPI_DOUBLE_PRECISION, MPI_SUM, &
        & MPI_COMM_WORLD, ierr)
         magmt = buffer4d
         Deallocate (buffer4d)
#endif
!--------------------------------!
!   interstitial magnetization   !
!--------------------------------!
         Call MPI_barrier (MPI_COMM_WORLD, ierr)
#ifndef MPI1
         Call MPI_allreduce (mpi_in_place, magir, ngrtot*ndmag, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#ifdef MPI1
         Allocate (buffer2d(ngrtot, ndmag))
         buffer2d = 0.d0
         Call MPI_allreduce (magir, buffer2d, ngrtot*ndmag, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         magir = buffer2d
         Deallocate (buffer2d)
#endif
      End If
!--------------------------------!
!   muffin-tin partial charges   !
!--------------------------------!
      Call MPI_barrier (MPI_COMM_WORLD, ierr)
#ifndef MPI1
      call mpi_allreduce(mpi_in_place,chgpart,lmmaxvr*natmtot*nstsv, &
      & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#ifdef MPI1
      Allocate (buffer3d(lmmaxvr,natmtot,nstsv))
      buffer3d = 0.d0
      call mpi_allreduce(chgpart,buffer3d,lmmaxvr*natmtot*nstsv, &
      & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      chgpart = buffer3d
      Deallocate (buffer3d)
#endif
      If (rank .Eq. 0) write (60,*) " MPI: summation (MPI_Allreduce) over processes for density, "// &
       & "magnetisation, and partial charges done"
#endif
End Subroutine mpisumrhoandmag
!EOC
