!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine terminate
      Use modmpi
      Implicit None
  ! abort MPI if necessary
#ifdef MPI
      Call mpi_abort (mpi_comm_world, 1, ierr)
      If (ierr .Eq. 0) Then
         Write (*, '(a)') 'MPI abort'
      Else
         Write (*, '(a)') 'MPI abort with errors - zombie processes mig&
        &ht remain!'
      End If
#endif
#ifndef MPI
      Write (*, '(a)') 'Abort'
#endif
  ! stop program
      Stop
End Subroutine terminate
