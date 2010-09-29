!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_putpmat
      Implicit None
Contains
!
!
      Subroutine putpmat (ik, tarec, filnam, pm)
         Use modmain
         Use modmpi
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: ik
    ! true if absolut record position is ik
         Logical, Intent (In) :: tarec
         Character (*), Intent (In) :: filnam
         Complex (8), Intent (In) :: pm (:, :, :)
    ! local variables
         Integer :: un, recl, ikr
         Logical :: tarect
#ifdef MPI
         Integer :: iproc, tag, status (MPI_STATUS_SIZE)
#endif
         tarect = tarec
         ikr = ik
         Inquire (IoLength=Recl) vkl (:, ik), nstsv, pm
         Call getunit (un)
#ifdef MPI
         tag = 77
         If (rank .Ne. 0) Call mpi_send (pm, size(pm), &
        & MPI_DOUBLE_COMPLEX, 0, tag, MPI_COMM_WORLD, ierr)
         If (rank .Eq. 0) Then
            Do iproc = 0, lastproc (ik, nkpt)
               ikr = firstofset (iproc, nkpt) - 1 + ik
               If (iproc .Ne. 0) Then
             ! receive data from slaves
                  Call mpi_recv (pm, size(pm), MPI_DOUBLE_COMPLEX, &
                 & iproc, tag, MPI_COMM_WORLD, status, ierr)
               End If
#endif
          ! only master is performing I/O
               Open (Unit=un, File=trim(filnam), Form='unformatted', &
              & Action='write', Access='direct', Recl=Recl)
               Write (un, Rec=ikr) vkl (:, ikr), nstsv, pm
               Close (un)
#ifdef MPI
            End Do
         End If
#endif
      End Subroutine putpmat
!
End Module m_putpmat
