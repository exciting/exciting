!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_putemat
      Implicit None
Contains
!
!
      Subroutine putemat (iq, ik, tarec, filnam, l1, h1, l2, h2, x12, &
     & l3, h3, l4, h4, x34)
         Use modmain
         Use modmpi
         Use modxs
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik
         Logical, Intent (In) :: tarec
         Character (*), Intent (In) :: filnam
         Integer, Intent (In) :: l1, h1, l2, h2
         Complex (8), Intent (InOut) :: x12 (:, :, :)
         Integer, Optional, Intent (In) :: l3, h3, l4, h4
         Complex (8), Optional, Intent (InOut) :: x34 (:, :, :)
    ! local variables
         Integer :: un, recl, ikr
         Logical :: tarec_
#ifdef MPI
         Integer :: iproc, tag1, tag2, status (MPI_STATUS_SIZE)
#endif
    ! check if all optional variables are present if any is present
         If ((present(l3) .Or. present(h3) .Or. present(l4) .Or. &
        & present(h4) .Or. present(x34)) .And. ( .Not. (present(l3) &
        & .And. present(h3) .And. present(l4) .And. present(h4) .And. &
        & present(x34)))) Then
            Write (*,*)
            Write (*, '("Error(putemat): optional parameters not comple&
           &te - check routine")')
            Write (*,*)
            Call terminate
         End If
         tarec_ = tarec
         ikr = ik
         Call getunit (un)
         If (present(x34)) Then
#ifdef MPI
            tag1 = 77
            tag2 = 78
            If (rank .Ne. 0) Call mpi_send (x12, size(x12), &
           & MPI_DOUBLE_COMPLEX, 0, tag1, MPI_COMM_WORLD, ierr)
            If (rank .Ne. 0) Call mpi_send (x34, size(x34), &
           & MPI_DOUBLE_COMPLEX, 0, tag2, MPI_COMM_WORLD, ierr)
            If (rank .Eq. 0) Then
               Do iproc = 0, lastproc (ik, nkpt)
                  ikr = firstofset (iproc, nkpt) - 1 + ik
                  If (iproc .Ne. 0) Then
                ! receive data from slaves
                     Call mpi_recv (x12, size(x12), MPI_DOUBLE_COMPLEX, &
                    & iproc, tag1, MPI_COMM_WORLD, status, ierr)
                     Call mpi_recv (x34, size(x34), MPI_DOUBLE_COMPLEX, &
                    & iproc, tag2, MPI_COMM_WORLD, status, ierr)
                  End If
#endif
             ! I/O record length
                  Inquire (IoLength=Recl) vql (:, iq), vkl (:, ikr), &
                 & nstsv, ngq (iq), l1, h1, l2, h2, l3, h3, l4, h4, &
                 & x12, x34
                  Open (Unit=un, File=trim(filnam), Form='unformatted', &
                 & Action='write', Access='direct', Recl=Recl)
                  Write (un, Rec=ikr) vql (:, iq), vkl (:, ikr), nstsv, &
                 & ngq (iq), l1, h1, l2, h2, l3, h3, l4, h4, x12, x34
#ifdef MPI
               End Do
            End If
#endif
         Else
#ifdef MPI
            tag1 = 77
            If (rank .Ne. 0) Call mpi_send (x12, size(x12), &
           & MPI_DOUBLE_COMPLEX, 0, tag1, MPI_COMM_WORLD, ierr)
            If (rank .Eq. 0) Then
               Do iproc = 0, lastproc (ik, nkpt)
                  ikr = firstofset (iproc, nkpt) - 1 + ik
                  If (iproc .Ne. 0) Then
                ! receive data from slaves
                     Call mpi_recv (x12, size(x12), MPI_DOUBLE_COMPLEX, &
                    & iproc, tag1, MPI_COMM_WORLD, status, ierr)
                  End If
#endif
             ! I/O record length
                  Inquire (IoLength=Recl) vql (:, iq), vkl (:, ikr), &
                 & nstsv, ngq (iq), l1, h1, l2, h2, x12
                 write(*,*) "write k",ikr, "in file",trim(filnam)
                 write(*,*)"recl",recl
                 write(*,*)"shape(l1)",shape(l1), "shape(h1)",shape(h1)
                 write(*,*)"shape(l2)", shape(l2),"shape(h2)",shape(h2)
                 write(*,*)"shape(x12)",shape(x12)

                  Open (Unit=un, File=trim(filnam), Form='unformatted', &
                 & Action='write', Access='direct', Recl=Recl)
                  Write (un, Rec=ikr) vql (:, iq), vkl (:, ikr), nstsv, &
                 & ngq (iq), l1, h1, l2, h2, x12
#ifdef MPI
               End Do
            End If
#endif
         End If
         Close (un)
      End Subroutine putemat
!
End Module m_putemat
