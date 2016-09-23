!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putbsemat (fname, zmat, ikkp, iknr, jknr, iq, iqr, n1, n2, &
& n3, n4)
      Use modmain
      Use modmpi
      Use m_getunit
      Implicit None
  ! arguments
      Character (*), Intent (In) :: fname
      Complex (8), Intent (In) :: zmat (n1, n2, n3, n4)
      Integer, Intent (In) :: ikkp, iknr, jknr, iq, iqr, n1, n2, n3, n4
  ! local variables
      Integer :: recl, un, iknrr, jknrr, ikkpr, iq_r, iqr_r
#ifdef MPI
      Integer :: nkkp, iproc, tag, status (MPI_STATUS_SIZE)
#endif
      ikkpr = ikkp
      iknrr = iknr
      jknrr = jknr
      iq_r = 0
      iqr_r = 0
      Call getunit (un)
      Inquire (IoLength=Recl) ikkp, iknr, jknr, iq, iqr, n1, n2, n3, &
     & n4, zmat
#ifdef MPI
      tag = 77
      nkkp = (nkptnr*(nkptnr+1)) / 2
      If (rank .Ne. 0) Call mpi_send (zmat, size(zmat), &
     & MPI_DOUBLE_COMPLEX, 0, tag, MPI_COMM_WORLD, ierr)
      If (rank .Eq. 0) Then
         Do iproc = 0, lastproc (ikkp, nkkp)
            ikkpr = firstofset (iproc, nkkp) - 1 + ikkp
            Call kkpmap (ikkpr, nkptnr, iknrr, jknrr)
!
            If (iproc .Ne. 0) Then
           ! receive data from slaves
               Call mpi_recv (zmat, size(zmat), MPI_DOUBLE_COMPLEX, &
              & iproc, tag, MPI_COMM_WORLD, status, ierr)
            End If
#endif
        ! only the master is performing file I/O
            Open (Unit=un, File=trim(fname), Form='unformatted', &
           & Action='write', Access='direct', Recl=Recl)
            Write (un, Rec=ikkpr) ikkpr, iknrr, jknrr, iq_r, iqr_r, n1, &
           & n2, n3, n4, zmat
            Close (un)
#ifdef MPI
         End Do
      End If
#endif
End Subroutine putbsemat
