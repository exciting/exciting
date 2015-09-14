!
!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dfgather
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use m_filedel
      Use m_getx0
      Use m_putx0
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'dfgather'
      Integer :: n, iq, iw, iproc
      Complex (8), Allocatable :: chi0 (:, :), chi0wg (:, :, :), chi0hd &
     & (:, :)
      Logical :: tq0
      Logical, External :: tqgamma
  ! loop over q-points
      Do iq = 1, nqpt
         tq0 = tqgamma (iq)
     ! calculate k+q and G+k+q related variables
         Call init1offs (qvkloff(1, iq))
     ! size of local field effects
         n = ngq (iq)
     ! allocate
         Allocate (chi0(n, n), chi0wg(n, 2, 3), chi0hd(3, 3))
     ! file extension for q-point
         Do iproc = 0, procs - 1
            Call genfilname (basename='X0', bzsampl=bzsampl, &
           & acont=input%xs%tddft%acont, nar= .Not. &
           & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, &
           & markfxcbse=tfxcbse, iqmt=iq, procs=procs, rank=iproc, &
           & filnam=fnchi0_t)
            wpari = firstofset (iproc, nwdf)
            wparf = lastofset (iproc, nwdf)
            Do iw = wpari, wparf
           ! exponential factor matrix elements
               if(iproc.eq.rank)then
               Call getx0 (tq0, iq, iw-wpari+1, trim(fnchi0_t), '', &
              & chi0, chi0wg, chi0hd)
               endif
#ifdef MPI
               call mpi_bcast(chi0,n*n,MPI_DOUBLE_COMPLEX,iproc,mpi_comm_world,ierr)
               call mpi_bcast(chi0wg,n*2*3,MPI_DOUBLE_COMPLEX,iproc,mpi_comm_world,ierr)
               call mpi_bcast(chi0hd,3*3,MPI_DOUBLE_COMPLEX,iproc,mpi_comm_world,ierr)
#endif
               if(rank.eq.0.or.(firstinnode.and. .not.input%sharedfs))then
               Call putx0 (tq0, iq, iw, trim(fnchi0), '', chi0, chi0wg, &
              & chi0hd)
              endif
            End Do
         End Do
         Call genfilname (basename='X0', iqmt=iq, procs=procs, &
           & rank=rank, filnam=fnchi0_t)
         Deallocate (chi0, chi0wg, chi0hd)
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): Kohn Sham&
        & response function gathered for q - point:', iq
      End Do
End Subroutine dfgather
