!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine tetgather
      Use modmain
      Use modxs
      Use modmpi
      Use m_gettetcw
      Use m_puttetcw
      Use m_filedel
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'tetgather'
      Character (256) :: filnam, filnam_t
      Integer :: iq, iproc, ik, i1, i2, nwdfp
      Real (8), Allocatable :: cw (:), cwa (:), cwsurf (:)
      Real (8), Allocatable :: cwp (:), cwap (:), cwsurfp (:)
      Allocate (cw(nwdf), cwa(nwdf), cwsurf(nwdf))
  ! loop over q-points
      Do iq = 1, nqpt
     ! calculate k+q and G+k+q related variables
         Call init1offs (qvkloff(1, iq))
     ! file name for output file
         Call genfilname (basename='TETW', iqmt=iq, filnam=filnam)
         Do ik = 1, nkpt
            Do i1 = 1, nst1
               Do i2 = 1, nst2
              ! collect weights from processes
                  Do iproc = 0, procs - 1
                 ! filename for input file
                     Call genfilname (basename='TETW', iqmt=iq, &
                    & rank=iproc, procs=procs, filnam=filnam_t)
                     wpari = firstofset (iproc, nwdf)
                     wparf = lastofset (iproc, nwdf)
                     nwdfp = wparf - wpari + 1
                     Allocate (cwp(nwdfp), cwap(nwdfp), cwsurfp(nwdfp))
                     Call gettetcw (iq, ik, i1, i2, nst1, nst2, nwdfp, &
                    & trim(filnam_t), cwp, cwap, cwsurfp)
                     cw (wpari:wparf) = cwp (:)
                     cwa (wpari:wparf) = cwap (:)
                     cwsurf (wpari:wparf) = cwsurfp (:)
                     Deallocate (cwp, cwap, cwsurfp)
                  End Do ! iproc
              ! write weights
                  Call puttetcw (iq, ik, i1, i2, nst1, nst2, &
                 & trim(filnam), cw, cwa, cwsurf)
               End Do
            End Do
        ! end loop over k-points
         End Do
         Do iproc = 0, procs - 1
            Call genfilname (basename='TETW', iqmt=iq, rank=rank, &
           & procs=procs, filnam=filnam_t)
            Call filedel (trim(filnam_t))
         End Do
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): weights f&
        &or tetrahedron method gathered for q - point:', iq
     ! end loop over q-points
      End Do
      Deallocate (cw, cwa, cwsurf)
End Subroutine tetgather
