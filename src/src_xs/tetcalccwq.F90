!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine tetcalccwq (iq)
      Use modmain
      Use modinput
      Use modxs
      Use modtetra
      Use modmpi
      Use m_genwgrid
      Use m_puttetcw
      Use m_getunit
      Use m_filedel
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Character (*), Parameter :: thisnam = 'tetcalccwq'
      Character (256) :: filnam, filnamt
      Complex (8), Allocatable :: w (:)
      Real (8), Parameter :: epstetra = 1.d-8
      Real (8), Allocatable :: eb (:, :)
      Real (8), Allocatable :: wreal (:)
      Real (8), Allocatable :: cwsurft2 (:, :), cwt2 (:, :), cwat2 (:, &
     & :)
      Real (8), Allocatable :: cwsurft1 (:), cwt1 (:), cwat1 (:)
      Real (8), Allocatable :: cwsurf (:, :, :), cw (:, :, :), cwa (:, &
     & :, :)
      Real (8) :: wt
      Integer :: ik, ist1, ist2
      Integer :: iw, wi, wf, nwdfp, un, recl, irec
  ! calculate k+q and G+k+q related variables
      Call init1offs (qvkloff(1, iq))
  ! generate link array for tetrahedra
      Call gentetlinkp (vql(1, iq), input%xs%tetra%qweights)
  ! initial and final w-point
      wi = wpari
      wf = wparf
      nwdfp = wf - wi + 1
      If (tscreen) Then
     ! generate filenames
         Call genfilname (basename='TETW', iq=iq, rank=rank, &
        & procs=procs, appfilext=.True., filnam=filnam)
         Call genfilname (basename='TETWT', iq=iq, rank=rank, &
        & procs=procs, appfilext=.True., filnam=filnamt)
      Else
     ! set q-dependent file extension
         Call genfilname (iqmt=iq, setfilext=.True.)
     ! generate filenames
         Call genfilname (basename='TETW', iqmt=iq, rank=rank, &
        & procs=procs, filnam=filnam)
         Call genfilname (basename='TETWT', iqmt=iq, rank=rank, &
        & procs=procs, filnam=filnamt)
      End If
  ! find highest (partially) occupied and lowest (partially) unoccupied states
      Call findocclims (iq, istocc0, istocc, istunocc0, istunocc, &
     & isto0, isto, istu0, istu)
  ! find band combinations
      Call ematbdcmbs (input%xs%emattype)
  ! allocate arrays
      Allocate (eb(nstsv, nkpt))
      Allocate (cw(nstsv, nstsv, nkpt))
      Allocate (cwa(nstsv, nstsv, nkpt))
      Allocate (cwsurf(nstsv, nstsv, nkpt))
      Allocate (cwt2(nst1, nst2), cwat2(nst1, nst2), cwsurft2(nst1, &
     & nst2))
      Allocate (w(nwdf))
      Allocate (wreal(nwdfp))
  ! get the eigenvalues from file
      Do ik = 1, nkpt
         Call getevalsv (vkl(1, ik), evalsv(1, ik))
      End Do
      eb (:, :) = evalsv (:, :)
  ! scissors shift
      Where (eb .Gt. efermi) eb = eb + input%xs%scissor
  ! generate complex energy grid
      Call genwgrid (nwdf, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_cmplx=w)
      wreal (:) = dble (w(wi:wf))
      If (wreal(1) .Lt. epstetra) wreal (1) = epstetra
      Call getunit (un)
      Inquire (IoLength=Recl) cwt2, cwat2, cwsurft2
  ! open temporary file for writing
      Open (un, File=trim(filnamt), Form='unformatted', Action='write', &
     & Status='replace', Access='direct', Recl=Recl)
  ! calculate weights
      Do iw = 1, nwdfp
         If ((modulo(iw, Max(nwdfp/10, 1)) .Eq. 0) .Or. (iw .Eq. &
        & nwdfp)) write (*, '("Info(tetcalccwq): tetrahedron weights fo&
        &r ", I6, " of ", I6, " w-points")') iw, nwdfp
         wt = wreal (iw)
     ! switch 2 below in tetcw defines bulk integration for real part
     ! resonant contribution
         Call tetcwifc (nkpt, nstsv, eb, efermi, wt, 2, cw)
     ! anti-resonant contribution
         Call tetcwifc (nkpt, nstsv, eb, efermi,-wt, 2, cwa)
     ! switch 4 below in tetcw defines surface integration for imag. part
         Call tetcwifc (nkpt, nstsv, eb, efermi, wt, 4, cwsurf)
         Do ik = 1, nkpt
            irec = (ik-1) * nwdfp + iw
            cwsurft2 (:, :) = cwsurf (istl1:istu1, istl2:istu2, ik)
            cwt2 (:, :) = cw (istl1:istu1, istl2:istu2, ik)
            cwat2 (:, :) = cwa (istl1:istu1, istl2:istu2, ik)
            Write (un, Rec=irec) cwt2, cwat2, cwsurft2
         End Do
     ! synchronize for common number of w-points to all processes
         If (iw <= nwdf/procs) Call barrier
      End Do
      Close (un)
      Deallocate (cw, cwa, cwsurf)
      Allocate (cw(nwdfp, nst1, nst2))
      Allocate (cwa(nwdfp, nst1, nst2))
      Allocate (cwsurf(nwdfp, nst1, nst2))
      Allocate (cwsurft1(nwdfp), cwt1(nwdfp), cwat1(nwdfp))
  ! open temporary file for reading
      Open (un, File=trim(filnamt), Form='unformatted', Action='read', &
     & Status='old', Access='direct', Recl=Recl)
      irec = 0
      Do ik = 1, nkpt
         Do iw = 1, nwdfp
            irec = irec + 1
            Read (un, Rec=irec) cwt2, cwat2, cwsurft2
            cw (iw, :, :) = cwt2 (:, :)
            cwa (iw, :, :) = cwat2 (:, :)
            cwsurf (iw, :, :) = cwsurft2 (:, :)
         End Do
         Do ist1 = 1, nst1
            Do ist2 = 1, nst2
               cwsurft1 (:) = cwsurf (:, ist1, ist2)
               cwt1 (:) = cw (:, ist1, ist2)
               cwat1 (:) = cwa (:, ist1, ist2)
           ! routine cares for record position
               Call puttetcw (iq, ik, ist1, ist2, nst1, nst2, filnam, &
              & cwt1, cwat1, cwsurft1)
            End Do
         End Do
      End Do
      Close (un)
      Call filedel (trim(filnamt))
      Deallocate (cwt2, cwat2, cwsurft2)
      Deallocate (cw, cwa, cwsurf, eb)
      Deallocate (cwt1, cwat1, cwsurft1)
      Deallocate (w, wreal)
      Write (unitout, '(a)') 'Info(' // trim (thisnam) // '): weights f&
     &or tetrahedron method finished.'
End Subroutine tetcalccwq
