!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematq (iq)
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use m_writegqpts
      Use m_filedel
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Character (*), Parameter :: thisnam = 'ematq'
      Integer :: ik
  ! filenames
      Call genfilname (basename='EMAT', iqmt=iq, &
     & etype=input%xs%emattype, filnam=fnemat)
      Call genfilname (basename='EMAT', iqmt=iq, &
     & etype=input%xs%emattype, procs=procs, rank=rank, &
     & filnam=fnemat_t)
      Call genfilname (nodotpar=.True., basename='EMAT_TIMING', &
     & iqmt=iq, etype=input%xs%emattype, procs=procs, rank=rank, &
     & filnam=fnetim)
  ! file extension for q-point
      Call genfilname (iqmt=iq, setfilext=.True.)
  ! calculate k+q and G+k+q related variables
      Call init1offs (qvkloff(1, iq))
  ! write G+q-vectors
      If (rank .Eq. 0) Then
         Call writegqpts (iq, filext)
         Call writekmapkq (iq)
      End If
  ! find highest (partially) occupied and lowest (partially) unoccupied states
      call findocclims(iq, ikmapikq(:,iq), istocc0, istunocc0, isto0, isto, istu0, istu)
      istunocc = istunocc0
      istocc = istocc0
      Call ematbdlims (1, nst1, istl1, istu1, nst2, istl2, istu2)
  ! generate radial integrals wrt. sph. Bessel functions
      Call ematrad (iq)
  ! delete timing information of previous runs
      Call filedel (trim(fnetim))
  ! write information
      Write (unitout, '(a, i6)') 'Info(' // thisnam // '): number of G+&
     &q vectors:', ngq (iq)
      Call ematqalloc
  ! loop over k-points
      Do ik = kpari, kparf
         Call ematqk1 (iq, ik)
      End Do

      Call ematqdealloc

      if (.not. input%sharedfs)call cpFileTonodes(fnemat)
      Call barrier
End Subroutine ematq
