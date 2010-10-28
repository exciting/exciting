
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: exccoulint
! !INTERFACE:
Subroutine exccoulint
! !USES:
      Use modmain
      Use modinput
      Use modmpi
      Use modxs
      Use ioarray
      Use m_xsgauntgen
      Use m_findgntn0
      Use m_writegqpts
      Use m_genfilname
      Use m_getunit
! !DESCRIPTION:
!   Calculates the exchange term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!      level for the treatment of core excitations (using local orbitals).
!      October 2010 (Weine Olovsson)
!EOP
!BOC      
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'exccoulint'
      Character (256) :: fnexcli
      Integer, Parameter :: iqmt = 1
      Real (8), Parameter :: epsortho = 1.d-12
      Integer :: iknr, jknr, iqr, iq, igq1, n, un
      Integer :: iv (3), j1, j2
      Integer :: ist1, ist2, ist3, ist4, nst12, nst34, nst13, nst24, &
     & ikkp, nkkp
      Integer :: sta1, sto1, sta2, sto2, rnst1, rnst2, rnst3, rnst4
      Real (8), Allocatable :: potcl (:)
      Complex (8), Allocatable :: exclit (:, :), excli (:, :, :, :)
      Complex (8), Allocatable :: emat12 (:, :), emat34 (:, :)
      Complex (8), Allocatable :: emat12k (:, :, :, :)
  !---------------!
  !   main part   !
  !---------------!
      input%xs%emattype = 1
      Call init0
      Call init1
      Call init2
  ! set the range of valence/core and conduction states to use
      sta1 = input%xs%bse%nstlbsemat(1)
      sto1 = input%xs%bse%nstlbsemat(2)
      sta2 = input%xs%bse%nstlbsemat(3)
      sto2 = input%xs%bse%nstlbsemat(4)      
      rnst1 = sto1-sta1+1
      rnst2 = sto2-sta2+1
      rnst3 = sto2-sta2+1
      rnst4 = sto1-sta1+1
  ! read Fermi energy from file
      Call readfermi
  ! save variables for the Gamma q-point
      Call xssave0
  ! generate Gaunt coefficients
      Call xsgauntgen (Max(input%groundstate%lmaxapw, lolmax), &
     & input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
  ! find indices for non-zero Gaunt coefficients
      Call findgntn0 (Max(input%xs%lmaxapwwf, lolmax), &
     & Max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
      Write (unitout, '(a,3i8)') 'Info(' // thisnam // '): Gaunt coeffi&
     &cients generated within lmax values:', input%groundstate%lmaxapw, &
     & input%xs%lmaxemat, input%groundstate%lmaxapw
      Write (unitout, '(a, i6)') 'Info(' // thisnam // '): number of q-&
     &points: ', nqpt
      Call flushifc (unitout)
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      Call findocclims (0, istocc0, istocc, istunocc0, istunocc, isto0, &
     & isto, istu0, istu)
  ! only for systems with a gap in energy
      If ( .Not. ksgap) Then
         Write (*,*)
         Write (*, '("Warning(",a,"): There is no KS-gap present&
        &")') trim (thisnam)
         Write (*,*)
      End If
  ! check number of empty states
      If (input%xs%screening%nempty .Lt. input%groundstate%nempty) Then
         Write (*,*)
         Write (*, '("Error(",a,"): too few empty states in screening e&
        &igenvector file - the screening should include many empty stat&
        &es (BSE/screening)", 2i8)') trim (thisnam), &
        & input%groundstate%nempty, input%xs%screening%nempty
         Write (*,*)
         Call terminate
      End If
      Call ematbdcmbs (input%xs%emattype)
      nst12 = rnst1 * rnst2
      nst34 = rnst3 * rnst4
      nst13 = rnst1 * rnst3
      nst24 = rnst2 * rnst4
      Call genfilname (dotext='_SCI.OUT', setfilext=.True.)
      If (rank .Eq. 0) Then
         Call writekpts
         Call writeqpts
      End If
      n = ngq (iqmt)
      Call ematrad (iqmt)
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      Allocate (potcl(n))
      Allocate (excli(rnst1, rnst2, rnst1, rnst2))
      Allocate (exclit(nst12, nst34))
      Allocate (emat12k(nst1, nst2, n, nkptnr))
      Allocate (emat12(nst12, n), emat34(nst34, n))
      potcl (:) = 0.d0
      excli (:, :, :, :) = zzero
  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
      Call genparidxran ('k', nkptnr)
      Call init1offs (qvkloff(1, iqmt))
      Call ematqalloc
      Do iknr = kpari, kparf
         Call chkpt (3, (/ task, 1, iknr /), 'task,sub,k-point; matrix &
        &elements of plane wave')
     ! matrix elements for k and q=0
         Call ematqk1 (iqmt, iknr)
         emat12k (:, :, :, iknr) = xiou (:, :, :)
         Deallocate (xiou, xiuo)
      End Do
  ! communicate array-parts wrt. k-points
      call mpi_allgatherv_ifc(nkptnr,nst1*nst2*n,zbuf=emat12k)
      input%xs%emattype = 1
      Call ematbdcmbs (input%xs%emattype)
  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
      nkkp = (nkptnr*(nkptnr+1)) / 2
      Call genparidxran ('p', nkkp)
      Call genfilname (basename='EXCLI', asc=.True., filnam=fnexcli)
      Call getunit (un)
      If (rank .Eq. 0) open (un, file=trim(fnexcli), form='formatted', &
     & action='write', status='replace')
!
      Do ikkp = ppari, pparf
         Call chkpt (3, (/ task, 2, ikkp /), 'task,sub,(k,kp)-pair; exc&
        &hange term of BSE Hamiltonian')
         Call kkpmap (ikkp, nkptnr, iknr, jknr)
         iv (:) = ivknr (:, jknr) - ivknr (:, iknr)
         iv (:) = modulo (iv(:), input%groundstate%ngridk(:))
     ! q-point (reduced)
         iqr = iqmapr (iv(1), iv(2), iv(3))
     ! q-point (non-reduced)
         iq = iqmap (iv(1), iv(2), iv(3))
!
     ! set G=0 term of Coulomb potential to zero [Ambegaokar-Kohn]
         potcl (1) = 0.d0
     ! set up Coulomb potential
         Do igq1 = 2, n
            Call genwiqggp (0, iqmt, igq1, igq1, potcl(igq1))
         End Do
!
         Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
         j1 = 0
         Do ist2 = sta2, sto2
            Do ist1 = sta1, sto1
               j1 = j1 + 1
               emat12 (j1, :) = emat12k (ist1, ist2, :, iknr)
            End Do
         End Do
         j2 = 0
         Do ist4 = sta2, sto2
            Do ist3 = sta1, sto1
               j2 = j2 + 1
               emat34 (j2, :) = emat12k (ist3, ist4, :, jknr) * potcl &
              & (:)
            End Do
         End Do
!
     ! calculate exchange matrix elements: V_{1234} = M_{12}^* M_{34}^T
         emat12 = conjg (emat12)
         Call zgemm ('n', 't', nst12, nst12, n, zone/omega/nkptnr, &
        & emat12, nst12, emat34, nst12, zzero, exclit, nst12)
!
     ! map back to individual band indices
         j2 = 0
         Do ist4 = 1, rnst2
            Do ist3 = 1, rnst1
               j2 = j2 + 1
               j1 = 0
               Do ist2 = 1, rnst2
                  Do ist1 = 1, rnst1
                     j1 = j1 + 1
                     excli (ist1, ist2, ist3, ist4) = exclit (j1, j2)
                  End Do
               End Do
            End Do
         End Do
!
         If ((rank .Eq. 0) .And. (ikkp .Le. 3)) Then
            Do ist1 = 1, rnst1
               Do ist2 = 1, rnst2
                  Do ist3 = 1, rnst1
                     Do ist4 = 1, rnst2
                        Write (un, '(i5,3x,3i4,2x,3i4,2x,4e18.10)') &
                       & ikkp, iknr, ist1+sta1-1, ist2+sta2-1, jknr,&
                       & ist3+sta2-1, ist4+sta1-1, &
                       & excli (ist1, ist2, ist3, ist4), Abs &
                       & (excli(ist1, ist2, ist3, ist4))
                     End Do
                  End Do
               End Do
            End Do
         End If
!
     ! parallel write
         Call putbsemat ('EXCLI.OUT', excli, ikkp, iknr, jknr, iq, iqr, &
        & rnst1, rnst2, rnst4, rnst3)
         Call genfilname (dotext='_SCI.OUT', setfilext=.True.)
!
     ! end loop over (k,kp) pairs
      End Do
      If (rank .Eq. 0) write (un, '("# ikkp, iknr,ist1,ist3, jknr,ist2,&
     &ist4,    Re(V),            Im(V),             |V|^2")')
      If (rank .Eq. 0) close (un)
!
      Call barrier
      Call findgntn0_clear
      Deallocate (emat12k, exclit, emat12, emat34)
      Deallocate (potcl, excli)
!
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): Exchange &
     &Coulomb interaction finished"
End Subroutine exccoulint
!EOC

