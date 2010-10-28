
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: scrcoulint
! !INTERFACE:
Subroutine scrcoulint
! !USES:
      Use modmain
      Use modinput
      Use modmpi
      Use modxs
      Use summations
      Use m_xsgauntgen
      Use m_findgntn0
      Use m_writevars
      Use m_genfilname
      Use m_getunit
! !DESCRIPTION:
!   Calculates the direct term of the Bethe-Salpeter Hamiltonian.
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
      Character (*), Parameter :: thisnam = 'scrcoulint'
      Character (256) :: fnsccli, fnscreeninv
      Real (8), Parameter :: epsortho = 1.d-12
      Integer :: ikkp, iknr, jknr, iqr, iq, iqrnr, jsym, jsymi, igq1, &
     & igq2, n, recl, un
      Integer :: nsc, iv (3), ivgsym (3), j1, j2, nkkp
      Integer :: ist1, ist2, ist3, ist4, nst12, nst34, nst13, nst24
      Integer :: sta1, sto1, sta2, sto2, rnst1, rnst2, rnst3, rnst4
      Logical :: tq0, tphf
      Real (8) :: vqr (3), vq (3), t1
      Integer :: sc (maxsymcrys), ivgsc (3, maxsymcrys)
      Integer, Allocatable :: igqmap (:)
      Complex (8) :: zt1
      Complex (8), Allocatable :: scclit (:, :), sccli (:, :, :, :), &
     & scclid (:, :)
      Complex (8), Allocatable :: scieffg (:, :, :), tm (:, :), tmi (:, &
     & :), bsedt (:, :)
      Complex (8), Allocatable :: phf (:, :), emat12 (:, :), emat34 (:, &
     & :)
  ! external functions
      Integer, External :: idxkkp
      Logical, External :: tqgamma
  !---------------!
  !   main part   !
  !---------------!
!
      input%xs%emattype = 2
      Call init0
      Call init1
      Call init2
 ! set the range of valence/core and conduction states to use
      sta1 = input%xs%bse%nstlbsemat(1)
      sto1 = input%xs%bse%nstlbsemat(2)
      sta2 = input%xs%bse%nstlbsemat(3)
      sto2 = input%xs%bse%nstlbsemat(4)
      rnst1 = sto1-sta1+1
      rnst2 = sto1-sta1+1
      rnst3 = sto2-sta2+1
      rnst4 = sto2-sta2+1
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
!
  ! local arrays
      Allocate (phf(ngqmax, ngqmax))
      Allocate (sccli(rnst1, rnst3, rnst2, rnst4), scclid(rnst1, rnst3))
      Allocate (scieffg(ngqmax, ngqmax, nqptr))
      sccli (:, :, :, :) = zzero
      scieffg (:, :, :) = zzero
!
  ! set file extension
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
!
  !-----------------------------------!
  !     loop over reduced q-points    !
  !-----------------------------------!
      Call getunit (un)
      Call genparidxran ('q', nqptr)
!
      Do iqr = qpari, qparf
         Call genfilname (basename='W_SCREEN', iq=iqr, &
        & filnam=fnscreeninv)
         Open (un, File=trim(fnscreeninv), Form='formatted', Action='wr&
        &ite', Status='replace')
         Call chkpt (3, (/ task, 1, iqr /), 'task,sub,reduced q-point; &
        &generate effective screened Coulomb potential')
     ! locate reduced q-point in non-reduced set
         iqrnr = iqmap (ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))
         n = ngq (iqrnr)
!
     ! calculate effective screened Coulomb interaction
         Call genscclieff (iqr, ngqmax, n, scieffg(1, 1, iqr))
         Do igq1 = 1, n
            Do igq2 = 1, n
               Write (un, '(2i8,3g18.10)') igq1, igq2, scieffg (igq1, &
              & igq2, iqr), Abs (scieffg(igq1, igq2, iqr))
            End Do
         End Do
         Call writevars (un, iqr, 0)
         Close (un)
!
     ! generate radial integrals for matrix elements of plane wave
         Call putematrad (iqr, iqrnr)
      End Do
  ! communicate array-parts wrt. q-points
      call mpi_allgatherv_ifc(nqptr,ngqmax*ngqmax,zbuf=scieffg)
      Call barrier
!
  ! information on size of output file
      nkkp = (nkptnr*(nkptnr+1)) / 2
      Inquire (IoLength=Recl) ikkp, iknr, jknr, iq, iqr, nst1, nst2, &
     & nst3, nst4, sccli
      Write (unitout,*)
      Write (unitout, '(a,f12.3)') 'file size for screened Coulomb inte&
     &raction (GB):', recl * nkkp / 1024.d0 ** 3
      Write (unitout,*)
!
  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
      Call genfilname (basename='SCCLI', asc=.True., filnam=fnsccli)
      Call getunit (un)
      If (rank .Eq. 0) open (un, file=trim(fnsccli), form='formatted', &
     & action='write', status='replace')
      Call genparidxran ('p', nkkp)
      Allocate (bsedt(3, 0:procs-1))
      bsedt (1, :) = 1.d8
      bsedt (2, :) = - 1.d8
      bsedt (3, :) = zzero
  ! loop over combinations of k-points
      Do ikkp = ppari, pparf
         Call chkpt (3, (/ task, 2, ikkp /), 'task,sub,(k,kp)-pair; dir&
        &ect term of BSE Hamiltonian')
         Call kkpmap (ikkp, nkptnr, iknr, jknr)
     ! k-point difference
         iv (:) = ivknr (:, jknr) - ivknr (:, iknr)
         iv (:) = modulo (iv(:), ngridq(:))
     ! q-point (reduced)
         iqr = iqmapr (iv(1), iv(2), iv(3))
         vqr (:) = vqlr (:, iqr)
     ! q-point (non-reduced)
         iq = iqmap (iv(1), iv(2), iv(3))
         vq (:) = vql (:, iq)
     ! local field effects size
         tq0 = tqgamma (iq)
         n = ngq (iq)
!
         Allocate (igqmap(n), emat12(nst12, n), emat34(nst34, n))
         Allocate (tm(n, n), tmi(n, n))
         Allocate (scclit(nst12, nst34))
!
     ! find symmetry operations that reduce the q-point to the irreducible
     ! part of the Brillouin zone
         Call findsymeqiv (input%xs%BSE%fbzq, vq, vqr, nsc, sc, ivgsc)
!
     ! find the map that rotates the G-vectors
         Call findgqmap (iq, iqr, nsc, sc, ivgsc, n, jsym, jsymi, &
        & ivgsym, igqmap)
     ! generate phase factor for dielectric matrix due to non-primitive
     ! translations
         Call genphasedm (iq, jsym, ngqmax, n, phf, tphf)
!
     ! get radial integrals
         Call getematrad (iqr, iq)
     ! rotate radial integrals
         Call rotematrad (n, igqmap)
     ! rotate inverse of screening, Coulomb potential and radial integrals
         tmi (:, :) = phf (:n, :n) * scieffg (igqmap, igqmap, iqr)
!
     ! calculate matrix elements of the plane wave
         input%xs%emattype = 2
         Call ematbdcmbs (input%xs%emattype)
         Call ematqalloc
         Call ematqk1 (iq, iknr)
         input%xs%emattype = 2
         Call ematbdcmbs (input%xs%emattype)
         Call chkpt (3, (/ task, 2, ikkp /), 'task,sub,(k,kp)-pair; dir&
        &ect term of BSE Hamiltonian')
!
     ! select screening level
         tm (:, :) = zzero
         Select Case (trim(input%xs%screening%screentype))
         Case ('longrange')
        ! constant screening (q=0 average tensor)
            Forall (igq1=1:n)
               tm (igq1, igq1) = fourpi * dot_product (vgqc(:, igq1, &
              & iq), matmul(dielten, vgqc(:, igq1, iq)))
            End Forall
         Case ('diag')
        ! only diagonal of screening
            Forall (igq1=1:n)
               tm (igq1, igq1) = tmi (igq1, igq1)
            End Forall
         Case ('full')
        ! full screening
            tm (:, :) = tmi (:, :)
         End Select
!
     ! combine indices for matrix elements of plane wave
         j1 = 0
         Do ist2 = sta1, sto1
            Do ist1 = sta1, sto1
               j1 = j1 + 1
               emat12 (j1, :) = xiou (ist1, ist2, :)
            End Do
         End Do
         j2 = 0
         Do ist4 = sta2, sto2
            Do ist3 = sta2, sto2
               j2 = j2 + 1
               emat34 (j2, :) = xiuo (ist3, ist4, :)
            End Do
         End Do
!
     ! matrix elements of direct term (as in BSE-code of Peter and
     ! in the SELF-documentation of Andrea Marini)
         scclit = matmul (conjg(emat12), matmul(tm, transpose(emat34))) &
        & / omega / nkptnr
!
     ! map back to individual band indices
         j2 = 0
         Do ist4 = 1, rnst4
            Do ist3 = 1, rnst3
               j2 = j2 + 1
               j1 = 0
               Do ist2 = 1, rnst2
                  Do ist1 = 1, rnst1
                     j1 = j1 + 1
                     sccli (ist1, ist3, ist2, ist4) = scclit (j1, j2)
                  End Do
               End Do
            End Do
         End Do
         If ((rank .Eq. 0) .And. (ikkp .Le. 3)) Then
        ! write to ASCII file
            Do ist1 = 1, rnst1
               Do ist3 = 1, rnst3
                  Do ist2 = 1, rnst2
                     Do ist4 = 1, rnst4
                        zt1 = sccli (ist1, ist3, ist2, ist4)
                        Write (un, '(i5,3x,3i4,2x,3i4,2x,4e18.10)') &
                       & ikkp, iknr, ist1+sta1-1, ist3+sta2-1, jknr,& 
                       & ist2+sta1-1, ist4+sta2-1, zt1, &
                       & Abs (zt1) ** 2, Atan2 (aimag(zt1), dble(zt1)) &
                       & / pi
                     End Do
                  End Do
               End Do
            End Do
         End If
     ! analyze BSE diagonal
         If (iknr .Eq. jknr) Then
            Do ist1 = 1, rnst1
               Do ist3 = 1, rnst3
                  zt1 = sccli (ist1, ist3, ist1, ist3)
                  scclid (ist1, ist3) = zt1
                  t1 = dble (zt1)
                  bsedt (1, rank) = Min (dble(bsedt(1, rank)), t1)
                  bsedt (2, rank) = Max (dble(bsedt(2, rank)), t1)
                  bsedt (3, rank) = bsedt (3, rank) + zt1 / (rnst1*rnst3)
               End Do
            End Do
         End If
!
     ! parallel write
         Call putbsemat ('SCCLI.OUT', sccli, ikkp, iknr, jknr, iq, iqr, &
        & rnst1, rnst3, rnst2, rnst4)
!
         Deallocate (igqmap, emat12, emat34)
         Deallocate (tm, tmi)
         Deallocate (scclit)
!
     ! end loop over (k,kp)-pairs
      End Do
      If (rank .Eq. 0) write (un, '("# ikkp, iknr,ist1,ist3, jknr,ist2,&
     &ist4,    Re(W),            Im(W),             |W|^2,           an&
     &g/pi")')
      If (rank .Eq. 0) close (un)
      Call barrier
  ! communicate array-parts wrt. q-points
      call mpi_allgatherv_ifc(procs,3,zbuf=bsedt)
  ! BSE kernel diagonal parameters
      bsedl = minval (dble(bsedt(1, :)))
      bsedu = maxval (dble(bsedt(2, :)))
      bsedd = bsedu - bsedl
      bsed = sum (bsedt(3, :)) / nkptnr
      Deallocate (bsedt, scclid)
  ! write BSE kernel diagonal parameters
      If (rank .Eq. 0) Call putbsediag ('BSEDIAG.OUT')
!
      Call findgntn0_clear
!
      Write (unitout, '("Info(scrcoulint): Screened Coulomb interaction&
     & finished")')
!
End Subroutine scrcoulint
!EOC

