!
!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kernxc_bse
! !INTERFACE:
!
!
Subroutine kernxc_bse
! !USES:
      Use modinput
      Use modmain
#ifdef TETRA      
      Use modtetra
#endif
      Use modxs
      Use m_xsgauntgen
      Use m_findgntn0
      Use m_writegqpts
      Use m_genwgrid
      Use m_xszoutpr3
      Use m_getpemat
      Use m_getunit
      Use m_genfilname
! !INPUT/OUTPUT PARAMETERS:
!   oct   : optical diagonal tensor component (in,integer)
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'kernxc_bse'
      Integer, Parameter :: iqmt = 1, noptc = 3
      Character (256) :: filnam2, filnam3, filnam4
      Integer :: iw, wi, wf, nwdfp, n, recl, un, un2, un3, j1, j2, oct
      Integer :: ikkp, iknr, jknr, iknrq, jknrq, igq1, igq2
      Integer :: ist1, ist2, ist3, ist4, nst12, nst34, nst13, nst24
      Real (8) :: t1, brd
      Real (8) :: cpu_init1offs, cpu_ematrad, cpu_ematqalloc, &
     & cpu_ematqk1
      Real (8) :: cpu_ematqdealloc, cpu_clph, cpu_suma, cpu_write
      Complex (8) :: zt1
  ! allocatable arrays
      Real (8), Allocatable :: dek (:, :), dok (:, :), scisk (:, :)
      Real (8), Allocatable :: dekp (:, :), dokp (:, :), sciskp (:, :)
      Real (8), Allocatable :: deval (:, :, :), docc (:, :, :), scis &
     & (:, :, :)
      Real (8), Allocatable :: dde (:, :)
      Complex (8), Allocatable :: zmr (:, :), zmq (:, :), zmra (:, :), &
     & zmqa (:, :)
      Complex (8), Allocatable :: scclit (:, :), sccli (:, :, :, :), &
     & scclih (:, :, :, :)
      Complex (8), Allocatable :: emat (:, :, :, :), emata (:, :, :, :)
      Complex (8), Allocatable :: den1 (:), den2 (:), den1a (:), den2a &
     & (:)
      Complex (8), Allocatable :: emat12p (:, :), emat12pa (:, :)
      Complex (8), Allocatable :: emat12k (:, :, :), emat12kp (:, :, :)
      Complex (8), Allocatable :: emat12ka (:, :, :), emat12kpa (:, :, &
     & :)
      Complex (8), Allocatable :: residr (:, :), residq (:, :), osca &
     & (:, :), oscb (:, :)
      Complex (8), Allocatable :: residra (:, :), residqa (:, :), oscaa &
     & (:, :), oscba (:, :)
      Complex (8), Allocatable :: fxc (:, :, :), w (:), bsedg (:, :), &
     & bufou (:, :, :), bufuo (:, :, :), pufou (:, :, :), pufuo (:, :, &
     & :)
  ! external functions
      Integer, External :: idxkkp, l2int
      Logical, External :: tqgamma
  ! check that if Kohn-Sham response is time-ordered, so is the setting for the
  ! kernel
      if (input%xs%tddft%tordfxc .neqv. input%xs%tddft%torddf) then
         write(*,*)
         write(*,'("Error(kernxc_bse): Both, the Kohn-Sham response function")')
         write(*,'(" and the BSE-derived xc kernel have to be either causal or time-ordered.")')
         write(*,*)
         call terminate
      end if
      brd = input%xs%broad
      input%xs%emattype = 2
      Call init0
      Call init1
      Call init2
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
      Write (unitout, '(a, 3i8)') 'Info(' // thisnam // '): Gaunt coeff&
     &icients generated within lmax values:', &
     & input%groundstate%lmaxapw, input%xs%lmaxemat, &
     & input%groundstate%lmaxapw
      Write (unitout, '(a, i6)') 'Info(' // thisnam // '): number of q-&
     &points: ', nqpt
      Call flushifc (unitout)
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      Call findocclims (0, istocc0, istocc, istunocc0, istunocc, isto0, &
     & isto, istu0, istu)
  ! only for systems with a gap in energy
      If ( .Not. ksgap) Then
         Write (*,*)
         Write (*, '("Error(", a, "): screened Coulomb interaction work&
        &s only for systems with KS-gap.")') trim (thisnam)
         Write (*,*)
         Call terminate
      End If
  ! check number of empty states
      If (input%xs%screening%nempty .Lt. input%groundstate%nempty) Then
         Write (*,*)
         Write (*, '("Error(", a, "): too few empty states in screening&
        & eigenvector file - the screening should include many empty st&
        &ates (BSE/screening)", 2i8)') trim (thisnam), &
        & input%groundstate%nempty, input%xs%screening%nempty
         Write (*,*)
         Call terminate
      End If

      if (.true.) then
         Call ematbdcmbs (input%xs%emattype)
      else
         ! band limits below (as in bse.F90) should replace those generated by ematbdcmbs above.
         ! However, consistency with reading of pmat and other arrays below should be checked first.
         nst1 = input%xs%bse%nstlbse(2) - input%xs%bse%nstlbse(1) + 1
         nst2 = nst1
         nst3 = input%xs%bse%nstlbse(4) - input%xs%bse%nstlbse(3) + 1
         nst4 = nst3

         istl1 = input%xs%bse%nstlbse(1)
         istu1 = input%xs%bse%nstlbse(2)
         istl2 = input%xs%bse%nstlbse(1)
         istu2 = input%xs%bse%nstlbse(2)
         istl3 = input%xs%bse%nstlbse(3) + istocc0 - 1
         istu3 = input%xs%bse%nstlbse(4) + istocc0 - 1
         istl4 = istl3
         istu4 = istu3

      end if

      nst12 = nst1 * nst2
      nst34 = nst3 * nst4
      nst13 = nst1 * nst3
      nst24 = nst2 * nst4

!
      Call genparidxran ('w', nwdf)
  ! sampling type for Brillouin zone sampling
      bzsampl = l2int (input%xs%tetra%tetradf)
  ! limits for w-points
      wi = wpari
      wf = wparf
      nwdfp = wparf - wpari + 1
  ! matrix size for local field effects (first q-point is Gamma-point)
      n = ngq (iqmt)
!
      Allocate (bufou(nst1, nst3, n))
      Allocate (bufuo(nst3, nst1, n))
      Allocate (pufou(3, nst1, nst3))
      Allocate (pufuo(3, nst3, nst1))
  ! allocate global arrays
      If (allocated(xiou)) deallocate (xiou)
      Allocate (xiou(nst1, nst3, n))
      If (allocated(xiuo)) deallocate (xiuo)
      Allocate (xiuo(nst3, nst1, n))
      If (allocated(pmou)) deallocate (pmou)
      Allocate (pmou(3, nst1, nst3))
      If (allocated(pmuo)) deallocate (pmuo)
      Allocate (pmuo(3, nst3, nst1))
      If (allocated(deou)) deallocate (deou)
      Allocate (deou(nst1, nst3))
      If (allocated(deuo)) deallocate (deuo)
      Allocate (deuo(nst3, nst1))
      If (allocated(docc12)) deallocate (docc12)
      Allocate (docc12(nst1, nst3))
      If (allocated(docc21)) deallocate (docc21)
      Allocate (docc21(nst3, nst1))
  ! allocate local arrays
      Allocate (emat12p(nst13,-3:n), zmr(nst13, nst13), zmq(nst13, &
     & nst13))
      Allocate (emat12pa(nst13,-3:n), zmra(nst13, nst13), zmqa(nst13, &
     & nst13))
      Allocate (dek(nst1, nst3), dekp(nst1, nst3), dde(nst1, nst3))
      Allocate (dok(nst1, nst3), dokp(nst1, nst3))
      Allocate (scisk(nst1, nst3), sciskp(nst1, nst3))
      Allocate (fxc(-3:n,-3:n, nwdf))
      Allocate (sccli(nst1, nst3, nst2, nst4))
      Allocate (scclih(nst1, nst3, nst2, nst4))
      Allocate (scclit(nst13, nst13))
      Allocate (emat12k(-3:n, nst1, nst3), emat12kp(nst1, nst3,-3:n))
      Allocate (residr(nst13,-3:n), residq(nst13,-3:n))
      Allocate (emat12ka(-3:n, nst3, nst1), emat12kpa(nst3, nst1,-3:n))
      Allocate (residra(nst13,-3:n), residqa(nst13,-3:n))
      Allocate (w(nwdf))
      Allocate (osca(-3:n,-3:n), oscb(-3:n,-3:n))
      Allocate (oscaa(-3:n,-3:n), oscba(-3:n,-3:n))
      Allocate (den1(nwdf), den2(nwdf), den1a(nwdf), den2a(nwdf))
      fxc (:, :, :) = zzero
      sccli (:, :, :, :) = zzero
      Allocate (emat(nst1, nst3, n, nkptnr))
      Allocate (emata(nst3, nst1, n, nkptnr))
      Allocate (deval(nst1, nst3, nkptnr))
      Allocate (docc(nst1, nst3, nkptnr))
      Allocate (scis(nst1, nst3, nkptnr))
      Allocate (bsedg(nst1, nst3))
!
!
      If ((input%xs%tddft%fxctypenumber .Eq. 7) .Or. &
     & (input%xs%tddft%fxctypenumber .Eq. 8)) Then
         Call getbsediag
         Write (unitout, '("Info(", a, "): read diagonal of BSE kernel"&
        &)') trim (thisnam)
         Write (unitout, '(" mean value : ", 2g18.10)') bsed
      End If
!
  ! generate energy grid
      Call genwgrid (nwdf, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_cmplx=w)
!
  ! precalculate matrix elements
      input%xs%emattype = 1
      Call ematbdcmbs (input%xs%emattype)
      Call ematrad (iqmt)
      Call ematqalloc
!
  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
      Do iknr = 1, nkptnr
         Call chkpt (3, (/ task, 1, iknr /), 'task, sub, k - point; gen&
        &erate matrix elements of plane wave')
         iknrq = ikmapikq (iknr, iqmt)
     ! matrix elements for k and q=0
         Call ematqk1 (iqmt, iknr)
         emat (:, :, :, iknr) = xiou (:, :, :)
         emata (:, :, :, iknr) = xiuo (:, :, :)
         Deallocate (xiou, xiuo)
         Call getdevaldoccsv (iqmt, iknr, iknrq, istl1, istu1, istl2, &
        & istu2, deou, docc12, scisk)
         Call getdevaldoccsv (iqmt, iknr, iknrq, istl2, istu2, istl1, &
        & istu1, deuo, docc21, sciskp)
         deval (:, :, iknr) = deou (:, :)
         docc (:, :, iknr) = docc12 (:, :)
         scis (:, :, iknr) = scisk (:, :)
      End Do
      input%xs%emattype = 2
      Call ematbdcmbs (input%xs%emattype)
!
  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
      If (allocated(xiou)) deallocate (xiou)
      If (allocated(xiuo)) deallocate (xiuo)
      ikkp = 0
  ! first k-point
      Do iknr = 1, nkptnr
         Call chkpt (3, (/ task, 3, iknr /), 'task, sub, k-point; BSE-f&
        &xc-kernel')
         iknrq = ikmapikq (iknr, iqmt)
         bsedg (:, :) = bsed
         input%xs%emattype = 1
         Call ematbdcmbs (input%xs%emattype)
         Allocate (xiou(nst1, nst3, n))
         Allocate (xiuo(nst3, nst1, n))
         xiou (:, :, :) = emat (:, :, :, iknr)
         xiuo (:, :, :) = emata (:, :, :, iknr)
         deou (:, :) = deval (:, :, iknr)
         docc12 (:, :) = docc (:, :, iknr)
     ! apply gauge wrt. symmetrized Coulomb potential
         Call getpemat (iqmt, iknr, 'PMAT_SCR.OUT', '', m12=bufou, &
        & p12=pufou, m34=bufuo, p34=pufuo)
         dek (:, :) = deou (:, :)
         dok (:, :) = docc12 (:, :)
     ! add BSE diagonal
         scisk (:, :) = scis (:, :, iknr) + bsedg (:, :)
!
     ! assign optical components
         Do oct = 1, noptc
            emat12k (-oct, :, :) = pufou (oct, :, :)
            emat12ka (-oct, :, :) = pufuo (oct, :, :)
         End Do
         Do igq1 = 1, n
            emat12k (igq1, :, :) = bufou (:, :, igq1)
            emat12ka (igq1, :, :) = bufuo (:, :, igq1)
         End Do
!
         Deallocate (xiou, xiuo)
         input%xs%emattype = 2
         Call ematbdcmbs (input%xs%emattype)
!
         residr (:, :) = zzero
         residq (:, :) = zzero
         residra (:, :) = zzero
         residqa (:, :) = zzero
     ! second k-point
         Do jknr = 1, nkptnr
            jknrq = ikmapikq (jknr, iqmt)
!
            cpu_init1offs = 0.d0
            cpu_ematrad = 0.d0
            cpu_ematqalloc = 0.d0
            cpu_ematqk1 = 0.d0
            cpu_ematqdealloc = 0.d0
            cpu_clph = 0.d0
            cpu_suma = 0.d0
            cpu_write = 0.d0
!
            If (iknr .Le. jknr) Then
           ! index for upper triangle
               ikkp = idxkkp (iknr, jknr, nkptnr)
            Else
           ! swapped index for lower triangle
               ikkp = idxkkp (jknr, iknr, nkptnr)
            End If
!
            input%xs%emattype = 1
            Call ematbdcmbs (input%xs%emattype)
            Allocate (xiou(nst1, nst3, n))
            Allocate (xiuo(nst3, nst1, n))
            xiou (:, :, :) = emat (:, :, :, jknr)
            xiuo (:, :, :) = emata (:, :, :, jknr)
            deou (:, :) = deval (:, :, jknr)
            docc12 (:, :) = docc (:, :, jknr)
!
        ! apply gauge wrt. symmetrized Coulomb potential
            Call getpemat (iqmt, jknr, 'PMAT_SCR.OUT', '', m12=bufou, &
           & p12=pufou, m34=bufuo, p34=pufuo)
            dekp (:, :) = deou (:, :)
            dokp (:, :) = docc12 (:, :)
            sciskp (:, :) = scis (:, :, jknr)
        ! assign optical component
            Do oct = 1, noptc
               emat12kp (:, :,-oct) = pufou (oct, :, :)
               emat12kpa (:, :,-oct) = pufuo (oct, :, :)
            End Do
            emat12kp (:, :, 1:) = bufou (:, :, :)
            emat12kpa (:, :, 1:) = bufuo (:, :, :)
!
            Deallocate (xiou, xiuo)
            input%xs%emattype = 2
            Call ematbdcmbs (input%xs%emattype)
!
        ! get screened Coulomb interaction
            If (iknr .Le. jknr) Then
               Call getbsemat ('SCCLI.OUT', ikkp, nst1, nst3, sccli)
            Else
               Call getbsemat ('SCCLI.OUT', ikkp, nst1, nst3, scclih)
           ! use Hermitian property for lower triangle
               Do ist1 = 1, nst1
                  Do ist3 = 1, nst3
                     Do ist2 = 1, nst1
                        Do ist4 = 1, nst3
                           sccli (ist1, ist3, ist2, ist4) = conjg &
                          & (scclih(ist2, ist4, ist1, ist3))
                        End Do
                     End Do
                  End Do
               End Do
            End If
        ! proper sign of screened Coulomb interaction
            sccli = - sccli
        ! set diagonal of Bethe-Salpeter kernel to zero
        ! (cf. A. Marini, PRL 2003)
            If (iknr .Eq. jknr) Then
               Do ist3 = 1, nst3
                  Do ist1 = 1, nst1
                     sccli (ist1, ist3, ist1, ist3) = zzero
                  End Do
               End Do
            End If
            j1 = 0
            Do ist2 = 1, nst3
               Do ist1 = 1, nst1
                  j1 = j1 + 1
                  emat12p (j1, :) = conjg (emat12kp(ist1, ist2, :))
                  emat12pa (j1, :) = conjg (emat12kpa(ist2, ist1, :))
               End Do
            End Do
        ! map
            j2 = 0
            Do ist3 = 1, nst3
               Do ist1 = 1, nst1
                  j2 = j2 + 1
                  j1 = 0
                  Do ist4 = 1, nst3
                     Do ist2 = 1, nst1
                        j1 = j1 + 1
                        zt1 = sccli (ist1, ist3, ist2, ist4)
                    ! four point energy difference
                        t1 = dekp (ist2, ist4) - dek (ist1, ist3)
                    ! arrays for R- and Q-residuals
                        If (Abs(t1) .Ge. input%xs%tddft%fxcbsesplit) &
                       & Then
                           zmr (j2, j1) = zt1 / t1
                           zmq (j2, j1) = zzero
                           zmra (j2, j1) = conjg (zt1) / t1
                           zmqa (j2, j1) = zzero
                        Else
                           zmr (j2, j1) = zzero
                           zmq (j2, j1) = zt1
                           zmra (j2, j1) = zzero
                           zmqa (j2, j1) = conjg (zt1)
                        End If
!
                     End Do
                  End Do
               End Do
            End Do
        ! calculate residual "R"; partial fraction decomposition without
	! double poles
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))
            residr = residr + matmul (zmr, emat12p)
            residra = residra + matmul (zmra, emat12pa)
        ! calculate residual "Q"; double poles part
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))
            residq = residq + matmul (zmq, emat12p)
            residqa = residqa + matmul (zmqa, emat12pa)
!
        ! end inner loop over k-points
         End Do
!
     !--------------------------!
     !     set up BSE-kernel    !
     !--------------------------!
         t1 = 1.d0 / (nkptnr*omega)
         Do ist3 = 1, nst3
            Do ist1 = 1, nst1
               osca (:, :) = zzero
               oscb (:, :) = zzero
               oscaa (:, :) = zzero
               oscba (:, :) = zzero
               j1 = ist1 + (ist3-1) * nst1
           ! set up inner part of kernel
           ! generate oscillators
               Call xszoutpr3 (n+noptc+1, n+noptc+1, zone, emat12k(:, &
              & ist1, ist3), residr(j1, :), osca)
               Call xszoutpr3 (n+noptc+1, n+noptc+1, zone, emat12ka(:, &
              & ist3, ist1), residra(j1, :), oscaa)
!
           ! add Hermitian transpose
               Forall (igq1=-3:n, igq2=-3:n)
                  osca (igq1, igq2) = osca (igq1, igq2) + conjg &
                 & (osca(igq2, igq1))
                  oscaa (igq1, igq2) = oscaa (igq1, igq2) + conjg &
                 & (oscaa(igq2, igq1))
               End Forall
!
               Call xszoutpr3 (n+noptc+1, n+noptc+1, zone, emat12k(:, &
              & ist1, ist3), residq(j1, :), oscb)
               Call xszoutpr3 (n+noptc+1, n+noptc+1, zone, emat12ka(:, &
              & ist3, ist1), residqa(j1, :), oscba)
          ! set up energy denominators
               den1 (:) = 2.d0 * t1 / (w(:)+scisk(ist1, ist3)+dek(ist1, &
              & ist3)+zi*brd)
               den2 (:) = 2.d0 * t1 / (w(:)+scisk(ist1, ist3)+dek(ist1, &
              & ist3)+zi*brd) ** 2
               den1a (:) = 2.d0 * t1 / (w(:)+scisk(ist1, &
              & ist3)-dek(ist1, ist3)+torfxc*zi*brd)
               den2a (:) = - 2.d0 * t1 / (w(:)+scisk(ist1, &
              & ist3)-dek(ist1, ist3)+torfxc*zi*brd) ** 2
           ! update kernel
               Do iw = 1, nwdf
              ! resonant contributions
                  fxc (:, :, iw) = fxc (:, :, iw) + osca (:, :) * den1 &
                 & (iw) + oscb (:, :) * den2 (iw)
              ! antiresonant contributions
                  If (input%xs%tddft%aresfxc) fxc (:, :, iw) = fxc (:, &
                 & :, iw) + oscaa (:, :) * den1a (iw) + oscba (:, :) * &
                 & den2a (iw)
               End Do
           ! end loop over states #1
            End Do
        ! end loop over states #3
         End Do
     ! end outer loop over k-points
      End Do
!
  ! filename for xc-kernel (ASCII)
      Call genfilname (basename='KERNXC_BSE', asc=.True., bzsampl=bzsampl, &
     & acont=input%xs%tddft%acont, nar= .Not. input%xs%tddft%aresfxc, &
     & tord=input%xs%tddft%tordfxc, iqmt=iqmt, filnam=filnam2)
      Call getunit (un)
      Open (un, File=trim(filnam2), Form='formatted', Action='write', &
     & Status='replace')
!
  ! filename for xc-kernel
      Call genfilname (basename='FXC_BSE', asc=.False., &
     & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
     & input%xs%tddft%aresfxc, tord=input%xs%tddft%tordfxc, iqmt=iqmt, &
     & filnam=filnam3)
      Inquire (IoLength=Recl) n, fxc (-3:-1,-3:-1, 1), fxc (-3:-1, 1:, &
     & 1), fxc (1:,-3:-1, 1), fxc (1:, 1:, 1)
      Call getunit (un2)
      Open (un2, File=trim(filnam3), Form='unformatted', Action='write',&
     &  Status='replace', Access='direct', Recl=Recl)
!
  ! filename for xc-kernel
      Call genfilname (basename='FXC_BSE_HEAD', asc=.False., &
     & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
     & input%xs%tddft%aresfxc, tord=input%xs%tddft%tordfxc, iqmt=iqmt, &
     & filnam=filnam4)
      Call getunit (un3)
      Open (un3, File=trim(filnam4), Form='formatted', Action='write', &
     & Status='replace')
!
      Do iw = 1, nwdf
         Write (un2, Rec=iw) n, fxc (-3:-1,-3:-1, iw), fxc (-3:-1, 1:, &
        & iw), fxc (1:,-3:-1, iw), fxc (1:, 1:, iw)
         Write (un3, '(i6, 2x, g18.10, 2x, 6g18.10)') iw, dble (w(iw)), &
        & (fxc(-oct,-oct, iw), oct=1, noptc)
      End Do
      Do iw = 1, nwdf, 10
         Do igq1 = - noptc, n
            Do igq2 = - noptc, n
               Write (un, '(3i6, 3g18.10)') iw, igq1, igq2, fxc (igq1, &
              & igq2, iw), Abs (fxc(igq1, igq2, iw))
            End Do
         End Do
      End Do
      Close (un)
      Close (un2)
      Close (un3)
!
  ! deallocate
      Deallocate (den1, den2, den1a, den2a)
      Deallocate (emat12p, zmr, zmq, dek, dekp, dde, dok, dokp, scisk, &
     & fxc)
      Deallocate (sccli, scclih, scclit, emat12k, emat12kp, residr, &
     & residq, w, osca, oscb)
      Deallocate (emat, deval, docc, scis)
  ! deallocate antiresonant parts
      Deallocate (emata, emat12pa, emat12ka, emat12kpa, residra, &
     & residqa, zmra, zmqa)
      Deallocate (oscaa, oscba)
!
      Deallocate (bsedg)
      Deallocate (bufou, bufuo, pufou, pufuo)
!
End Subroutine kernxc_bse
!EOC
