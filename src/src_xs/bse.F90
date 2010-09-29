!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bse
! !INTERFACE:
!
!
Subroutine bse
! !USES:
      Use modinput
      Use modmain
      Use modxs
      Use m_genwgrid
      Use m_getpmat
      Use m_genfilname
      Use m_getunit
      Use m_writeeps
! !DESCRIPTION:
!   Solves the Bethe-Salpeter equation (BSE). The BSE is treated as equivalent
!   effective eigenvalue problem (thanks to the spectral theorem that can
!   be applied to the original BSE in the case of a statically screened Coulomb
!   interaction). The effective BSE-Hamiltonian consists of three parts
!   originating from different sources. It reads
!   $$ H^{\rm eff} = H^{\rm diag} + 2H^{\rm x} + H^{\rm c}, $$
!   where $H^{\rm diag}$ is the diagonal part stemming from the independent
!   particle transitions, $H^{\rm x}$ denotes the exchange-term caused by the
!   unscreened (bare) Coulomb interaction, whereas $H^{\rm c}$ accounts for the
!   particle-hole correlations and is originating from the screened Coulomb
!   interaction.
!   For the purpose of describing independent particle transitions with the
!   BSE only the diagonal term is referred to:
!   $$ H^{\rm eff} = H^{\rm diag}. $$
!   By neglecting the correlation part in the effective Hamiltonian we arrive
!   at the {\it random phase approximation} (RPA)
!   $$ H^{\rm eff} = H^{\rm diag} + 2H^{\rm x}. $$
!   Investigations on the spin-structure of the BSE-Hamiltonian show that there
!   are tow channels, namely the {\it singlet}-channel as solution to the
!   Hamiltonian
!   $$  H^{\rm eff} = H^{\rm diag} + 2H^{\rm x} + H^{\rm c} $$
!   and a {\it triplet} channel with the exchange-part being absent.
!   $$ H^{\rm eff} = H^{\rm diag} + H^{\rm c}. $$
!   The equation of the eigenvalue problem is given by
!   $$ \sum_{v'c'{\bf k'}} H^{\rm eff}_{vc{\bf k},v'c'{\bf k'}}
!       A^{\lambda}_{v'c'{\bf k'}}
!       =  \varepsilon_{\lambda} A^{\lambda}_{vc{\bf k}}. $$
!   For the diagonalization of the Hamiltonian, a LAPACK-routine ({\tt zheevx})
!   is invoked to obtain the eigenvalues $\varepsilon_{\lambda}$ and
!   eigenvectors $A^{\lambda}_{vc{\bf k}}$ (alternatively, a time-evolution
!   method shall be implemented to obtain the macroscopic dielectric function
!   directly).
!   Consequently, the transition amplitudes $t_{\lambda}$ are calculated
!   according to
!   $$ t^{i}_{\lambda} = \left|\sum_{vc{\bf k}} A^{\lambda}_{vc{\bf k}}
!      \frac{ p^{i}_{vc{\bf k}} }{ \varepsilon_{c{\bf k}}-
!                                  \varepsilon_{v{\bf k}} } \right|^2. $$
!   Here, the index $i$ labels the polarization and the matrix elements
!   $p^{i}_{vc{\bf k}}$ are the ones for the $i$-th component of the momentum
!   operator in Cartesian coordinates.
!   The macroscopic dielectric function (MDF) is obtained by the realation
!   $$ {\rm Im}\; \epsilon^{i}_{\rm M}(\omega) = \frac{8\pi^2}{V}
!                     \sum_{\lambda} t^{i}_{\lambda}
!                     \delta(\omega-\varepsilon_{\lambda}+\Delta),$$
!   where $\epsilon^{i}_{\rm M}$ is the MDF for the $i$-th polarization, $V$
!   denotes the crystal volume and $\Delta$ is a constant shift of the
!   conduction bands (scissors shift). The delta-function in the latter
!   expression is convoluted with a (symmetrized) Lorentzian
!   $$ \pi\delta(\omega-\omega_0) = \lim_{\eta\rightarrow 0} \left[
!                         \frac{\eta}{(\omega-\omega_0)^2+\eta^2} +
!                         \frac{\eta}{(-\omega-\omega_0)^2-\eta^2} \right] =
!     \pi\delta(\omega-\omega_0) +  \pi\delta(\omega+\omega_0)       $$
!   which is true for $\omega \ge 0$ if $\omega_0>0$. In doing so, the analytic
!   property ${\rm Im}\epsilon_{\rm M}(0)=0$ is fulfilled.
!   The broadening $\eta$ in the latter expression is adjusted by the
!   {\tt broad} parameter. (All parts of the documentation written by
!   S. Sagmeister are part of the author's PhD-thesis.)
!
! !REVISION HISTORY:
!   Created June 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Integer, Parameter :: iqmt = 0
      Integer, Parameter :: noptcmp = 3
      Real (8), Parameter :: epsortho = 1.d-12
      Character (256) :: fnexc, fnexcs
      Integer :: iknr, jknr, iqr, iq, iw, iv2 (3), s1, s2, hamsiz, &
     & nexc, ne
      Integer :: unexc, ist1, ist2, ist3, ist4, ikkp, oct, iv, ic, &
     & nvdif, ncdif
      Integer :: rnst1, rnst2, rnst3, rnst4 !(wol)
      Integer :: nrnst1, nrnst2, nrnst3, nrnst4 !(wol)
      Integer :: sta1, sto1, sta2, sto2 !(wol)
      Integer :: nsta1, nsto1, nsta2, nsto2 !(wol)
      Integer :: ist, jst !(wol)
      Real (8) :: de, egap, ts0, ts1
  ! allocatable arrays
      Integer, Allocatable :: sor (:)
      Real (8), Allocatable :: beval (:), w (:), oszsa (:)
      Real (8), Allocatable :: docc(:,:), kdocc (:) !(wol)
      Complex (8), Allocatable :: excli (:, :, :, :), sccli (:, :, :, &
     & :), ham (:, :)
      Complex (8), Allocatable :: bevec (:, :), pm (:, :, :), pmat (:), &
     & oszs (:), spectr (:)
  ! external functions
      Integer, External :: l2int
  ! *** TODO: symmetrize head of DM for spectrum
  !---------------------------!
  !     exciton variables     !   USE this ****************************
  !---------------------------!
  !!if (allocated(excite)) deallocate(excite)
  !!allocate(excite(nexcitmax,3))
  !!excite(:,:)=0.d0
  !!if (allocated(excito)) deallocate(excito)
  !!allocate(excito(nexcitmax,3))
  !!excito(:,:)=0.d0
      Call init0
      Call init1
      Call init2
      Call xssave0
  ! read Fermi energy from file
      Call readfermi
  ! initialize states below and above the Fermi energy
! initialize ranges for valence and conduction states !(wol)
      sta1 = input%xs%bse%nstlbse(1)
      sto1 = input%xs%bse%nstlbse(2)
      sta2 = input%xs%bse%nstlbse(3)
      sto2 = input%xs%bse%nstlbse(4)
      nsta1 = input%xs%bse%nstlbse2(1)
      nsto1 = input%xs%bse%nstlbse2(2)
      nsta2 = input%xs%bse%nstlbse2(3)
      nsto2 = input%xs%bse%nstlbse2(4)
      rnst1 = sto1-sta1+1
      rnst2 = sto1-sta1+1
      rnst3 = sto2-sta2+1
      rnst4 = sto2-sta2+1
      nrnst1 = nsto1-nsta1+1
      nrnst2 = nsto1-nsta1+1
      nrnst3 = nsto2-nsta2+1
      nrnst4 = nsto2-nsta2+1 
!(wol)--------end---------------------- 
!(wol)      nbfbse = input%xs%bse%nstlbse(1)
!(wol)      nafbse = input%xs%bse%nstlbse(2)
!(wol)- change this routine!      Call initocc (nbfbse, nafbse)
  ! use eigenvector files from screening-calculation
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      Call findocclims (iqmt, istocc0, istocc, istunocc0, istunocc, &
     & isto0, isto, istu0, istu)
!(wol)      nvdif = nstocc0 - nbfbse
!(wol)      ncdif = nstunocc0 - nafbse
!
!(wol)      Write (*,*) 'nbfbse, nafbse', nbfbse, nafbse
!(wol)      Write (*,*) 'nvdif, ncdif', nvdif, ncdif
!(wol)      Write (*,*) 'sta1, sto1, sta2, sto2', sta1, sto1, sta2, sto2 !(wol)
!(wol)      Write (*,*) 'nsta1, nsto1, nsta2, nsto2', nsta1, nsto1, nsta2, nsto2 !(wol)
!(wol)      Write (*,*) 'nkptnr', nkptnr
!(wol)      Write (*,*) 'nst1, nst2, nst3, nst4', nst1, nst2, nst3, nst4 !(wol)
!(wol)      Write (*,*) 'istocc0, istocc', istocc0, istocc !(wol)
!(wol)      Write (*,*) 'istunocc0, istunocc', istunocc0, istunocc !(wol)
!
  ! ****************************************************
      input%xs%emattype = 2
      Call ematbdcmbs (input%xs%emattype)
      Write (unitout,*)
      Write (unitout, '("Info(bse): information on number of states:")')
      Write (unitout, '(" number of states below Fermi energy in Hamilt&
     &onian:", i6)') nrnst1
!(wol)     &onian:", i6)') nbfbse
      Write (unitout, '(" number of states above Fermi energy in Hamilt&
     &onian:", i6)') nrnst3
!(wol)     &onian:", i6)') nafbse
      Write (unitout, '(" ranges of states according to BSE matrix:")')
      Write (unitout, '("  range of first index and number  :", 2i6, 3x&
     &, i6)') istl1, istu1, nst1
      Write (unitout, '("  range of second index and number :", 2i6, 3x&
     &, i6)') istl2, istu2, nst2
      Write (unitout, '("  range of third index and number  :", 2i6, 3x&
     &, i6)') istl3, istu3, nst3
      Write (unitout, '("  range of fourth index and number :", 2i6, 3x&
     &, i6)') istl4, istu4, nst4
!(wol)      If ((nvdif .Lt. 0) .Or. (ncdif .Lt. 0)) Then
!(wol)         Write (unitout,*)
!(wol)         Write (unitout, '("Error(bse): inconsistency in ranges of stat&
!(wol)        &es - check  routine for nvdif, ncdif")')
!(wol)         Write (unitout,*)
!(wol)         Call terminate
!(wol)      End If
      Write (*,*) 'nst1, nst2, nst3, nst4', nst1, nst2, nst3, nst4 !(wol)
  ! size of BSE-Hamiltonian
!(wol)      hamsiz = nbfbse * nafbse * nkptnr
      hamsiz = nrnst1 * nrnst3 * nkptnr
      Write (*,*) 'hamsiz', hamsiz !(wol)
  ! allocate arrays for Coulomb interactons
!(wol)      Allocate (sccli(nst1, nst3, nst2, nst4))
!(wol)      Allocate (excli(nst1, nst3, nst2, nst4))
      Allocate (sccli(rnst1, rnst3, rnst2, rnst4))
      Allocate (excli(rnst1, rnst3, rnst2, rnst4))
  ! allocate array for occupation number difference if needed !(wol) 
      Write (*,*) 'istl3, nstocc0, nsta2', istl3, nstocc0, nsta2 !(wol)
      	Allocate (docc (rnst1, rnst3))
      	Allocate (kdocc (hamsiz))
  ! allocate BSE-Hamiltonian (large matrix, up to several GB)
      Allocate (ham(hamsiz, hamsiz))
      ham (:, :) = zzero
  ! read in energies
      Do iknr = 1, nkptnr
         Call getevalsv (vkl(1, iknr), evalsv(1, iknr))
      End Do
  ! read mean value of diagonal of direct term
      bsed = 0.d0
      If ((trim(input%xs%bse%bsetype) .Eq. 'singlet') .Or. &
     & (trim(input%xs%bse%bsetype) .Eq. 'triplet')) Then
         If (input%xs%bse%bsedirsing) Then
            Call getbsediag
            Write (unitout, '("Info(bse): read diagonal of BSE kernel")&
           &')
            Write (unitout, '(" mean value : ", 2g18.10)') bsed
         End If
      End If
  ! determine gap
      egap = 1.d8
      Do iknr = 1, nkptnr
!(wol)         Do ist1 = 1 + nvdif, nst1
         Do ist1 = nsta1, nsto1
!(wol)            Do ist3 = 1, nst3 - ncdif
            Do ist3 = nsta2, nsto2
               egap = Min (egap, evalsv(ist3+istocc, iknr)-evalsv(ist1, &
              & iknr)+input%xs%scissor)
            End Do
         End Do
      End Do
      Write (unitout, '("Info(bse): gap:", g18.10)') egap
      If (egap .Lt. input%groundstate%epspot) Then
         Write (unitout,*)
!(wol)         Write (unitout, '("Error(bse): BSE needs system with gap")')
         Write (unitout, '("Warning(bse): system has no gap")')
         Write (unitout,*)
!(wol)         Call terminate
      End If
!TEST!!!!!!!
!	  If (istl3 - nsta2 + 1 .Le. nstocc0) Then
!	  egap = 0.0 !(wol)
!	  bsed = 0.0 !(wol)
!	  End If
!TEST!!!!!!!	  
  ! set up BSE-Hamiltonian
      ikkp = 0
      Do iknr = 1, nkptnr
         Do jknr = iknr, nkptnr
            ikkp = ikkp + 1
            iv2 (:) = ivknr (:, jknr) - ivknr (:, iknr)
            iv2 (:) = modulo (iv2(:), input%groundstate%ngridk(:))
        ! q-point (reduced)
            iqr = iqmapr (iv2(1), iv2(2), iv2(3))
        ! q-point (non-reduced)
            iq = iqmap (iv2(1), iv2(2), iv2(3))
            Select Case (trim(input%xs%bse%bsetype))
            Case ('singlet', 'triplet')
           ! read screened Coulomb interaction
!(wol)               Call getbsemat ('SCCLI.OUT', ikkp, nst1, nst3, sccli)
               Call getbsemat ('SCCLI.OUT', ikkp, rnst1, rnst3, sccli)
            End Select
        ! read exchange Coulomb interaction
            Select Case (trim(input%xs%bse%bsetype))
            Case ('rpa', 'singlet')
!(wol)               Call getbsemat ('EXCLI.OUT', ikkp, nst1, nst3, excli)
               Call getbsemat ('EXCLI.OUT', ikkp, rnst1, rnst3, excli)
            End Select
            !Write (*,*) 'nst1, nst2, nst3, nst4', nst1, nst2, nst3, nst4 !(wol)    
            !Write (*,*) 'nstsv', nstsv !(wol) 
            !Write (*,*) 'nstocc0, nstunocc0', nstocc0, nstunocc0 !(wol)   
            !Write (*,*) 'istocc0, istocc, istunocc0, istunocc, isto0, isto,&
            !& istu0, istu', istocc0, istocc, istunocc0, istunocc, isto0, isto,&
            !& istu0, istu !(wol)
            !If (istl3 - nsta2 + 1 .Le. nstocc0) Then !(wol)
            !	Write (*,*) 'within range!' !(wol)
            !End If !(wol)    			
      		Call getdocc (iq, iknr, jknr, nsta1, nsto1, istl3+nsta2-1,& 
      	   & istl3+nsto2-1, docc)
        ! set up matrix
!(wol)            Do ist1 = 1 + nvdif, nst1
            Do ist1 = nsta1, nsto1
!(wol)               Do ist3 = 1, nst3 - ncdif
               Do ist3 = nsta2, nsto2
!(wol)                  Do ist2 = 1 + nvdif, nst2
                  Do ist2 = nsta1, nsto1
!(wol)                     Do ist4 = 1, nst4 - ncdif
                     Do ist4 = nsta2, nsto2
!(wol)                        s1 = hamidx (ist1-nvdif, ist3, iknr, nbfbse, &
!(wol)                       & nafbse)
!(wol)                        s2 = hamidx (ist2-nvdif, ist4, jknr, nbfbse, &
!(wol)                       & nafbse)
                        s1 = hamidx (ist1-nsta1+1, ist3-nsta2+1, iknr, nrnst1, &
                       & nrnst3)
                        s2 = hamidx (ist2-nsta1+1, ist4-nsta2+1, jknr, nrnst2, &
                       & nrnst4)
					    kdocc (s1) = docc (ist1-nsta1+1, ist3-nsta2+1) !(wol)
					   !Write (*,*) 'ist1, ist3, ist2, ist4, kdocc(s1)',&
					   !& ist1, ist3, ist2, ist4, kdocc (s1)
                       !Write (*,*) 'iknr, jknr', iknr, jknr !(wol)
                       !Write (*,*) 'ist1, ist2, ist3, ist4', ist1, ist2, ist3, ist4 !(wol)
                       !Write (*,*) 's1, s2', s1, s2 !(wol)
                     ! add diagonal term
!(wol)                        If (s1 .Eq. s2) Then
!(wol)                           de = evalsv (ist3+istocc, iknr) - evalsv &
!(wol)                          & (ist1, iknr) + input%xs%scissor
!(wol)                           ham (s1, s2) = ham (s1, s2) + de - egap + &
!(wol)                          & bsed
                        If (s1 .Eq. s2) Then
                           de = evalsv (ist3+istl3-1, iknr) - evalsv &
                          & (ist1, iknr) + input%xs%scissor                                                                                                       
                            ham (s1, s2) = ham (s1, s2) + de - egap + &
                          & bsed
!(wol)---partial-start-----------------------------------------------------
!							If (kdocc (s1) .Eq. 0) Then
!								ham (s1,s2) = 0
!								Write (*,*) 'kdocc kill! > s1', s1
!                       	    End If
!(wol)---stop--------------------------------------------------------------
                        End If
						If (kdocc (s1) .Gt. 0 .And. kdocc (s1) .Lt. 2) Then
							Write (*,*) 'kdocc s1', kdocc(s1), s1
						End If
                    ! add exchange term
                        Select Case (trim(input%xs%bse%bsetype))
                        Case ('rpa', 'singlet')
                        	ham (s1, s2) = ham (s1, s2) + 2.0d0 * excli &
                           & (ist1, ist3, ist2, ist4) * &
                           & kdocc (s1) * 0.5
                        End Select
                    ! add correlation term
                        Select Case (trim(input%xs%bse%bsetype))
                        Case ('singlet', 'triplet')
                        	ham (s1, s2) = ham (s1, s2) - sccli (ist1, &
                           & ist3, ist2, ist4) * &
                           & kdocc (s1)	* 0.5						
                        End Select
                     !Write (*,*) 'check' !(wol)
                     End Do
                  End Do
               End Do
            End Do
        ! end loop over (k,kp)-pairs
         End Do
      End Do
	  Write (*,*) 'ikkp, s1, s2', ikkp, s1, s2 !(wol)
      Deallocate (excli, sccli, docc)
      Write (unitout,*)
      Write (unitout, '("Info(bse): invoking Lapack routine ZHEEVX")')
      Write (unitout, '(" size of BSE-Hamiltonian	   : ", i8)') hamsiz
      Write (unitout, '(" number of requested solutions : ", i8)') &
     & input%xs%bse%nexcitmax
  ! allocate eigenvector and eigenvalue arrays
      Allocate (beval(hamsiz), bevec(hamsiz, hamsiz))
  ! set number of excitons
      ne = hamsiz
      Call timesec (ts0)
  ! diagonalize Hamiltonian
      Call bsesoldiag (hamsiz, ne, ham, beval, bevec)
      Call timesec (ts1)
  ! deallocate BSE-Hamiltonian
      Deallocate (ham)
      Write (unitout, '(" timing (in seconds)	   :", f12.3)') ts1 - ts0
  ! number of excitons to consider
      nexc = hamsiz
      Allocate (oszs(nexc), oszsa(nexc), sor(nexc), pmat(hamsiz))
      Allocate (w(input%xs%dosWindow%points), &
     & spectr(input%xs%dosWindow%points))
      Call genwgrid (nwdf, input%xs%dosWindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_real=w)
      Do oct = 1, noptcmp
         oszs (:) = zzero
         Call genfilname (basename='EPSILON', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%tddft%aresdf, filnam=fneps)
         Call genfilname (basename='EXCITON', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%tddft%aresdf, filnam=fnexc)
         Call genfilname (basename='EXCITON_SORTED', tq0=.True., &
        & oc1=oct, oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%tddft%aresdf, filnam=fnexcs)
		 Write (*,*) 'after Hmat' !(wol)
     ! read momentum matrix elements
         Allocate (pm(3, nstsv, nstsv))
         Do iknr = 1, nkptnr
            Call getpmat (iknr, vkl, 1, nstsv, 1, nstsv, .True., 'PMAT_XS.OUT',&
            & pm)
!(wol)            Do ist1 = 1 + nvdif, nstsv - nstunocc0
            Do ist1 = nsta1, nsto1
!(wol)               Do ist2 = nstocc0 + 1, nstsv - ncdif
               Do ist2 = istl3 + nsta2 - 1, istl3 + nsto2 - 1
!(wol)                  s1 = hamidx (ist1-nvdif, ist2-nstocc0, iknr, nbfbse, &
!(wol)                 & nafbse)
                  s1 = hamidx (ist1-nsta1+1, ist2-istl3-nsta2+2,&
                 & iknr, nrnst1, nrnst3)
!(wol)---partial-start-----------------------------------------------------
!				  If (istl3 - nsta2 + 1 .Le. nstocc0 &							
!				 & .And. istl3 + ist2 - 1 .Le. nstocc0) Then
!                  	pmat (s1) = pm (oct, ist1, ist2) * &
!                   & docc (ist1, ist2-istl3-nsta2+2) * 0.5		
!                  Else                 
                  	pmat (s1) = pm (oct, ist1, ist2)
!				  End If
               End Do
            End Do
         End Do
         Deallocate (pm)
		 Write (*,*) 'after Hmat' !(wol)
		 Write (*,*) 'nexc', nexc !(wol)
     ! calculate oscillators for spectrum
         Do s1 = 1, nexc
            Do iknr = 1, nkptnr
!(wol)               Do iv = 1, nbfbse
               Do iv = 1, nrnst1
!(wol)                  Do ic = 1, nafbse
                  Do ic = 1, nrnst3
!(wol)                     s2 = hamidx (iv, ic, iknr, nbfbse, nafbse)
!(wol)                     oszs (s1) = oszs (s1) + bevec (s2, s1) * pmat (s2) &
!(wol)                    & / (evalsv(ic+istocc, iknr)-evalsv(iv+nvdif, &
!(wol)                    & iknr))					
                     s2 = hamidx (iv, ic, iknr, nrnst1, nrnst3)
                     oszs (s1) = oszs (s1) + bevec (s2, s1) * pmat (s2) &
                    & * kdocc (s2) * 0.5 & 
                    & / (evalsv(ic+istl3+nsta2-1, iknr)-evalsv(iv+nsta1-1, &
                    & iknr))
                  End Do
               End Do
            End Do
         End Do
		 Write (*,*) 'after Hmat2' !(wol)
         spectr (:) = zzero
         Do iw = 1, input%xs%dosWindow%points
            Do s1 = 1, nexc
           ! Lorentzian lineshape
               spectr (iw) = spectr (iw) + Abs (oszs(s1)) ** 2 * &
              & (1.d0/(w(iw)-(beval(s1)+egap-bsed)+zi*input%xs%broad))
               If (input%xs%tddft%aresdf) spectr (iw) = spectr (iw) + &
              & Abs (oszs(s1)) ** 2 * &
              & (1.d0/(-w(iw)-(beval(s1)+egap-bsed)-zi*input%xs%broad))
            End Do
         End Do
         spectr (:) = l2int (oct .Eq. oct) * 1.d0 - spectr (:) * 8.d0 * &
        & pi / omega / nkptnr
		write (*,*) 'omega, nkptnr, egap, bsed', omega, nkptnr,egap, bsed !(wol)
     ! write BSE spectrum
         Call writeeps (iqmt, oct, oct, w, spectr, fneps)
     ! oscillator strengths
         Call getunit (unexc)
         Open (unexc, File=fnexc, Form='formatted', Action='write', &
        & Status='replace')
         Do s2 = 1, hamsiz
            Write (unexc, '(i8, 5g18.10)') s2, &
           & (beval(s2)+egap-dble(bsed)) * escale, &
           & (beval(s2)+dble(bsed)) * escale, Abs (oszs(s2))
         End Do
         Write (unexc, '("# Nr.  E		      E - E_gap        |Osc.Str.|")&
        &')
         Write (unexc, '("# E_gap : ", g18.10)') egap * escale
         If (input%xs%tevout) write (unexc, '("# energies are in electr&
        &on volts")')
         Close (unexc)
     ! oscillator strengths sorted
         oszsa = Abs (oszs)
         Call sortidx (hamsiz, oszsa, sor)
         sor = sor (hamsiz:1:-1)
         Open (unexc, File=fnexcs, Form='formatted', Action='write', &
        & Status='replace')
         Do s1 = 1, hamsiz
            s2 = sor (s1)
            Write (unexc, '(i8, 4g18.10)') s1, &
           & (beval(s2)+egap-dble(bsed)) * escale, &
           & (beval(s2)+dble(bsed)) * escale, Abs (oszs(s2))
         End Do
         Write (unexc, '("#	  Nr.	E		  E - E_gap	   |Osc.Str.|")')
         Write (unexc, '("# E_gap : ", g18.10)') egap * escale
         If (input%xs%tevout) write (unexc, '("# energies are in electr&
        &on volts")')
         Close (unexc)
     ! end loop over optical components
      End Do
Contains
!
      Integer Function hamidx (i1, i2, ik, n1, n2)
         Implicit None
         Integer, Intent (In) :: i1, i2, ik, n1, n2
         hamidx = i2 + n2 * (i1-1) + n1 * n2 * (ik-1)
      End Function hamidx
!
End Subroutine bse
!EOC
