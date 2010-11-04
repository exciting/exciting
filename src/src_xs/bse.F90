
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bse
! !INTERFACE:
Subroutine bse
! !USES:
      Use modinput
      Use modmain
      use modmpi
      Use modxs
      Use m_genwgrid
      Use m_getpmat
      Use m_genfilname
      Use m_getunit
      Use m_genloss
      Use m_gensigma
      Use m_gensumrls
      Use m_writeeps
      Use m_writeloss
      Use m_writesigma
      Use m_writesumrls
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
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!      level for the treatment of core excitations (using local orbitals).
!      October 2010 (Weine Olovsson)
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
     & nvdif, ncdif, optcompt(3), ist, jst
      Integer :: sta1, sto1, sta2, sto2, nsta1, nsto1, nsta2, nsto2, &
     & rnst1, rnst2, rnst3, rnst4, nrnst1, nrnst2, nrnst3, nrnst4
      Real (8) :: de, egap, ts0, ts1, sumrls(3)
  ! allocatable arrays
      Integer, Allocatable :: sor (:)
      Real (8), Allocatable :: beval (:), w (:), oszsa (:), loss(:), &
     & docc(:,:), kdocc (:)
      Complex (8), Allocatable :: excli (:, :, :, :), sccli (:, :, :, &
     & :), ham (:, :)
      Complex (8), Allocatable :: bevec (:, :), pm (:, :, :), pmat (:), &
     & oszs (:), spectr (:), sigma(:), buf(:,:,:)
  ! external functions
      Integer, External :: l2int
  ! routine not yet parallelized
  if (rank .ne. 0) goto 10
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
  ! initialize the selected ranges of valence/core and conduction states (Note
  ! that core states are always understood to be local orbitals set by the user
  ! and actually treated as valence states)
      sta1 = input%xs%bse%nstlbsemat(1)
      sto1 = input%xs%bse%nstlbsemat(2)
      sta2 = input%xs%bse%nstlbsemat(3)
      sto2 = input%xs%bse%nstlbsemat(4)
      nsta1 = input%xs%bse%nstlbse(1)
      nsto1 = input%xs%bse%nstlbse(2)
      nsta2 = input%xs%bse%nstlbse(3)
      nsto2 = input%xs%bse%nstlbse(4)
      rnst1 = sto1-sta1+1
      rnst2 = sto1-sta1+1
      rnst3 = sto2-sta2+1
      rnst4 = sto2-sta2+1
      nrnst1 = nsto1-nsta1+1
      nrnst2 = nsto1-nsta1+1
      nrnst3 = nsto2-nsta2+1
      nrnst4 = nsto2-nsta2+1 
  ! use eigenvector files from screening-calculation
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      Call findocclims (iqmt, istocc0, istocc, istunocc0, istunocc, &
     & isto0, isto, istu0, istu)
      input%xs%emattype = 2
      Call ematbdcmbs (input%xs%emattype)
      Write (unitout,*)
      Write (unitout, '("Info(bse): information on number of states:")')
      Write (unitout, '(" number of states below Fermi energy in Hamilt&
     &onian:", i6)') nrnst1
      Write (unitout, '(" number of states above Fermi energy in Hamilt&
     &onian:", i6)') nrnst3
      Write (unitout, '(" ranges of states according to BSE matrix:")')
      Write (unitout, '("  range of first index and number  :", 2i6, 3x&
     &, i6)') istl1, istu1, nst1
      Write (unitout, '("  range of second index and number :", 2i6, 3x&
     &, i6)') istl2, istu2, nst2
      Write (unitout, '("  range of third index and number  :", 2i6, 3x&
     &, i6)') istl3, istu3, nst3
      Write (unitout, '("  range of fourth index and number :", 2i6, 3x&
     &, i6)') istl4, istu4, nst4
      If ((nsta1 .Lt. sta1) .Or. (nsto1 .Gt. sto1) .Or. (nsta2 .Lt. sta2)&
     & .Or. (nsto2 .Gt. sto2)) Then
         Write (unitout,*)
         Write (unitout, '("Error(bse): inconsistency in ranges of stat&
        &es - nstlbse must be within the range of nstlbsemat")')
         Write (unitout,*)
         Call terminate
      End If
  ! size of BSE-Hamiltonian
      hamsiz = nrnst1 * nrnst3 * nkptnr
  ! allocate arrays for Coulomb interactons
      Allocate (sccli(rnst1, rnst3, rnst2, rnst4))
      Allocate (excli(rnst1, rnst3, rnst2, rnst4))
  ! allocate array for occupation number difference (future use )
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
         Do ist1 = nsta1, nsto1
            Do ist3 = nsta2, nsto2
               egap = Min (egap, evalsv(ist3+istocc, iknr)-evalsv(ist1, &
              & iknr)+input%xs%scissor)
            End Do
         End Do
      End Do
      Write (unitout, '("Info(bse): gap:", g18.10)') egap
      If (egap .Lt. input%groundstate%epspot) Then
         Write (unitout,*)
         Write (unitout, '("Warning(bse): the system has no gap")')
         Write (unitout,*)
      End If	  
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
               Call getbsemat ('SCCLI.OUT', ikkp, rnst1, rnst3, sccli)
            End Select
        ! read exchange Coulomb interaction
            Select Case (trim(input%xs%bse%bsetype))
            Case ('RPA', 'singlet')
               Call getbsemat ('EXCLI.OUT', ikkp, rnst1, rnst3, excli)
            End Select
        ! get occupation numbers (future use for partial occupancy) 			
      		Call getdocc (iq, iknr, jknr, nsta1, nsto1, istl3+nsta2-1,& 
      	   & istl3+nsto2-1, docc)
        ! set up matrix
            Do ist1 = nsta1, nsto1
               Do ist3 = nsta2, nsto2
                  Do ist2 = nsta1, nsto1
                     Do ist4 = nsta2, nsto2
                        s1 = hamidx (ist1-nsta1+1, ist3-nsta2+1, iknr, nrnst1, &
                       & nrnst3)
                        s2 = hamidx (ist2-nsta1+1, ist4-nsta2+1, jknr, nrnst2, &
                       & nrnst4)
					    kdocc (s1) = docc (ist1-nsta1+1, ist3-nsta2+1)
                     ! add diagonal term
                        If (s1 .Eq. s2) Then
                           de = evalsv (ist3+istl3-1, iknr) - evalsv &
                          & (ist1, iknr) + input%xs%scissor                                                                                                       
                            ham (s1, s2) = ham (s1, s2) + de - egap + &
                          & bsed
                        End If
                    ! write partially occupied states in the output    
						If (kdocc (s1) .Gt. input%groundstate%epsocc .And. kdocc (s1) &
						  .Lt. 2.d0-input%groundstate%epsocc) Then
							Write (*,*) 'kdocc s1', kdocc(s1), s1
						End If
                    ! add exchange term
                        Select Case (trim(input%xs%bse%bsetype))
                        Case ('RPA', 'singlet')
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
                     End Do
                  End Do
               End Do
            End Do
        ! end loop over (k,kp)-pairs
         End Do
      End Do
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
      Allocate (w(input%xs%energywindow%points), spectr(input%xs%energywindow%points))
      Allocate (buf(3,3,input%xs%energywindow%points))
      Allocate (loss(input%xs%energywindow%points), sigma(input%xs%energywindow%points))
      Call genwgrid (input%xs%energywindow%points, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_real=w)
     buf(:,:,:)=zzero
      Do oct = 1, noptcmp
         optcompt = (/ oct, oct, 0 /)
         oszs (:) = zzero
         Call genfilname (basename='EXCITON', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnexc)
         Call genfilname (basename='EXCITON_SORTED', tq0=.True., &
        & oc1=oct, oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnexcs)
     ! read momentum matrix elements
         Allocate (pm(3, nstsv, nstsv))
         Do iknr = 1, nkptnr
            Call getpmat (iknr, vkl, 1, nstsv, 1, nstsv, .True., 'PMAT_XS.OUT',&
            & pm)
            Do ist1 = nsta1, nsto1
               Do ist2 = istl3 + nsta2 - 1, istl3 + nsto2 - 1
                  s1 = hamidx (ist1-nsta1+1, ist2-istl3-nsta2+2,&
                 & iknr, nrnst1, nrnst3)              
                  	pmat (s1) = pm (oct, ist1, ist2)
               End Do
            End Do
         End Do
         Deallocate (pm)
     ! calculate oscillators for spectrum
         Do s1 = 1, nexc
            Do iknr = 1, nkptnr
               Do iv = 1, nrnst1
                  Do ic = 1, nrnst3				
                     s2 = hamidx (iv, ic, iknr, nrnst1, nrnst3)
                     oszs (s1) = oszs (s1) + bevec (s2, s1) * pmat (s2) &
                    & * kdocc (s2) * 0.5 & 
                    & / (evalsv(ic+istl3-1+nsta2-1, iknr)- &
                    evalsv(iv+nsta1-1, iknr))
                  End Do
               End Do
            End Do
         End Do
         spectr (:) = zzero
         Do iw = 1, input%xs%energywindow%points
            Do s1 = 1, nexc
           ! Lorentzian lineshape
               spectr (iw) = spectr (iw) + Abs (oszs(s1)) ** 2 * &
              & (1.d0/(w(iw)-(beval(s1)+egap-bsed)+zi*input%xs%broad))             
               If (input%xs%bse%aresbse) spectr (iw) = spectr (iw) + &
              & Abs (oszs(s1)) ** 2 * &
              & (1.d0/(-w(iw)-(beval(s1)+egap-bsed)-zi*input%xs%broad))
            End Do
         End Do
         spectr (:) = l2int (oct .Eq. oct) * 1.d0 - spectr (:) * 8.d0 * &
        & pi / omega / nkptnr
         buf(oct,oct,:)=spectr(:)
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
      Do oct = 1, noptcmp
         optcompt = (/ oct, oct, 0 /)
         Call genfilname (basename='EPSILON', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fneps)
         Call genfilname (basename='LOSS', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnloss)
         Call genfilname (basename='SIGMA', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnsigma)
         Call genfilname (basename='SUMRULES', tq0=.True., oc1=oct, &
        & oc2=oct, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnsumrules)
     ! symmetrize the macroscopic dielectric function tensor
         Call symt2app (oct, oct, input%xs%energywindow%points, symt2, buf, spectr)
     ! generate optical functions
         Call genloss (spectr, loss)
         Call gensigma (w, spectr, optcompt, sigma)
         Call gensumrls (w, spectr, sumrls)
     ! write optical functions to file
         Call writeeps (iq, oct, oct, w, spectr, trim(fneps))
         Call writeloss (iq, w, loss, trim(fnloss))
         Call writesigma (iq, w, sigma, trim(fnsigma))
         Call writesumrls (iq, sumrls, trim(fnsumrules))
      end do
      deallocate(beval,bevec,oszs,oszsa,sor,pmat,w,spectr,loss,sigma,buf)
10 continue
      call barrier      
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

