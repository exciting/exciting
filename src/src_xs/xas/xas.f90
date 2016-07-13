
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: xas
! !INTERFACE:
Subroutine xas
! !USES:
      Use modinput
      Use modmain
      use modmpi
      Use modxs
      Use m_genwgrid
      Use m_getpmatxas
      Use m_genfilname
      Use m_getunit
      Use m_genloss
      Use m_gensigma
      Use m_gensumrls
      Use m_writeeps
      Use m_writeloss
      Use m_writesigma
      Use m_writesumrls
	  Use modxas
! !DESCRIPTION:
!   Solves the Bethe-Salpeter equation (BSE). The BSE is treated as equivalent
!   effective eigenvalue problem (thanks to the spectral theorem that can
!   be applied to the original BSE in the case of a statically screened Coulomb
!   interaction). The effective BSE-Hamiltonian consists of three parts
!   originating from different sources. It reads
!   $$ H^{\rm eff} = H^{\rm diag} + H^{\rm x} + H^{\rm c}, $$
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
!   $$ H^{\rm eff} = H^{\rm diag} + H^{\rm x}. $$
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
!   Added possibility to compute off-diagonal optical components, Dec 2013 (Stefan Kontur, STK)
!   Adaption to the calculation of core-level excitations, June 2015 (Christian Vorwerk)
!EOP
!BOC
      Implicit None
  ! local variables
      
      ! CC
      integer :: iostat, ievec
      logical :: exist
      Character (256) :: locext

      Integer, Parameter :: iqmt = 0
      Integer, Parameter :: noptcmp = 3
      Real (8), Parameter :: epsortho = 1.d-12
      Character (256) :: fnexc, fnexcs, dotext
      Integer :: iknr, jknr, iqr, iq, iw, iv2 (3), s1, s2, hamsiz, &
     & nexc, ne
      Integer :: unexc, ist1, ist2, ist3, ist4, ikkp, iv, ic, &
     & nvdif, ncdif, optcompt(3), ist, jst, xasspecies, xasatom, icg
      integer :: oct1, oct2, octu, octl
      Integer ::  nrnst1, nrnst2, nrnst3, nrnst4
      Real (8) :: de, egap, ts0, ts1, sumrls(3)
  ! allocatable arrays
      Integer, Allocatable :: sor (:)
      Real (8), Allocatable :: beval (:), w (:), oszsa (:), loss(:, :, :), &
     &  eval0(:,:)
      Complex (8), Allocatable :: excli (:, :, :, :), sccli (:, :, :, &
     & :), ham (:, :)
      Complex (8), Allocatable :: bevec (:, :), pm (:, :, :), pmat (:, :), &
     & oszs (:, :), spectr (:), sigma(:), buf(:,:,:)
  ! external functions
      Integer, External :: l2int

      integer :: Recl, nstsv_
      real(8) :: vkl_(3)
      
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
  ! Get the species for the XAS calculation
	  xasspecies=input%xs%bse%xasspecies
	  xasatom=input%xs%bse%xasatom
  ! initialize the selected ranges of valence/core and conduction states (Note
  ! that core states are always understood to be local orbitals set by the user
  ! and actually treated as valence states)
      sta1 = xasstart                          ! band index for core states(counted from lowest energy eigenvalue) 
      sto1 = xasstop
      sta2 = input%xs%bse%nstlxas(1) !  ---""---  for conduction states (counted from Fermi energy)
      sto2 = input%xs%bse%nstlxas(2)
      nrnst1 = sto1-sta1+1
      nrnst2 = sto1-sta1+1
      nrnst3 = sto2-sta2+1
      nrnst4 = sto2-sta2+1 
  ! use eigenvector files from screening-calculation
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      Call findocclims (iqmt, istocc0, istocc, istunocc0, istunocc, &
     & isto0, isto, istu0, istu)
      input%xs%emattype = 2
      Call ematbdcmbs (input%xs%emattype)
      Write (unitout,*)
      Write (unitout, '("Info(xas): information on number of states:")')
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
  ! size of BSE-Hamiltonian
      hamsiz = nrnst1 * nrnst3 * nkptnr
  ! allocate arrays for Coulomb interactons
      Allocate (sccli(nrnst1, nrnst3, nrnst2, nrnst4))
      Allocate (excli(nrnst1, nrnst3, nrnst2, nrnst4))
  ! allocate BSE-Hamiltonian (large matrix, up to several GB)
      Allocate (ham(hamsiz, hamsiz))
      ham (:, :) = zzero

  ! read KS energies
      Do iknr = 1, nkptnr
        Call getevalsv (vkl(1, iknr), evalsv(1, iknr))
      End Do
       if (associated(input%gw)) then                                         ! GW Part
        ! to KS eigenvalues to use them later for renormalizing PMAT
        allocate(eval0(nstsv,nkptnr))
        eval0(:,:)=evalsv(:,:)
        ! if scissor correction is presented, one should nullify it
        input%xs%scissor=0.0d0
        ! Read QP Fermi energies and eigenvalues from file
        call getevalqp(nkptnr,vkl,evalsv)
        Write(unitout,'("  Quasi particle energies are read from EVALQP.OUT")')
      end if ! GW															 ! End of GW Part	

  ! read mean value of diagonal of direct term
      bsed = 0.d0
      If ((trim(input%xs%bse%bsetype) .Eq. 'singlet') .Or. &
      &   (trim(input%xs%bse%bsetype) .Eq. 'triplet')) Then
         If (input%xs%bse%bsedirsing) Then
            Call getbsediag
            Write (unitout, '("Info(xas): read diagonal of BSE kernel")')
            Write (unitout, '(" mean value : ", 2g18.10)') bsed
         End If
      End If
  ! determine gap
      egap = 1.d8
      Do iknr = 1, nkptnr
         Do ist1 = sta1, sto1
            Do ist3 = sta2, sto2
               egap = Min (egap, evalsv(ist3+istocc, iknr)-ecore(ist1)+input%xs%scissor)
            End Do
         End Do
      End Do
      Write (unitout, '("Info(xas): gap:", g18.10)') egap
      If (egap .Lt. input%groundstate%epspot) Then
         Write (unitout,*)
         Write (unitout, '("Warning(xas): the system has no gap")')
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
              Call getbsemat ('SCCLI.OUT', ikkp, nrnst1, nrnst3, sccli)
            End Select
        ! read exchange Coulomb interaction
            Select Case (trim(input%xs%bse%bsetype))
            Case ('RPA', 'singlet')
               Call getbsemat ('EXCLI.OUT', ikkp, nrnst1, nrnst3, excli)
            End Select
        ! set up matrix
            Do ist1 = sta1, sto1
               Do ist3 = sta2, sto2
                  Do ist2 = sta1, sto1
                     Do ist4 = sta2, sto2
                        s1 = hamidx (ist1-sta1+1, ist3-sta2+1, iknr, nrnst1, nrnst3)
                        s2 = hamidx (ist2-sta1+1, ist4-sta2+1, jknr, nrnst2, nrnst4)
                     ! add diagonal term
                        If (s1 .Eq. s2) Then
                           de = evalsv (ist3+istl3-sta2, iknr) - ecore &
                          & (ist1) + input%xs%scissor                                                                                                    
                            ham (s1, s2) = ham (s1, s2) + de - egap + &
                          & bsed
                        End If
                    ! add exchange term
                        Select Case (trim(input%xs%bse%bsetype))
                        Case ('RPA', 'singlet')
                           ham (s1, s2) = ham (s1, s2) + excli &
                           & (ist1-sta1+1, ist3-sta2+1, ist2-sta1+1, ist4-sta2+1)
                        End Select
                    ! add correlation term
                        Select Case (trim(input%xs%bse%bsetype))
                        Case ('singlet', 'triplet')	
                           ham (s1, s2) = ham (s1, s2) - sccli (ist1-sta1+1, &
                           & ist3-sta2+1, ist2-sta1+1, ist4-sta2+1)
                        End Select
                     End Do
                  End Do
               End Do
            End Do
        ! end loop over (k,kp)-pairs
         End Do
      End Do
      Deallocate (excli, sccli)
      Write (unitout,*)
      Write (unitout, '("Info(xas): invoking Lapack routine ZHEEVX")')
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
      Allocate (oszs(nexc, 3), oszsa(nexc), sor(nexc), pmat(hamsiz, 3))
      Allocate (w(input%xs%energywindow%points), spectr(input%xs%energywindow%points))
      Allocate (buf(3,3,input%xs%energywindow%points))
      Allocate (loss(3, 3, input%xs%energywindow%points), sigma(input%xs%energywindow%points))
      Call genwgrid (input%xs%energywindow%points, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_real=w)
      buf(:,:,:)=zzero

      do oct1 = 1, noptcmp
     ! read momentum matrix elements
         Allocate (pm(3, ncg, nstsv))
         Do iknr = 1, nkptnr
            Call getpmatxas (iknr, vkl, 1, ncg, 1, nstsv, .True., &
           & 'PMAT_XS.OUT', pm)
            Do ist1 = sta1, sto1
               Do ist2 = istl3, istl3 + nrnst3 - 1
! DIN: Renormalize pm according to DelSole PRB48, 11789 (1993)
                  if (associated(input%gw)) then
                     pm(1:3,ist1,ist2)=pm(1:3,ist1,ist2)  * &
                    &  (evalsv(ist2,iknr)-evalsv(ist1,iknr)) / &
                    &   (eval0(ist2,iknr)- eval0(ist1,iknr))
                  end if 
! DIN
                  s1 = hamidx (ist1-sta1+1, ist2-istl3+1,&
                 & iknr, nrnst1, nrnst3)    
                  pmat (s1, oct1) = pm (oct1, ist1, ist2)
               End Do
            End Do
         End Do
         Deallocate (pm)

     ! calculate oscillators for spectrum
         oszs (:, oct1) = zzero
         Do s1 = 1, nexc
            Do iknr = 1, nkptnr
               Do iv = 1, nrnst1
                  Do ic = 1, nrnst3
                     s2 = hamidx (iv, ic, iknr, nrnst1, nrnst3)
                     oszs (s1, oct1) = oszs (s1, oct1) + bevec (s2, s1) * pmat (s2, oct1) &
                    &  / (evalsv(ic+istl3-1, iknr)- &
                    & ecore(iv+sta1-1))
                  End Do
               End Do
            End Do
         End Do
     ! STK: add case of double grid
         if (dgrid) then 
            Write (dotext, '("_SG", I3.3, ".OUT")') iksubpt
            Call genfilname (basename='EXCITON', tq0=.True., oc1=oct1, &
        &    oc2=oct1, bsetype=input%xs%bse%bsetype, &
        &    scrtype=input%xs%screening%screentype, nar= .Not. &
        &    input%xs%bse%aresbse, dotext=dotext, filnam=fnexc)
            Call genfilname (basename='EXCITON_SORTED', tq0=.True., &
        &    oc1=oct1, oc2=oct1, bsetype=input%xs%bse%bsetype, &
        &    scrtype=input%xs%screening%screentype, nar= .Not. &
        &    input%xs%bse%aresbse, dotext=dotext, filnam=fnexcs)
         else
            Call genfilname (basename='EXCITON', tq0=.True., oc1=oct1, &
        &    oc2=oct1, bsetype=input%xs%bse%bsetype, &
        &    scrtype=input%xs%screening%screentype, nar= .Not. &
        &    input%xs%bse%aresbse, filnam=fnexc)
            Call genfilname (basename='EXCITON_SORTED', tq0=.True., &
        &    oc1=oct1, oc2=oct1, bsetype=input%xs%bse%bsetype, &
        &    scrtype=input%xs%screening%screentype, nar= .Not. &
        &    input%xs%bse%aresbse, filnam=fnexcs)
         endif
     ! oscillator strengths
         Call getunit (unexc)
         Open (unexc, File=fnexc, Form='formatted', Action='write', &
        & Status='replace')
         Do s2 = 1, hamsiz
            Write (unexc, '(i8, 6g18.10)') s2, &
           & (beval(s2)+egap-dble(bsed)) * escale, &
           & (beval(s2)+dble(bsed)) * escale, abs(oszs(s2, oct1)), dble(oszs(s2, oct1)), aimag(oszs(s2, oct1))
         End Do
         Write (unexc, '("# Nr.  E		      E - E_gap        |Osc.Str.|      Re      Im")&
        &')
         Write (unexc, '("# E_gap : ", g18.10)') egap * escale
         If (input%xs%tevout) write (unexc, '("# energies are in electron volts")')
         Close (unexc)
     ! oscillator strengths sorted
         oszsa(:) = Abs (oszs(:, oct1))
         Call sortidx (hamsiz, oszsa, sor)
         sor = sor (hamsiz:1:-1)
         Open (unexc, File=fnexcs, Form='formatted', Action='write', &
        & Status='replace')
         Do s1 = 1, hamsiz
            s2 = sor (s1)
            Write (unexc, '(i8, 4g18.10)') s1, &
           & (beval(s2)+egap-dble(bsed)) * escale, &
           & (beval(s2)+dble(bsed)) * escale, Abs (oszs(s2, oct1))
         End Do
         Write (unexc, '("#	  Nr.	E		  E - E_gap	   |Osc.Str.|")')
         Write (unexc, '("# E_gap : ", g18.10)') egap * escale
         If (input%xs%tevout) write (unexc, '("# energies are in electr&
        &on volts")')
         Close (unexc)
      enddo

! STK: if run is within a double grid loop stop here
      if (dgrid) goto 10

      Do oct1 = 1, noptcmp
! STK: compute off-diagonal optical components if requested
       If (input%xs%dfoffdiag) Then
             octl = 1
             octu = noptcmp
       Else
             octl = oct1
             octu = oct1
       End If
        Do oct2 = octl, octu
         optcompt = (/ oct1, oct2, 0 /)
         spectr (:) = zzero
         Do iw = 1, input%xs%energywindow%points
            Do s1 = 1, nexc
           ! Lorentzian lineshape
               spectr (iw) = spectr (iw) + & !Abs (oszs(s1)) ** 2 * &
! STK
              & oszs(s1, oct1) * conjg(oszs(s1, oct2)) * &
              & (1.d0/(w(iw)-(beval(s1)+egap-bsed)+zi*input%xs%broad))             
               If (input%xs%bse%aresbse) spectr (iw) = spectr (iw) + &
! STK
              & oszs(s1, oct1) * conjg(oszs(s1, oct2)) * &
!             & Abs (oszs(s1)) ** 2 * &
              & (1.d0/(-w(iw)-(beval(s1)+egap-bsed)-zi*input%xs%broad))
            End Do
         End Do
! STK
         spectr (:) = l2int (oct1 .Eq. oct2) * 1.d0 - spectr (:) * 8.d0 * &
        & pi / omega / nkptnr
         buf(oct1,oct2,:)=spectr(:)
     ! end loops over optical components
        enddo
      End Do
! STK
      call genloss (buf, loss, noptcmp)
!
      Do oct1 = 1, noptcmp
       If (input%xs%dfoffdiag) Then
             octl = 1
             octu = noptcmp
       Else
             octl = oct1
             octu = oct1
       End If
        Do oct2 = octl, octu
         optcompt = (/ oct1, oct2, 0 /)
         Call genfilname (basename='EPSILON', tq0=.True., oc1=oct1, &
        & oc2=oct2, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fneps)
         Call genfilname (basename='LOSS', tq0=.True., oc1=oct1, &
        & oc2=oct2, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnloss)
         Call genfilname (basename='SIGMA', tq0=.True., oc1=oct1, &
        & oc2=oct2, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnsigma)
         Call genfilname (basename='SUMRULES', tq0=.True., oc1=oct1, &
        & oc2=oct2, bsetype=input%xs%bse%bsetype, &
        & scrtype=input%xs%screening%screentype, nar= .Not. &
        & input%xs%bse%aresbse, filnam=fnsumrules)
     ! symmetrize the macroscopic dielectric function tensor
         Call symt2app (oct1, oct2, input%xs%energywindow%points, symt2, buf, spectr)
     ! generate optical functions
! STK
!        Call genloss (spectr, loss)
         Call gensigma (w, spectr, optcompt, sigma)
         Call gensumrls (w, spectr, sumrls)
     ! write optical functions to file
         Call writeeps (iq, oct1, oct2, w, spectr, trim(fneps))
         Call writeloss (iq, w, loss(oct1, oct2, :), trim(fnloss))
         Call writesigma (iq, w, sigma, trim(fnsigma))
         Call writesumrls (iq, sumrls, trim(fnsumrules))
        enddo
      end do

      if (associated(input%xs%storeexcitons)) then
        !-----------------------------------------------------
        ! upon request, store array with exciton coefficients
        !-----------------------------------------------------
        if ( (input%xs%storeexcitons%MinNumberExcitons .lt. 1) .or. &
          &  (input%xs%storeexcitons%MinNumberExcitons .gt. hamsiz) .or. &
          &  (input%xs%storeexcitons%MaxNumberExcitons .lt. 1) .or. &
          &  (input%xs%storeexcitons%MaxNumberExcitons .gt. hamsiz) .or. &
          &  (input%xs%storeexcitons%MinNumberExcitons .gt. input%xs%storeexcitons%MaxNumberExcitons) ) then
          write(*,*)
          write(*,'("Error(bse): wrong range of exciton indices: ", 2I5)') &
                & input%xs%storeexcitons%MinNumberExcitons, input%xs%storeexcitons%MaxNumberExcitons
          write(*,*)
          stop
        end if  
        ! write bin
        open(50,File='EXCCOEFF.bin', & 
             Action='WRITE',Form='UNFORMATTED', IOstat=iostat)
        if ((iostat/=0) .and. (rank==0)) then
          write(*,*) iostat
          write(*,'("Error(bse): error creating EXCCOEFF.bin")')
          write(*,*)
          stop
        end if
        ! write
        write(50) input%xs%storeexcitons%MinNumberExcitons, input%xs%storeexcitons%MaxNumberExcitons, & 
        &         nkptnr, istl3, sta1, sta2, nrnst1, nrnst3, hamsiz
        do ievec = input%xs%storeexcitons%MinNumberExcitons, input%xs%storeexcitons%MaxNumberExcitons
           write(50) beval(ievec), bevec(1:hamsiz,ievec)
        end do
        close(50)
      end if

      deallocate(beval,bevec,oszs,oszsa,sor,pmat,w,spectr,loss,sigma,buf)
      if (associated(input%gw)) deallocate(eval0)

      10 continue
      call barrier

Contains
!
      Integer Function hamidx (i1, i2, ik, n1, n2)
         Implicit None
         Integer, Intent (In) :: i1, i2, ik, n1, n2
         hamidx = i2 + n2 * (i1-1) + n1 * n2 * (ik-1)
      End Function hamidx

End Subroutine xas
!EOC

