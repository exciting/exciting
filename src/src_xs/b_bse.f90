! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_bse
! !INTERFACE:
subroutine b_bse
! !USES:
  use modinput, only: input
  use mod_constants, only: zzero, h2ev
  use mod_kpoint, only: nkptnr, vkl
  use mod_eigenvalue_occupancy, only: evalsv, nstsv
  use modmpi, only: rank, barrier
  use modxs, only: unitout, bcbs, bsed
  use modbse
  use m_genwgrid
  use m_genfilname
  use m_diagfull
  use m_writeoscillator
  use m_writecmplxparts
  use m_dhesolver
  use m_dzgemm
  use modsclbse
! !DESCRIPTION:
!   Solves the Bethe-Salpeter equation(BSE). The BSE is treated as equivalent
!   effective eigenvalue problem(thanks to the spectral theorem that can
!   be applied to the original BSE in the case of a statically screened Coulomb
!   interaction). The effective BSE-Hamiltonian consists of three parts
!   originating from different sources. It reads
!   $$ H^{\rm eff} = H^{\rm diag} + 2H^{\rm x} + H^{\rm c}, $$
!   where $H^{\rm diag}$ is the diagonal part stemming from the independent
!   particle transitions, $H^{\rm x}$ denotes the exchange-term caused by the
!   unscreened(bare) Coulomb interaction, whereas $H^{\rm c}$ accounts for the
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
!   For the diagonalization of the Hamiltonian, a LAPACK-routine({\tt zheevx})
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
!   The macroscopic dielectric function(MDF) is obtained by the realation
!   $$ {\rm Im}\; \epsilon^{i}_{\rm M}(\omega) = \frac{8\pi^2}{V}
!                     \sum_{\lambda} t^{i}_{\lambda}
!                     \delta(\omega-\varepsilon_{\lambda}+\Delta),$$
!   where $\epsilon^{i}_{\rm M}$ is the MDF for the $i$-th polarization, $V$
!   denotes the crystal volume and $\Delta$ is a constant shift of the
!   conduction bands(scissors shift). The delta-function in the latter
!   expression is convoluted with a(symmetrized) Lorentzian
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
!   Created June 2008(S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!      level for the treatment of core excitations(using local orbitals).
!      October 2010(Weine Olovsson)
!   Added possibility to compute off-diagonal optical components, Dec 2013(Stefan Kontur, STK)
!   Forked from bse.F90 and adapted for non TDA BSE. (Aurich)
!EOP
!BOC

  implicit none

  ! Local variables
  ! Parameters
  integer(4), parameter :: iqmt = 0
  ! Variables
  integer(4) :: iknr, iq, nw
  integer(4) :: hamsize, nexc
  real(8) :: ts0, ts1
  logical :: fcoup, fwp, fscal

  ! Allocatable arrays
  real(8), allocatable, dimension(:) :: bevalim, bevalre, w
  complex(8), allocatable, dimension(:,:) :: ham, bevecr, oszsr, oszsa
  complex(8), allocatable, dimension(:,:,:) :: symspectr
  ! Distributed arrays
  integer(4) :: ip
  complex(8), allocatable :: buff(:,:)
  type(dzmat) :: dham, dbevecr, doszsr


  ! Include coupling terms
  fcoup = input%xs%bse%coupling

  ! Use ScaLAPACK
  fscal = input%xs%bse%distribute

  if(fscal .and. fcoup) then 
    write(*,*) "Coupling and scalapack not possible"
    call terminate
  end if

  ! Non-parallelized code.
  if(.not. fscal .and. rank == 0) then 

    ! General init
    call init0
    ! K-grid init
    call init1
    ! Q-grid init
    call init2
    ! Saves all k grid related variables in modxs, so 
    ! that another k grid can be used for k+q.
    call xssave0

    ! Read Fermi energy from file
    call readfermi

    !! Initialize the selected ranges of valence/core and conduction states(note
    !! that core states are always understood to be local orbitals set by the user
    !! and actually treated as valence states)
    ! Use eigenvector files from screening-calculation 
    ! and set GS occupation limits in modxs (read from *_SCR.OUT)
    call genfilname(dotext='_SCR.OUT', setfilext=.true.)

    ! Set ist* variables in modxs using findocclims
    call setranges_modxs(iqmt)

    ! Set band combinations (modbse:bcou & modbse:bcouabs)
    call setbcbs_bse
    ! Number of kkp combinations with ikp >= ik (modbse:nk & modbse:nkkp)
    nk   = nkptnr
    nkkp = nkptnr*(nkptnr+1)/2

    ! Read KS energies (form EVALSV_SCR.OUT) and save them to modxs
    if(allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,nkptnr))
    do iknr = 1, nkptnr
      call getevalsv(vkl(1:3, iknr), evalsv(1:nstsv, iknr))
    end do

    ! If on top of GW
    if(associated(input%gw)) then
      ! Save KS eigenvalues to use them later for renormalizing PMAT
      allocate(eval0(nstsv,nkptnr))
      eval0(:,:)=evalsv(:,:)
      ! If scissor correction is presented, one should nullify it
      input%xs%scissor=0.0d0
      ! Read QP Fermi energies and eigenvalues from file
      call getevalqp(nkptnr,vkl,evalsv)
      ! evalsv = evalsv+efermi
      write(unitout,'("  Quasi particle energies are read from EVALQP.OUT")')
    end if

    ! Read mean value of diagonal of direct term
    bsed = 0.d0
    if((trim(input%xs%bse%bsetype) .eq. 'singlet')&
      & .or. (trim(input%xs%bse%bsetype) .eq. 'triplet')) then

      ! False by default
      if(input%xs%bse%bsedirsing) then
        call getbsediag
        write(unitout, '("Info(bse): read diagonal of BSE kernel")')
        write(unitout, '(" mean value : ", 2g18.10)') bsed
      end if

    end if

    ! Determine the minimal optical gap (q=0), w.r.t. the selected bands
    egap = minval(evalsv(bcouabs%il2,1:nkptnr) - evalsv(bcouabs%iu1,1:nkptnr)&
      &+ input%xs%scissor)
    write(unitout,*)
    write(unitout, '("Info(bse): gap:", E23.16)') egap
    ! Warn if the system has no gap even with scissor (or no scissor and on top of GW)
    if(egap .lt. input%groundstate%epspot) then
      write(unitout,*)
      write(unitout, '("Warning(bse): the system has no gap, setting it to 0")')
      write(unitout,*)
      egap = 0.0d0
    end if  

    !!<-- Inspecting occupancies 
    ! Do not use state combination with zero occupancy difference 
    ! (Theory not complete here?). Also do not use state combinations
    ! where the "unoccupied" state has larger occupancy, since hat will
    ! make the BSE Hamiltonian non-hermitian, even in TDA and q=0.
    !
    ! Set gamma point (currently iq does nothing in checkoccupancies)
    iq = 1
    ! Get adjusted hamsize, index maps and occupation factors.
    ! Also write information about skipped combinations to file.
    call checkoccupancies(iq, hamsize)
    !!-->

    fwp = input%xs%bse%writeparts

    ! Write Info
    if(rank == 0) then
      write(unitout,*)
      write(unitout, '("Info(bse): Assembling BSE matrix")')
      write(unitout, '("RR/RA blocks of BSE-Hamiltonian: Shape=",i4, " Size=",i8)')&
        & hamsize, hamsize**2
      if(fcoup) then
        write(unitout, '(" Including coupling terms ")')
        write(unitout, '(" Full BSE-Hamiltonian: Shape=",i4, " Size=",i8)')&
          & 2*hamsize, 4*hamsize**2
      end if
      if(fwp) then
        write(unitout, '(" Writing real and imaginary parts of Hamiltonian to file ")')
      end if
    end if
    ! Assemble Hamiltonian matrix 
    call timesec(ts0)
    call setuphamiltonian(ham, fwp)
    call timesec(ts1)
    write(unitout, '(" Matrix build.")')
    write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0

    ! Number of excitons to consider
    nexc = input%xs%bse%nexc
    if(nexc > hamsize .or. nexc < 1) then
      nexc = hamsize
    end if
    if(fcoup) then
      ! Use all solutions, since zgeevx produces all anyways.
      nexc = 2*hamsize
    end if

    ! Write Info
    if(rank == 0) then
      write(unitout,*)
      if(fcoup) then
        write(unitout, '("Info(bse): Diagonalizing full non symmetric Hamiltonian")')
        write(unitout, '("Info(bse): Invoking lapack routine ZGEEVX")')
      else
        write(unitout, '("Info(bse): Diagonalizing RR Hamiltonian (TDA)")')
        write(unitout, '("Info(bse): Invoking lapack routine ZHEEVR")')
        write(unitout, '("  Number of requested solutions : ", i8)') nexc
      end if
    end if

    ! Allocate eigenvector and eigenvalue arrays
    if(fcoup) then 
      allocate(bevalre(2*hamsize), bevecr(2*hamsize, 2*hamsize))
      allocate(bevalim(2*hamsize))
    else
      allocate(bevalre(hamsize), bevecr(hamsize, nexc))
    end if

    ! Test write out Hamiltonian
    if(fwp) then
      call writecmplxparts('Ham', dble(ham), immat=aimag(ham))
    end if

    ! Diagonalize Hamiltonian (destroys the content of ham)
    call timesec(ts0)
    if(fcoup) then
      call diagfull(2*hamsize, ham, bevalre,&
        & evalim=bevalim, evecr=bevecr, fbalance=.false., frcond=.false.)
    else
      call hermitian_solver(nexc, hamsize, ham, bevalre, bevecr)
    end if
    call timesec(ts1)

    ! Test write out right-eigenvectors
    if(fwp) then
      call writecmplxparts('bevecr', dble(bevecr), immat=aimag(bevecr))
    end if

    ! Deallocate BSE-Hamiltonian
    deallocate(ham)

    write(unitout, '("  Eigen solutions found.")')
    write(unitout, '("  Timing (in seconds)	   :", f12.3)') ts1 - ts0
    write(unitout,*)

    ! Calculate oscillator strengths.
    allocate(oszsr(nexc,3))
    if(fcoup) then
      allocate(oszsa(nexc,3))
    end if

    write(unitout, '("Making oszillator strengths.")')
    call timesec(ts0)
    if(fcoup) then 
      call makeoszillatorstrength(oszsr, oszstra=oszsa)
    else
      call makeoszillatorstrength(oszsr)
    end if
    call timesec(ts1)
    write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0

    ! Write excition energies and oscillator strengths to 
    ! text file. 
    write(unitout, '("Writing excition energies and oszillator strengths.")')
    if(fcoup) then
      call writeoscillator(2*hamsize, nexc, egap, bevalre, oszsr,&
        & evalim=bevalim, oszstra=oszsa)
    else
      call writeoscillator(hamsize, nexc, egap, bevalre, oszsr)
    end if

    ! Calculate macroscopic dielectric tensor with finite broadening, i.e. "the spectrum"
    nw = input%xs%energywindow%points

    ! Allocate arrays used in spectrum construction
    allocate(w(nw))
    allocate(symspectr(3,3,nw))

    ! Generate an evenly spaced frequency grid 
    call genwgrid(nw, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)

    write(unitout, '("Making spectrum.")')
    call timesec(ts0)
    ! Calculate lattice symmetrized spectrum.
    if(fcoup) then 
      write(unitout, '("  Using general formula.")')
      call makespectrum(nw, w, symspectr)
    else
      write(unitout, '("  Using TDA formula.")')
      call makespectrum_tda(nw, w, symspectr)
    end if
    call timesec(ts1)
    write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0

    write(unitout, '("Writing derived quantities.")')
    ! Generate and write derived optical quantities
    iq = 1
    call writederived(iq, symspectr, nw, w)
    write(unitout, '("  Derived quantities written.")')

    ! Store excitonic energies and wave functions to file
    if(associated(input%xs%storeexcitons)) then
      if(fcoup) then
        call storeexcitons(2*hamsize,nexc,nkptnr,iuref,bcou,smap,bevalre,bevecr)
      else
        call storeexcitons(hamsize,nexc,nkptnr,iuref,bcou,smap,bevalre,bevecr)
      end if
    end if

    ! Clean up
    deallocate(bevalre, bevecr, oszsr, w, symspectr, evalsv)
    if(fcoup) then 
      deallocate(bevalim, oszsa)
    end if
    if(associated(input%gw)) deallocate(eval0)

    call barrier

  ! Parallel version 
  else if (fscal) then

    ! Set up process girds for BLACS 
    call setupblacs

    ! Only MPI ranks that have an associated 
    ! BLACS grid process do anything.
    if( .not. myprow1d_c < 0 ) then 

      ! General init
      call init0
      ! K-grid init
      call init1
      ! Q-grid init
      call init2
      ! Saves all k grid related variables in modxs, so 
      ! that another k grid can be used for k+q.
      call xssave0

      ! Read Fermi energy from file
      call readfermi

      !! Initialize the selected ranges of valence/core and conduction states(note
      !! that core states are always understood to be local orbitals set by the user
      !! and actually treated as valence states)
      ! Use eigenvector files from screening-calculation 
      ! and set GS occupation limits in modxs (read from *_SCR.OUT)
      call genfilname(dotext='_SCR.OUT', setfilext=.true.)

      ! Set ist* variables in modxs using findocclims
      call setranges_modxs(iqmt)

      ! Set band combinations (modbse:bcou & modbse:bcouabs)
      call setbcbs_bse

      ! Number of kkp combinations with ikp >= ik (modbse:nk & modbse:nkkp)
      nk = nkptnr
      nkkp = nkptnr*(nkptnr+1)/2

      ! Read KS energies (form EVALSV_SCR.OUT) and save them to modxs
      if(allocated(evalsv)) deallocate(evalsv)
      allocate(evalsv(nstsv,nkptnr))
      do iknr = 1, nkptnr
        call getevalsv(vkl(1:3, iknr), evalsv(1:nstsv, iknr))
      end do

      ! If on top of GW
      if(associated(input%gw)) then
        ! Save KS eigenvalues to use them later for renormalizing PMAT
        allocate(eval0(nstsv,nkptnr))
        eval0(:,:)=evalsv(:,:)
        ! If scissor correction is presented, one should nullify it
        input%xs%scissor=0.0d0
        ! Read QP Fermi energies and eigenvalues from file
        call getevalqp(nkptnr,vkl,evalsv)
        ! evalsv = evalsv+efermi
        write(unitout,'("  Quasi particle energies are read from EVALQP.OUT")')
      end if

      ! Read mean value of diagonal of direct term
      bsed = 0.d0
      if((trim(input%xs%bse%bsetype) .eq. 'singlet')&
        & .or. (trim(input%xs%bse%bsetype) .eq. 'triplet')) then

        ! False by default
        if(input%xs%bse%bsedirsing) then
          call getbsediag
          write(unitout, '("Info(bse): read diagonal of BSE kernel")')
          write(unitout, '(" mean value : ", 2g18.10)') bsed
        end if

      end if

      ! Determine the minimal optical gap (q=0), w.r.t. the selected bands
      egap = minval(evalsv(bcouabs%il2,1:nkptnr) - evalsv(bcouabs%iu1,1:nkptnr)&
        &+ input%xs%scissor)
      write(unitout,*)
      write(unitout, '("Info(bse): gap, gap+scissor:", E23.16,1x,E23.16)')&
        & egap-input%xs%scissor, egap
      ! Warn if the system has no gap even with scissor 
      ! (or no scissor and on top of GW)
      if(egap .lt. input%groundstate%epspot) then
        write(unitout,*)
        write(unitout, '("Warning(bse): the system has no gap, setting it to 0")')
        write(unitout,*)
        egap = 0.0d0
      end if  

      !!<-- Inspecting occupancies 
      ! Do not use state combination with zero occupancy difference 
      ! (Theory not complete here?). Also do not use state combinations
      ! where the "unoccupied" state has larger occupancy, since hat will
      ! make the BSE Hamiltonian non-hermitian, even in TDA and q=0.
      !
      ! Set gamma point (currently iq does nothing in checkoccupancies)
      iq = 1
      ! Get adjusted hamsize, index maps and occupation factors.
      ! Also write information about skipped combinations to file.
      call checkoccupancies(iq, hamsize)
      !!-->

      fwp = input%xs%bse%writeparts

      ! Define global distributed Hamiltonian matrix.
      call new_dzmat(dham, hamsize, hamsize, context=ictxt2d)

      ! Write Info
      if(rank == 0) then
        write(unitout,*)
        write(unitout, '("Info(bse): Assembling distributed BSE matrix")')
        write(unitout, '("  RR/RA blocks of global BSE-Hamiltonian:")')
        write(unitout, '("  Shape=",i8)') hamsize
        write(unitout, '("  nk=", i8)') nk
        write(unitout, '("  nkkp=", i8)') nkkp
        write(unitout, '("  Distributing matrix to ",i3," processes")') nproc
        write(unitout, '("  Local matrix schape ",i6," x",i6)')&
          & dham%nrows_loc, dham%ncols_loc
      end if

      ! Assemble Hamiltonian matrix 
      if(rank == 0) call timesec(ts0)
      call setup_distributed_bse(dham)
      if(rank == 0) then
        call timesec(ts1)
        write(unitout, '("All processes build their local matrix")')
        write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0
      end if

      if(fwp) then
        call dzmat_send2global_root(ham, dham)
        if(rank == 0) then
          call writecmplxparts("GlobalHam", dble(ham), immat=aimag(ham))
          deallocate(ham)
        end if
      end if

      ! Number of excitons to consider
      nexc = input%xs%bse%nexc ! (default -1)
      if(nexc > hamsize .or. nexc < 1) then
        nexc = hamsize
      end if

      ! Write Info
      if(rank == 0) then
        write(unitout,*)
        write(unitout, '("Info(bse): Diagonalizing RR Hamiltonian (TDA)")')
        write(unitout, '("Info(bse): Invoking scalapack routine PZHEEVX")')
        write(unitout, '("  Number of requested solutions : ", i8)') nexc
      end if

      ! Eigenvectors are distributed
      call new_dzmat(dbevecr, hamsize, nexc, context=ictxt2d)
      ! Eigenvalues global
      allocate(bevalre(hamsize))

      ! Diagonalize Hamiltonian (destroys the content of ham)
      if(rank == 0) call timesec(ts0)
        call dhesolver(nexc, dham, dbevecr, bevalre, eecs=5)
      if(rank == 0) call timesec(ts1)

      ! Deallocate BSE-Hamiltonian
      call del_dzmat(dham)

      if(rank == 0) then
        write(unitout, '("  Eigen solutions found.")')
        write(unitout, '("  Timing (in seconds)	   :", f12.3)') ts1 - ts0
        write(unitout,*)
      end if

      ! Calculate oscillator strengths.
      call new_dzmat(doszsr, nexc, 3, context=ictxt2d,&
        & rblck=mblck1d_c, cblck=nblck1d_c)
      if(rank == 0) then
        write(unitout, '("Making oszillator strengths (distributed).")')
        call timesec(ts0)
      end if

      call make_doszstren

      if(rank == 0) then 
        call timesec(ts1)
        write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0
      end if

      ! Every process gets a copy of the oscillator strength
      ! (actually only rank 0 writes them to file, but is is not much 
      !  memory and it make the setup for the spectrum calculation easier) 
      call dzmat_send2global_all(oszsr, doszsr)

      if(rank == 0) then
        ! Write excition energies and oscillator strengths to 
        ! text file. 
        write(unitout, '("Writing excition energies and oszillator strengths.")')
        call writeoscillator(hamsize, nexc, egap, bevalre, oszsr)
      end if

      ! Allocate arrays used in spectrum construction
      nw = input%xs%energywindow%points
      allocate(w(nw))
      ! Generate an evenly spaced frequency grid 
      call genwgrid(nw, input%xs%energywindow%intv,&
        & input%xs%tddft%acont, 0.d0, w_real=w)

      if(rank == 0) then
        allocate(symspectr(3,3,nw))
        write(unitout, '("Making spectrum.")')
        if(input%xs%dfoffdiag) then
          write(unitout, '("  Including off-diagonal terms")')
        end if
        call timesec(ts0)
        ! Calculate lattice symmetrized spectrum.
        write(unitout, '("  Using TDA formula.")')
      end if

      call make_dist_spectrum_tda(nw, w)

      if(rank == 0) then
        call timesec(ts1)
        write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0

        write(unitout, '("Writing derived quantities.")')
        ! Generate and write derived optical quantities
        iq = 1
        call writederived(iq, symspectr, nw, w)
        write(unitout, '("  Derived quantities written.")')
      end if

      ! Clean up
      deallocate(bevalre, oszsr, w, evalsv)
      if(associated(input%gw)) deallocate(eval0)
      if(rank == 0) then 
        deallocate(symspectr)
      end if

      call barrier

    else

      write(*,*) "(b_bse): MPI rank", rank, "is idle."

      call barrier

    end if

  else ! not fscal and not rank 0

    call barrier

  end if

contains

  subroutine setuphamiltonian(ham, writeparts)

    !! I/O
    logical, intent(in), optional :: writeparts
    complex(8), allocatable, intent(out) :: ham(:, :)

    !! Local variables
    ! Flags
    logical :: wp
    ! Indices
    integer(4) :: ikkp
    integer(4) :: ik1, ik2
    integer(4) :: s1l, s1u, s2l, s2u
    ! Work arrays
    complex(8), allocatable, dimension(:,:) :: rrham, raham
    complex(8), allocatable, dimension(:,:,:,:) :: excli, sccli
    complex(8), allocatable, dimension(:,:,:,:) :: exclic, scclic
    ! Arrays for optional output, only allocated when writeparts=.true.
    real(8), allocatable :: wre(:,:), wim(:,:), vre(:,:), vim(:,:), kstransen(:,:)
    real(8), allocatable :: wcre(:,:), wcim(:,:), vcre(:,:), vcim(:,:)

    ! Check if matrices should be written out to human readable files.
    if(present(writeparts)) then
      wp = writeparts
    else
      wp = .false.
    end if
    if(wp) then
      allocate(kstransen(hamsize,hamsize))
      kstransen=0.0d0
      ! Real and imaginary parts of resonant-resonant parts
      allocate(wre(hamsize,hamsize),wim(hamsize,hamsize))
      allocate(vre(hamsize,hamsize),vim(hamsize,hamsize))
      wre=0.0d0
      wim=0.0d0
      vre=0.0d0
      vim=0.0d0
      if(fcoup == .true.) then
        ! Real and imaginary parts of resonant-anti-resonant parts
        allocate(wcre(hamsize,hamsize),wcim(hamsize,hamsize))
        allocate(vcre(hamsize,hamsize),vcim(hamsize,hamsize))
        wcre=0.0d0
        wcim=0.0d0
        vcre=0.0d0
        vcim=0.0d0
      end if
    end if

    ! Allocate and zero work arrays
    allocate(rrham(hamsize, hamsize)) ! Resonant-resonant part of Hamiltonian
    allocate(excli(no,nu,no,nu))      ! RR part of V, to be read from file
    allocate(sccli(no,nu,no,nu))      ! RR part of W, to be read from file
    rrham = zzero
    excli = zzero
    sccli = zzero
    if(fcoup == .true.) then
      allocate(raham(hamsize, hamsize)) ! Resonant-anit-resonant part of Hamiltonian
      allocate(exclic(no,nu,no,nu))     ! RA part of V, to be read from file
      allocate(scclic(no,nu,no,nu))     ! RA part of W, to be read from file
      raham = zzero
      exclic = zzero
      scclic = zzero
    end if

    ! Set up kkp blocks of Hamiltonian
    !! Note: If the Hamilton matrix 
    !! has the elements H_{i,j} and the indices enumerate the
    !! states according to
    !! i = {o1u1k1, o1u2k1, ..., o1uMk1,
    !!      o2uMk1, ..., oMuMk1, o1u1k2, ..., oMuMkN} -> {1,...,M**2N}
    !! then because of H_{j,i} = H^*_{i,j} only kj = ki,..,kN is 
    !! needed.
    do ikkp = 1, nkkp

      ! Read corresponding ikkp blocks of W and V from file
      select case(trim(input%xs%bse%bsetype))
        case('singlet', 'triplet')
          ! Read RR part of screened coulomb interaction W_{ouki,o'u'kj}
          call getbsemat('SCCLI.OUT', ikkp, no, nu, sccli)
          if(fcoup == .true.) then
            ! Read RA part of screened coulomb interaction Wc_{ouki,o'u'kj}
            call getbsemat('SCCLIC.OUT', ikkp, no, nu, scclic)
          end if
      end select

      select case(trim(input%xs%bse%bsetype))
        case('RPA', 'singlet')
          ! Read RR part of exchange interaction v_{ouki,o'u'kj}
          call getbsemat('EXCLI.OUT', ikkp, no, nu, excli)
          if(fcoup == .true.) then
            ! Read RA part of exchange interaction vc_{ouki,o'u'kj}
            call getbsemat('EXCLIC.OUT', ikkp, no, nu, exclic)
          end if
      end select

      ! Get ik1 and ik2 from ikkp
      call kkpmap(ikkp, nkptnr, ik1, ik2)

      ! Calculate where to write the ikkp block in the BSE matrix
      s1l = sum(kousize(1:ik1-1)) + 1
      s1u = s1l + kousize(ik1) - 1 
      s2l = sum(kousize(1:ik2-1)) + 1
      s2u = s2l + kousize(ik2) - 1

      !! RR part
      ! For blocks on the diagonal, add the KS transition
      ! energies to the diagonal of the block.
      if(ik1 .eq. ik2) then
        if(wp) then
          call writekstrans(ik1, rrham(s1l:s1u,s2l:s2u), remat=kstransen(s1l:s1u,s2l:s2u))
        else
          call writekstrans(ik1, rrham(s1l:s1u,s2l:s2u))
        end if
      end if
        
      !! RR and RA part
      ! Add exchange term
      ! + 2* v_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
      ! * sqrt(abs(f_{o_{s1} k_{s1}} - f_{u_{s1} k_{s1}}))
      ! * sqrt(abs(f_{o_{s2} k_{s2}} - f_{u_{s2} k_{s2}}))
      select case(trim(input%xs%bse%bsetype))
        case('RPA', 'singlet')
          if(wp) then
            call addexchange(ik1, ik2,&
              & rrham(s1l:s1u,s2l:s2u), excli,&
              & ofac(s1l:s1u), ofac(s2l:s2u),&
              & remat=vre(s1l:s1u,s2l:s2u), immat=vim(s1l:s1u,s2l:s2u))
            if(fcoup == .true.) then
              call addexchange(ik1, ik2,&
                & raham(s1l:s1u,s2l:s2u), exclic,&
                & ofac(s1l:s1u), ofac(s2l:s2u),&
                & remat=vcre(s1l:s1u,s2l:s2u), immat=vcim(s1l:s1u,s2l:s2u))
            end if
          else
            call addexchange(ik1, ik2,&
              & rrham(s1l:s1u,s2l:s2u), excli,&
              & ofac(s1l:s1u), ofac(s2l:s2u))
            if(fcoup == .true.) then
              call addexchange(ik1, ik2,&
                & raham(s1l:s1u,s2l:s2u), exclic,&
                & ofac(s1l:s1u), ofac(s2l:s2u))
            end if
          end if
      end select
      ! Add correlation term
      ! - W_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
      ! * sqrt(abs(f_{o_{s1} k_{s1}} - f_{u_{s1} k_{s1}}))
      ! * sqrt(abs(f_{o_{s2} k_{s2}} - f_{u_{s2} k_{s2}}))
      select case(trim(input%xs%bse%bsetype))
        case('singlet', 'triplet')
          if(wp) then
            call adddirect(ik1, ik2,&
              & rrham(s1l:s1u,s2l:s2u), sccli,&
              & ofac(s1l:s1u), ofac(s2l:s2u),&
              & remat=wre(s1l:s1u,s2l:s2u), immat=wim(s1l:s1u,s2l:s2u))
            if(fcoup == .true.) then
              call adddirect(ik1, ik2,&
                & raham(s1l:s1u,s2l:s2u), scclic,&
                & ofac(s1l:s1u), ofac(s2l:s2u),&
                & remat=wcre(s1l:s1u,s2l:s2u), immat=wcim(s1l:s1u,s2l:s2u))
            end if
          else
            call adddirect(ik1, ik2,&
              & rrham(s1l:s1u,s2l:s2u), sccli, ofac(s1l:s1u), ofac(s2l:s2u))
            if(fcoup == .true.) then
              call adddirect(ik1, ik2,&
                & raham(s1l:s1u,s2l:s2u), scclic,&
                & ofac(s1l:s1u), ofac(s2l:s2u))
            end if
          end if
      end select

    end do

    ! Optional write (use only with small problem sizes!)
    if(wp) then
      call writecmplxparts('E',kstransen)
      call writecmplxparts('V',vre,immat=vim)
      call writecmplxparts('W',wre,immat=wim)
      if(fcoup == .true.) then
        call writecmplxparts('VC',vcre,immat=vcim)
        call writecmplxparts('WC',wcre,immat=wcim)
      end if
      ! Make some room
      deallocate(kstransen)
      deallocate(vre,vim)
      deallocate(wre,wim)
      if(fcoup == .true.) then
        deallocate(vcre,vcim)
        deallocate(wcre,wcim)
      end if
    end if

    !! Assemble full Hamiltonian from blocks
    if(.not. fcoup) then
      ! Just use the RR part, rrham is deallocated afterwards
      call move_alloc(rrham, ham)
    else
      call makefullham(ham, rrham, raham)
    end if


    ! Work arrays
    deallocate(excli)
    deallocate(sccli)
    if(fcoup == .true.) then
      deallocate(rrham)
      deallocate(raham)
      deallocate(exclic)
      deallocate(scclic)
    end if

  end subroutine setuphamiltonian

  subroutine makefullham(ham, rrham, raham)
    complex(8), allocatable, intent(out) :: ham(:,:)
    complex(8), intent(inout) :: rrham(hamsize,hamsize)
    complex(8), intent(inout) :: raham(hamsize,hamsize)

    integer(4) :: i, j

    ! Make rrham explicitly hermitian, since
    ! only the upper triangle was constructed.
    do i=1, hamsize
      ! Set imaginary part of diagonal exactly 0.
      ! (It should be zero anyways, but this is a precaution)
      rrham(i,i) = cmplx(dble(rrham(i,i)), 0.0d0, 8)
      do j=i+1, hamsize
        rrham(j,i) = conjg(rrham(i,j))
      end do
    end do

    ! Make raham explicitly symmetric.
    do i=1, hamsize
      do j=i, hamsize
        raham(j,i) = raham(i,j)
      end do
    end do

    ! Build the 4 blocks of the Hamiltonian
    allocate(ham(2*hamsize,2*hamsize))
    ham = zzero
    ! RR
    ham(1:hamsize, 1:hamsize) = rrham
    ! RA
    ham(1:hamsize, hamsize+1:2*hamsize) = raham
    ! AR
    ham(hamsize+1:2*hamsize, 1:hamsize) = -conjg(raham)
    ! AA
    ham(hamsize+1:2*hamsize, hamsize+1:2*hamsize) = -conjg(rrham)

  end subroutine makefullham

  subroutine writekstrans(ik, hamblock, remat)

    integer(4), intent(in) :: ik
    complex(8), intent(out) :: hamblock(:,:)

    real(8), intent(out), optional :: remat(:,:)

    integer(4) :: io, iu, iou, is
    real(8) :: devals(nou)

    real(8) :: tmp

    ! Calculate ks energy differences
    iou = 0 
    do io = bcouabs%il1, bcouabs%iu1
      do iu = bcouabs%il2, bcouabs%iu2
        iou = iou+1
        ! de = e_{u, k} - e_{o, k} + scissor
        devals(iou) = evalsv(iu,ik) - evalsv(io,ik)&
          &+ input%xs%scissor
      end do
    end do

    is=0
    do iou = 1, nou
      ! Only add if k o u combination was allowed
      if(kouflag(iou, ik) == .true.) then
        is = is +1 
        if(fcoup) then
          ! TEST-ish
          ! Currently not in binding energies, since then the 
          ! spectrum construction would be more tricky. 
          tmp = devals(iou) 
        else
          ! Subtract gap if present, so that the exciton
          ! energies will be binding energies.
          tmp = devals(iou) - egap + bsed
        end if
        if(present(remat)) then
          remat(is,is) = tmp
        end if
        hamblock(is, is) = cmplx(tmp,0.0d0,8)
      end if
    end do

  end subroutine writekstrans

  subroutine addexchange(ik1, ik2, hamblock, excli, ofac1, ofac2, remat, immat)
    integer(4), intent(in) :: ik1, ik2
    real(8), intent(in) :: ofac1(:), ofac2(:)
    complex(8), intent(in) :: excli(no,nu,no,nu)
    complex(8), intent(inout) :: hamblock(:,:)
    
    real(8), intent(out), optional :: remat(:,:), immat(:,:)

    integer(4) :: s1, s2
    integer(4) :: io1, iu1, iou1
    integer(4) :: io2, iu2, iou2
    complex(8) :: tmp

    s1=0
    do iou1 = 1, nou
      if(kouflag(iou1, ik1) == .false.) cycle
      s1 = s1 + 1
      call subhamidx_back(iou1, io1, iu1, bcou%n2)
      s2=0
      do iou2 = 1, nou
        if(kouflag(iou2, ik2) == .false.) cycle
        s2 = s2 + 1
        call subhamidx_back(iou2, io2, iu2, bcou%n2)
        tmp = 2.d0*excli(io1, iu1, io2, iu2) * ofac1(s1) * ofac2(s2)
        if(present(remat)) then
          remat(s1,s2) = dble(tmp)
        end if
        if(present(immat)) then
          immat(s1,s2) = aimag(tmp)
        end if
        hamblock(s1,s2) = hamblock(s1,s2) + tmp
      end do
    end do

  end subroutine addexchange

  subroutine adddirect(ik1, ik2, hamblock, sccli, ofac1, ofac2, remat, immat)
    integer(4), intent(in) :: ik1, ik2
    real(8), intent(in) :: ofac1(:), ofac2(:)
    complex(8), intent(in) :: sccli(no,nu,no,nu)
    complex(8), intent(inout) :: hamblock(:,:)

    real(8), intent(out), optional :: remat(:,:), immat(:,:)
    
    integer(4) :: s1, s2
    integer(4) :: io1, iu1, iou1
    integer(4) :: io2, iu2, iou2
    complex(8) :: tmp


    s1=0
    do iou1 = 1, nou
      if(kouflag(iou1, ik1) == .false.) cycle
      s1 = s1 + 1
      call subhamidx_back(iou1, io1, iu1, bcou%n2)
      s2=0
      do iou2 = 1, nou
        if(kouflag(iou2, ik2) == .false.) cycle
        s2 = s2 + 1
        call subhamidx_back(iou2, io2, iu2, bcou%n2)
        tmp = -sccli(io1, iu1, io2, iu2) * ofac1(s1) * ofac2(s2)

        if(present(remat)) then
          remat(s1,s2) = dble(tmp)
        end if
        if(present(immat)) then
          immat(s1,s2) = aimag(tmp)
        end if
        hamblock(s1,s2) = hamblock(s1,s2) + tmp

      end do
    end do

  end subroutine adddirect

  subroutine make_doszstren
    use mod_constants, only: zone
    use m_getpmat
    
    implicit none

    ! Local
    real(8) :: t1, t0
    type(dzmat) :: drmat

    ! Distributed position operator matrix
    call new_dzmat(drmat, hamsize, 3, context=ictxt2d,&
      & rblck=mblck1d_c, cblck=nblck1d_c)

    ! Build position operator matrix elements
    if(rank == 0) then 
      write(unitout, '("  Building Rmat.")')
      call timesec(t0)
    end if
    call setup_distributed_rmat(drmat)
    if(rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Building resonant oscillator strengths.")')
      call timesec(t0)
    end if

    !! Resonant oscillator strengths
    ! t^R_\lambda,i = < \tilde{R}^i | X_\lambda> =
    ! \Sum_{s1} f_{s1} R^*_{{u_{s1} o_{s1} k_{s1}},i} X_{o_{s1} u_{s1} k_{s1}},\lambda= 
    ! \Sum_{s1} f_{s1} P^*_{{u_{s1} o_{s1} k_{s1}},i} / 
    !   (e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) X_{o_{s1} u_{s1} k_{s1}},\lambda= 
    ! \Sum_{s1} f_{s1} P_{{o_{s1} u_{s1} k_{s1}},i} / 
    !   (e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) X_{o_{s1} u_{s1} k_{s1}},\lambda
    call dzgemm(dbevecr, drmat, doszsr, transa='T')
    if(rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if

    call del_dzmat(drmat)
  end subroutine make_doszstren

  
  subroutine makeoszillatorstrength(oszstrr, oszstra)
    use mod_constants, only: zone
    use m_getpmat
    
    implicit none

    ! I/O
    complex(8), intent(out) :: oszstrr(nexc,3)
    complex(8), intent(out), optional :: oszstra(nexc, 3)

    ! Local
    real(8) :: t1, t0
    complex(8), allocatable :: pmou(:,:,:), rmat(:,:)
    integer(4) :: io, iu, ioabs, iuabs, ik, ikprev
    integer(4) :: a1
    
    ! Allocate momentum matrix slice needed.
    allocate(pmou(3, no, nu))
    allocate(rmat(hamsize, 3))

    ! Loop over 3 directions
    ikprev = -1

    write(unitout, '("  Building Rmat.")')
    call timesec(t0)
    do a1 = 1, hamsize
      
      ! Get indices
      ik = smap(a1,1)
      io = smap(a1,2)
      iu = smap(a1,3)               ! Relative 
      ioabs = io + bcouabs%il1 - 1  ! Absolute
      iuabs = iu + bcouabs%il2 - 1  ! Absolute

      ! Read momentum matrix slice for given k-point
      ! if not already present (k index varies the slowest in the smap)
      if(ik /= ikprev) then 
        call getpmat(ik, vkl,& 
          & bcouabs%il1, bcouabs%iu1,&
          & bcouabs%il2, bcouabs%iu2,&
          & .true., 'PMAT_XS.OUT', pmou)
      end if
      ikprev = ik

      ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
      ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
      !   Where E are quasi-particle energies and e are KS energies.
      if(associated(input%gw)) then
         pmou(1:3,io,iu) = pmou(1:3,io,iu)&
           &* (evalsv(iuabs,ik)-evalsv(ioabs,ik))/(eval0(iuabs,ik)- eval0(ioabs,ik))
      end if 

      ! Build complex conjugate R-matrix from p-matrix
      ! \tilde{R}^*_{u_{s1},o_{s1},k_{s1}},i = 
      ! (f_{o_{s1},k_{s1}}-f_{u_{s1},k_{s1}}) *
      !   P_{o_{s1},u_{s1},k_{s1}},i /(e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}})
      rmat(a1, :) = ofac(a1) * pmou(:, io, iu)/(evalsv(ioabs, ik) - evalsv(iuabs, ik))

    end do
    ! Momentum matrix elements no longer needed
    deallocate(pmou)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Building resonant oscillator strengths.")')
    call timesec(t0)
    oszstrr = zzero
    !! Resonant oscillator strengths
    ! t^R_\lambda,i = < \tilde{R}^i | X_\lambda> =
    ! \Sum_{s1} f_{s1} R^*_{{u_{s1} o_{s1} k_{s1}},i} X_{o_{s1} u_{s1} k_{s1}},\lambda= 
    ! \Sum_{s1} f_{s1} P^*_{{u_{s1} o_{s1} k_{s1}},i} / 
    !   (e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) X_{o_{s1} u_{s1} k_{s1}},\lambda= 
    ! \Sum_{s1} f_{s1} P_{{o_{s1} u_{s1} k_{s1}},i} / 
    !   (e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) X_{o_{s1} u_{s1} k_{s1}},\lambda
    call zgemm('t','n', nexc, 3, hamsize,&
      & zone, bevecr(1:hamsize,1:nexc), hamsize, rmat, hamsize, zzero, oszstrr, nexc)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    if(fcoup .and. present(oszstra)) then
      write(unitout, '("  Building anti-resonant oscillator strengths.")')
      call timesec(t0)
      oszstra = zzero
      !! Anti-resonant oscillator strengths
      ! t^A_\lambda,i = < \tilde{R}^{i*} | Y_\lambda> =
      ! \Sum_{s1} f_{s1} R_{{u_{s1} o_{s1} k_{s1}},i} Y_{o_{s1} u_{s1} k_{s1}},\lambda= 
      ! \Sum_{s1} f_{s1} P_{{u_{s1} o_{s1} k_{s1}},i} / 
      !   (e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) Y_{o_{s1} u_{s1} k_{s1}},\lambda= 
      ! \Sum_{s1} f_{s1} P^*_{{o_{s1} u_{s1} k_{s1}},i} / i
      !   (e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) Y_{o_{s1} u_{s1} k_{s1}},\lambda
      call zgemm('t','n', nexc, 3, hamsize,&
        & zone, bevecr(hamsize+1:2*hamsize,1:nexc), hamsize,&
        & conjg(rmat), hamsize, zzero, oszstra, nexc)
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if

    deallocate(rmat)

  end subroutine makeoszillatorstrength

  subroutine makespectrum(nfreq, freq, spectrum)
    use mod_lattice, only: omega
    use mod_constants, only: zone, zi, pi
    use modxs, only: symt2
    use invert
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(3,3,nfreq)

    ! Local
    integer(4) :: i, j, o1, o2
    real(8) :: t1,t0
    complex(8) :: buf(3,3,nfreq)
    complex(8), allocatable :: ns_spectr(:,:)
    complex(8), allocatable :: overlap(:,:), invoverlap(:,:)
    complex(8), allocatable :: amat(:,:), bmat(:,:), enw(:,:), o(:,:), op(:,:)

    ! Calculate overlap matrix between right eigenvectors
    write(unitout, '("  Making overlap matrix between right EVs.")')
    call timesec(t0)
    allocate(overlap(nexc,nexc))
    call zgemm('C','N', nexc, nexc, 2*hamsize, zone,&
      & bevecr, 2*hamsize, bevecr, 2*hamsize, zzero, overlap, nexc)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    ! Invert overlap matrix
    write(unitout, '("  Inverting said overlap matrix.")')
    call timesec(t0)
    allocate(invoverlap(nexc,nexc))
    call zinvert_lapack(overlap, invoverlap)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    ! Test output of overlap matrices
    if(fwp) then
      call writecmplxparts('Q',dble(overlap),immat=aimag(overlap))
      call writecmplxparts('invQ',dble(invoverlap),immat=aimag(invoverlap))
    end if
    ! Acctual overlap no longer needed
    deallocate(overlap)
    
    ! Make combined resonan-anti-resonant oscillator strengths
    allocate(o(nexc,3), op(nexc,3))
    o = oszsr+oszsa
    op = conjg(oszsr-oszsa)

    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)
    ! Make energy denominator for each frequency 
    allocate(enw(nfreq, nexc))
    do i = 1, nexc
      do j = 1, nfreq
        enw(j,i) = zone/(bevalre(i)-freq(j)-zi*input%xs%broad)
      end do
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    
    write(unitout, '("  Making helper matrix amat.")')
    ! Make helper matrices 
    ! A_\lambda,i = \Sum_\lambda' Q^-1_\lambda,\lambda' \tilde{O}^*_\lambda',i
    call timesec(t0)
    allocate(amat(nexc,3))
    call zgemm('N','N', nexc, 3, nexc, zone, invoverlap, nexc,&
      & op, nexc, zzero, amat, nexc)
    ! Inverse of overlap no longer need
    deallocate(invoverlap)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Making helper matrix bmat.")')
    ! B_\lambda,ij = A_\lambda,i*O_\lambda,j
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      allocate(bmat(nexc,9))
      do o1 = 1, 3
        do o2 = 1, 3
          j = o2 + (o1-1)*3
          bmat(:,j) = amat(:,o1)*o(:,o2)
        end do
      end do
    else
      allocate(bmat(nexc,3))
      do o1 = 1, 3
        bmat(:,o1) = amat(:,o1)*o(:,o1)
      end do
    end if
    ! amat not needed anymore
    deallocate(amat)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
         
    write(unitout, '("  Calculating spectrum.")')
    ! Make non-lattice-symmetrized spectrum
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      allocate(ns_spectr(nfreq,9))
      call zgemm('N','N', nfreq, 9, nexc, zone, enw, nfreq,&
        & bmat, nexc, zzero, ns_spectr, nfreq)
      ns_spectr = ns_spectr*8.d0*pi/omega/nkptnr
      ns_spectr(:,1) = zone + ns_spectr(:,1)
      ns_spectr(:,4) = zone + ns_spectr(:,4)
      ns_spectr(:,7) = zone + ns_spectr(:,7)
      ! Write to buffer for symmetry routine
      buf=zzero
      do o1=1,3
        do o2=1,3
          j = o2 + (o1-1)*3
          buf(o1,o2,:) = ns_spectr(:,j)
        end do
      end do
    else
      allocate(ns_spectr(nfreq,3))
      call zgemm('N','N', nfreq, 3, nexc, zone, enw, nfreq,&
        & bmat, nexc, zzero, ns_spectr, nfreq)
      ns_spectr = zone + ns_spectr*8.d0*pi/omega/nkptnr
      ! Write to buffer for symmetry routine
      buf=zzero
      do o1=1,3
        buf(o1,o1,:) = ns_spectr(:,o1)
      end do
    end if
    ! Helper no longer needed
    deallocate(bmat,enw,ns_spectr)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Symmetrizing spectrum.")')
    ! Symmetrize spectrum 
    call timesec(t0)
    spectrum = zzero
    do o1=1,3
      do o2=1,3
        ! Symmetrize the macroscopic dielectric tensor
        call symt2app(o1, o2, nfreq, symt2, buf, spectrum(o1,o2,:))
      end do 
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
  end subroutine makespectrum

  subroutine make_dist_spectrum_tda(nfreq, freq)
    use mod_lattice, only: omega
    use mod_constants, only: zone, zi, pi
    use modxs, only: symt2
    use invert
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)

    ! Local
    integer(4) :: i, j, o1, o2, ig, jg, nopt
    real(8) :: t1, t0
    complex(8), allocatable :: buf(:,:,:)
    complex(8), allocatable :: ns_spectr(:,:)
    type(dzmat) :: denwr, denwa, dbmatr, dns_spectr

    if(rank == 0) then
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)
    end if

    ! Shift energies back to absolute values.
    bevalre = bevalre + egap - bsed

    ! Make energy denominator for each frequency (resonant & anti-resonant) 
    call new_dzmat(denwr, nfreq, nexc, context=ictxt2d)
    do i = 1, denwr%nrows_loc
      do j = 1, denwr%ncols_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = indxl2g(i, denwr%mblck, myprow, 0, nprow)
        jg = indxl2g(j, denwr%nblck, mypcol, 0, npcol)
        denwr%za(i,j) = zone/(bevalre(jg)-freq(ig)-zi*input%xs%broad)
#else
        denwr%za(i,j) = zone/(bevalre(j)-freq(i)-zi*input%xs%broad)
#endif
      end do
    end do
    call new_dzmat(denwa, nfreq, nexc, context=ictxt2d)
    do i = 1, denwa%nrows_loc
      do j = 1, denwa%ncols_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = indxl2g(i, denwa%mblck, myprow, 0, nprow)
        jg = indxl2g(j, denwa%nblck, mypcol, 0, npcol)
        denwa%za(i,j) = zone/(bevalre(jg)+freq(ig)+zi*input%xs%broad)
#else
        denwa%za(i,j) = zone/(bevalre(j)+freq(i)+zi*input%xs%broad)
#endif
      end do
    end do

    if(rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Making helper matrix bmat.")')
      call timesec(t0)
    end if
    
    if(input%xs%dfoffdiag) then
      nopt = 9
    else
      nopt = 3
    end if

    call new_dzmat(dbmatr, nexc, nopt, context=ictxt2d,&
      & rblck=mblck1d_c, cblck=nblck1d_c)
    do i = 1, dbmatr%nrows_loc
      do j = 1, dbmatr%ncols_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = indxl2g(i, dbmatr%mblck, myprow, 0, nprow)
        jg = indxl2g(j, dbmatr%nblck, mypcol, 0, npcol)
#else
        ig = i
        jg = j
#endif
        ! Get individual opical indices
        if(input%xs%dfoffdiag) then
          o1 = (jg-1)/3 + 1
          o2 = jg-(o1-1)*3
        else
          o1 = jg
          o2 = jg
        end if
        dbmatr%za(i,j) = oszsr(ig,o1)*conjg(oszsr(ig,o2))
      end do
    end do

    if(rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)
    end if
         
    ! Make non-lattice-symmetrized spectrum
    call new_dzmat(dns_spectr, nfreq, nopt, context=ictxt2d,&
      & rblck=mblck1d_c, cblck=nblck1d_c)

    call dzgemm(denwr, dbmatr, dns_spectr)

    dbmatr%za = conjg(dbmatr%za)

    call dzgemm(denwa, dbmatr, dns_spectr, beta=zone)

    ! Add 1 to diagonal elements o1=o2
    if(input%xs%dfoffdiag) then
      dns_spectr%za = dns_spectr%za*8.d0*pi/omega/nkptnr
      do i = 1, dns_spectr%nrows_loc
        do j = 1, dns_spectr%ncols_loc
#ifdef SCAL
          ig = indxl2g(i, dns_spectr%mblck, myprow, 0, nprow)
          jg = indxl2g(j, dns_spectr%nblck, mypcol, 0, npcol)
#else
          ig = i
          jg = j
#endif
          if(jg == 1 .or. jg == 4 .or. jg == 7) then
            dns_spectr%za(i,j) = zone + dns_spectr%za(i,j)
          end if
        end do
      end do
    else
      dns_spectr%za = zone + dns_spectr%za*8.d0*pi/omega/nkptnr
    end if
      
    call del_dzmat(denwr)
    call del_dzmat(denwa)
    call del_dzmat(dbmatr)

    ! Send spectrum to a global matrix at rank 0 to 
    ! interface with non parallel post-processing routines.
    call dzmat_send2global_root(ns_spectr, dns_spectr)
    call del_dzmat(dns_spectr)
    if(rank == 0) then 

      ! Write to buffer for symmetry routine
      if(allocated(buf)) deallocate(buf)
      allocate(buf(3,3,nfreq))
      buf=zzero

      do o1=1,3
        do o2=1,3
          j = o2 + (o1-1)*3
          buf(o1,o2,:) = ns_spectr(:,j)
        end do
      end do
      deallocate(ns_spectr)

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0

      write(unitout, '("  Symmetrizing spectrum.")')
      call timesec(t0)

      symspectr = zzero
      do o1=1,3
        do o2=1,3
          ! Symmetrize the macroscopic dielectric tensor
          call symt2app(o1, o2, nfreq, symt2, buf, symspectr(o1,o2,:))
        end do 
      end do
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    end if

  end subroutine make_dist_spectrum_tda

  subroutine makespectrum_tda(nfreq, freq, spectrum)
    use mod_lattice, only: omega
    use mod_constants, only: zone, zi, pi
    use modxs, only: symt2
    use invert
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(3,3,nfreq)

    ! Local
    integer(4) :: i, j, o1, o2
    real(8) :: t1, t0
    complex(8) :: buf(3,3,nfreq)
    complex(8), allocatable :: ns_spectr(:,:)
    complex(8), allocatable :: enwr(:,:), enwa(:,:), bmatr(:,:)

    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)
    ! Shift energies back to absolute values.
    bevalre = bevalre + egap - bsed
    ! Make energy denominator for each frequency (resonant & anti-resonant) 
    allocate(enwr(nfreq, nexc))
    do i = 1, nexc
      do j = 1, nfreq
        enwr(j,i) = zone/(bevalre(i)-freq(j)-zi*input%xs%broad)
      end do
    end do
    allocate(enwa(nfreq, nexc))
    do i = 1, nexc
      do j = 1, nfreq
        enwa(j,i) = zone/(bevalre(i)+freq(j)+zi*input%xs%broad)
      end do
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    
    write(unitout, '("  Making helper matrix bmat.")')
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      allocate(bmatr(nexc,9))
      do o1 = 1, 3
        do o2 = 1, 3
          j = o2 + (o1-1)*3
          bmatr(:,j) = oszsr(:,o1)*conjg(oszsr(:,o2))
        end do
      end do
    else
      allocate(bmatr(nexc,3))
      do o1 = 1, 3
        bmatr(:,o1) = oszsr(:,o1)*conjg(oszsr(:,o1))
      end do
    end if
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
         
    write(unitout, '("  Calculating spectrum.")')
    ! Make non-lattice-symmetrized spectrum
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      allocate(ns_spectr(nfreq,9))
      call zgemm('N','N', nfreq, 9, nexc, zone, enwr, nfreq,&
        & bmatr, nexc, zzero, ns_spectr, nfreq)
      call zgemm('N','N', nfreq, 9, nexc, zone, enwa, nfreq,&
        & conjg(bmatr), nexc, zone, ns_spectr, nfreq)
      ns_spectr = ns_spectr*8.d0*pi/omega/nkptnr
      ns_spectr(:,1) = zone + ns_spectr(:,1)
      ns_spectr(:,4) = zone + ns_spectr(:,4)
      ns_spectr(:,7) = zone + ns_spectr(:,7)
      ! Write to buffer for symmetry routine
      buf=zzero
      do o1=1,3
        do o2=1,3
          j = o2 + (o1-1)*3
          buf(o1,o2,:) = ns_spectr(:,j)
        end do
      end do
    else
      allocate(ns_spectr(nfreq,3))
      call zgemm('N','N', nfreq, 3, nexc, zone, enwr, nfreq,&
        & bmatr, nexc, zzero, ns_spectr, nfreq)
      call zgemm('N','N', nfreq, 3, nexc, zone, enwa, nfreq,&
        & conjg(bmatr), nexc, zone, ns_spectr, nfreq)
      ns_spectr = zone + ns_spectr*8.d0*pi/omega/nkptnr
      ! Write to buffer for symmetry routine
      buf=zzero
      do o1=1,3
        buf(o1,o1,:) = ns_spectr(:,o1)
      end do
    end if
    ! Helper no longer needed
    deallocate(bmatr, enwr, enwa, ns_spectr)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Symmetrizing spectrum.")')
    ! Symmetrize spectrum 
    call timesec(t0)
    spectrum = zzero
    do o1=1,3
      do o2=1,3
        ! Symmetrize the macroscopic dielectric tensor
        call symt2app(o1, o2, nfreq, symt2, buf, spectrum(o1,o2,:))
      end do 
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
  end subroutine makespectrum_tda

end subroutine b_bse
!EOC
