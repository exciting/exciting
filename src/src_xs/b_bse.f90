! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_bse
! !INTERFACE:
subroutine b_bse
! !USES:
  ! Basics
  use modinput, only: input
  use mod_constants, only: zzero, h2ev
  use mod_kpoint, only: nkptnr, vkl
  use mod_eigenvalue_occupancy, only: evalsv, nstsv
  ! MPI and BLACS/ScaLAPACK
  use modmpi
  use modscl
  ! XS
  use modxs, only: unitout, bcbs, bsed
  ! BSE
  use modbse
  use m_putgetbsemat
  use m_genwgrid
  use m_genfilname
  use m_diagfull
  use m_writeoscillator
  use m_writecmplxparts
  use m_dhesolver
  use m_hesolver
  use m_dzgemm
  use m_setup_bse
  use m_setup_rmat
  use m_storeexcitons
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
  ! Variables
  integer(4) :: iknr, iq, iqmt
  integer(4) :: nexc
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

    write(*,*) "b_bse: Running non parallel version."

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

    ! Set EVALSV_QMTXXX.OUT as basis for the occupation limits search
    !   Note: EVALSV_QMT000.OUT is always produced,
    !         EVALSV_QMT001.OUT has the same content, when
    !         the first entry in the q-point list is set 0 0 0
    ! To be exact the following genfilname set filext to _QMTXXX.OUT
    iqmt = 0
    call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

    ! Set ist* variables and ksgap in modxs using findocclims
    ! This also reads in 
    ! (QMTXXX)
    ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
    ! (QMT000)
    ! modxs:evalsv0, modxs:occsv0
    if(read_eval_occ_qmt) then
      call setranges_modxs(iqmt)
    end if

    ! Select relevant transitions for the construction
    ! of the BSE hamiltonian
    ! Also sets nkkp_bse, nk_bse 
    if(seltrans) then
      call select_transitions(iqmt)
    end if

    ! If on top of GW
    if(associated(input%gw)) then
      ! Save KS eigenvalues to use them later for renormalizing PMAT
      allocate(eval0(nstsv, nkptnr))
      eval0=evalsv
      ! If scissor correction is presented, one should nullify it
      input%xs%scissor=0.0d0
      ! Read QP Fermi energies and eigenvalues from file
      call getevalqp(nkptnr,vkl,evalsv)
      write(unitout,'("  Quasi particle energies are read from EVALQP.OUT")')
    end if

    ! Energy to shift the BSE eigenvalues by:
    evalshift = -egap+bsed+input%xs%scissor

    ! Allocate frequency array used in spectrum construction
    allocate(w(nw))
    ! Generate an evenly spaced frequency grid 
    call genwgrid(nw, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)

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

    write(unitout,*)
    write(unitout, '("Info(bse): gap, gap+scissor:", E23.16,1x,E23.16)')&
      & egap-input%xs%scissor, egap
    write(unitout, '("Info(bse): Shifting evals by:", E23.16)')&
      evalshift
    ! Warn if the system has no gap even with scissor (or no scissor and on top of GW)
    if(egap .lt. input%groundstate%epspot) then
      write(unitout,*)
      write(unitout, '("Warning(bse): the system has no gap, setting it to 0")')
      write(unitout,*)
      egap = 0.0d0
    end if  

    fwp = input%xs%bse%writeparts

    ! Write Info
    if(rank == 0) then
      write(unitout,*)
      write(unitout, '("Info(bse): Assembling BSE matrix")')
      write(unitout, '("  RR/RA blocks of global BSE-Hamiltonian:")')
      write(unitout, '("  Shape=",i8)') hamsize
      write(unitout, '("  nk_bse=", i8)') nk_bse
      write(unitout, '("  nkkp_bse=", i8)') nkkp_bse
      if(fcoup) then
        write(unitout, '(" Including coupling terms ")')
        write(unitout, '(" Full BSE-Hamiltonian:")')
        write(unitout, '("  Shape=",i8)') 2*hamsize
      end if
      if(fwp) then
        write(unitout, '(" Writing real and imaginary parts of Hamiltonian to file ")')
      end if
    end if
    ! Assemble Hamiltonian matrix 
    call timesec(ts0)
    if(fcoup) then 
      call setup_full_hamiltonian(ham)
    else
      call setup_bse(ham, iqmt, .false.)
    end if
    call timesec(ts1)
    write(unitout, '(" Matrix build.")')
    write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0

    ! Write Info
    if(rank == 0) then
      write(unitout,*)
      if(fcoup) then
        write(unitout, '("Info(bse): Diagonalizing full non symmetric Hamiltonian")')
        write(unitout, '("Info(bse): Invoking lapack routine ZGEEVX")')
      else
        write(unitout, '("Info(bse): Diagonalizing RR Hamiltonian (TDA)")')
        write(unitout, '("Info(bse): Invoking lapack routine ZHEEVR")')
      end if
    end if

    ! Allocate eigenvector and eigenvalue arrays
    if(fcoup) then 
      allocate(bevalre(2*hamsize), bevecr(2*hamsize, 2*hamsize))
      allocate(bevalim(2*hamsize))
    else
      allocate(bevalre(hamsize), bevecr(hamsize, hamsize))
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
      ! Set nexc to half the number of excitions
      nexc = hamsize
    else
      call hesolver(ham, bevecr, bevalre,&
       & v1=wl-ewidth-egap, v2=wu+ewidth-egap, found=nexc)
    end if
    call timesec(ts1)

    ! Test write out right-eigenvectors
    if(fwp) then
      call writecmplxparts('bevecr', dble(bevecr), immat=aimag(bevecr))
    end if

    ! Deallocate BSE-Hamiltonian
    deallocate(ham)

    if(fcoup) then 
      write(unitout, '("  All eigen solutions found.")')
    else
      write(unitout, '("  ",i8," eigen solutions found.")') nexc
    end if
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
    write(unitout,&
      & '("Writing excition energies and oszillator strengths to text file.")')
    if(fcoup) then
      call writeoscillator(2*hamsize, nexc, egap, bevalre, oszsr,&
        & evalim=bevalim, oszstra=oszsa)
    else
      call writeoscillator(hamsize, nexc, egap, bevalre, oszsr)
    end if

    ! Allocate arrays used in spectrum construction
    allocate(symspectr(3,3,nw))

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
    call writederived(iqmt+1, symspectr, nw, w)
    write(unitout, '("  Derived quantities written.")')

    ! Store excitonic energies and wave functions to file
    if(associated(input%xs%storeexcitons)) then
      call storeexcitons(nexc,bevalre,bevecr,fcoup)
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

    ! Set up process grids for BLACS 
    !   Make square'ish process grid
    call setupblacs(mpiglobal, '2d', bi2d)
    !   Also make 1d grids with the same number of processes
    call setupblacs(mpiglobal, '1dc', bi1dc, np=bi2d%nprocs)
    !call setupblacs(mpiglobal, '1dr', bi1dr, np=bi2d%nprocs)

    ! Only MPI ranks that have an associated 
    ! BLACS grid process do anything.
    if( .not. bi1dc%myprow < 0 ) then 

      ! General init
      call init0
      ! k-grid init
      call init1
      ! q-grid init
      call init2
      ! Saves all k grid related variables in modxs, so 
      ! that another k grid can be used for k+q.
      call xssave0

      ! Read Fermi energy from file
      call readfermi

      ! Set EVALSV_QMTXXX.OUT as basis for the occupation limits search
      !   Note: EVALSV_QMT000.OUT is always produced,
      !         EVALSV_QMT001.OUT has the same content, when
      !         the first entry in the q-point list is set 0 0 0
      ! To be exact the following genfilname set filext to _QMTXXX.OUT
      iqmt = 0
      call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

      ! Set ist* variables and ksgap in modxs using findocclims
      ! This also reads in 
      ! (QMTXXX)
      ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
      ! (QMT000)
      ! modxs:evalsv0, modxs:occsv0
      if(read_eval_occ_qmt) then
        call setranges_modxs(iqmt)
      end if

      ! Select relevant transitions for the construction
      ! of the BSE hamiltonian
      ! Also sets nkkp_bse, nk_bse 
      if(seltrans) then
        call select_transitions(iqmt)
      end if

      ! If on top of GW
      if(associated(input%gw)) then
        ! Save KS eigenvalues to use them later for renormalizing PMAT
        allocate(eval0(nstsv,nkptnr))
        eval0(:,:)=evalsv(:,:)
        ! If scissor correction is presented, one should nullify it
        input%xs%scissor=0.0d0
        ! Read QP Fermi energies and eigenvalues from file
        call getevalqp(nkptnr,vkl,evalsv)
        write(unitout,'("  Quasi particle energies are read from EVALQP.OUT")')
      end if

      ! Energy to shift the BSE eigenvalues by:
      evalshift = -egap+bsed+input%xs%scissor

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

      write(unitout,*)
      write(unitout, '("Info(bse): gap, gap+scissor:", E23.16,1x,E23.16)')&
        & egap-input%xs%scissor, egap
      write(unitout, '("Info(bse): Shifting evals by:", E23.16)')&
        evalshift
      ! Warn if the system has no gap even with scissor 
      ! (or no scissor and on top of GW)
      if(egap .lt. input%groundstate%epspot) then
        write(unitout,*)
        write(unitout, '("Warning(bse): the system has no gap, setting it to 0")')
        write(unitout,*)
        egap = 0.0d0
      end if  

      fwp = input%xs%bse%writeparts

      ! Define global distributed Hamiltonian matrix.
      call new_dzmat(dham, hamsize, hamsize, bi2d)

      ! Write Info
      if(mpiglobal%rank == 0) then
        write(unitout,*)
        write(unitout, '("Info (bse): Assembling distributed BSE matrix")')
        write(unitout, '("  RR/RA blocks of global BSE-Hamiltonian:")')
        write(unitout, '("  Shape=",i8)') hamsize
        write(unitout, '("  nk=", i8)') nk_bse
        write(unitout, '("  nkkp=", i8)') nkkp_bse
        write(unitout, '("  Distributing matrix to ",i3," processes")') bi2d%nprocs
        write(unitout, '("  Local matrix schape ",i6," x",i6)')&
          & dham%nrows_loc, dham%ncols_loc
      end if

      ! Assemble RR/RA block of Hamiltonian matrix
      if(mpiglobal%rank == 0) call timesec(ts0)
      call setup_distributed_bse(dham, iqmt, fcoup, bi2d)
      if(mpiglobal%rank == 0) then
        call timesec(ts1)
        write(unitout, '("All processes build their local matrix")')
        write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0
      end if
      
      if(fwp) then
        call dzmat_send2global_root(ham, dham, bi2d)
        if(mpiglobal%rank == 0) then
          call writecmplxparts("GlobalHam", dble(ham), immat=aimag(ham))
          deallocate(ham)
        end if
      end if

      ! Write Info
      if(rank == 0) then
        write(unitout,*)
        write(unitout, '("Info(bse): Diagonalizing RR Hamiltonian (TDA)")')
        write(unitout, '("Info(bse): Invoking scalapack routine PZHEEVX")')
        !write(unitout, '("  Number of requested solutions : ", i8)') nexc
      end if

      ! Eigenvectors are distributed
      ! must be NxN because Scalapack solver expects it
      call new_dzmat(dbevecr, hamsize, hamsize, bi2d) 
      ! Eigenvalues global
      allocate(bevalre(hamsize))

      ! Diagonalize Hamiltonian (destroys the content of ham)
      if(mpiglobal%rank == 0) call timesec(ts0)
        call dhesolver(dham, dbevecr, bevalre, bi2d,&
         & v1=wl-ewidth-egap, v2=wu+ewidth-egap, found=nexc,&
         & eecs=input%xs%bse%eecs)
      if(mpiglobal%rank == 0) call timesec(ts1)

      ! Deallocate BSE-Hamiltonian
      call del_dzmat(dham)

      if(mpiglobal%rank == 0) then
        write(unitout, '("  Eigen solutions found:", i8)') nexc
        write(unitout, '("  Timing (in seconds)	   :", f12.3)') ts1 - ts0
        write(unitout,*)
      end if

      ! Calculate oscillator strengths.
      call new_dzmat(doszsr, nexc, 3, bi2d,&
        & rblck=1, cblck=bi2d%nblck)
      if(mpiglobal%rank == 0) then
        write(unitout, '("Making oszillator strengths (distributed).")')
        call timesec(ts0)
      end if

      call make_doszstren

      if(mpiglobal%rank == 0) then 
        call timesec(ts1)
        write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0
      end if

      ! Eigen vectors no longer needed.
      call del_dzmat(dbevecr)

      ! Every process gets a copy of the oscillator strength
      ! (actually only rank 0 writes them to file, but is is not much 
      !  memory and it make the setup for the spectrum calculation easier) 
      call dzmat_send2global_all(oszsr, doszsr, bi2d)

      if(mpiglobal%rank == 0) then
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

      if(mpiglobal%rank == 0) then
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

      if(mpiglobal%rank == 0) then
        call timesec(ts1)
        write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0

        write(unitout, '("Writing derived quantities.")')
        ! Generate and write derived optical quantities
        call writederived(iqmt+1, symspectr, nw, w)
        write(unitout, '("  Derived quantities written.")')
      end if

      ! Clean up
      deallocate(bevalre, oszsr, w, evalsv)
      if(associated(input%gw)) deallocate(eval0)
      if(mpiglobal%rank == 0) then 
        deallocate(symspectr)
      end if

      call barrier

    else

      write(*,*) "(b_bse): MPI rank", mpiglobal%rank, "is idle."

      call barrier

    end if

  else ! not fscal and not rank 0

    call barrier

  end if

contains

  ! NOTE: only for qmt = 0
  subroutine setup_full_hamiltonian(ham)

    ! I/O
    complex(8), intent(inout) :: ham(:, :)

    integer(4) :: i, j

    ! RR
    call setup_bse(ham(1:hamsize,1:hamsize), iqmt, .false.)
    ! Make RR part of ham explicitly hermitian, since
    ! only the upper triangle was constructed.
    do i=1, hamsize
      ! Set imaginary part of diagonal exactly 0.
      ! (It should be zero anyways, but this is a precaution)
      ham(i,i) = cmplx(dble(ham(i,i)), 0.0d0, 8)
      do j=i+1, hamsize
        ham(j,i) = conjg(ham(i,j))
      end do
    end do
    ! RA
    call setup_bse(ham(1:hamsize,hamsize+1:hamsize*2), iqmt, .true.)
    ! Make RA part of ham explicitly symmetric.
    do i=1, hamsize
      do j=hamsize+i, 2*hamsize
        ham(j-hamsize,i+hamsize) = ham(i,j)
      end do
    end do
    ! AR
    ham(hamsize+1:2*hamsize, 1:hamsize) = -conjg(ham(1:hamsize,hamsize+1:hamsize*2))
    ! AA
    ham(hamsize+1:2*hamsize, hamsize+1:2*hamsize) = -conjg(ham(1:hamsize,1:hamsize))

  end subroutine setup_full_hamiltonian

  subroutine make_doszstren
    use mod_constants, only: zone
    use m_getpmat
    
    implicit none

    ! Local
    real(8) :: t1, t0
    type(dzmat) :: drmat

    ! Distributed position operator matrix
    call new_dzmat(drmat, hamsize, 3, bi2d,&
      & rblck=1, cblck=bi2d%nblck)

    ! Build position operator matrix elements
    if(mpiglobal%rank == 0) then 
      write(unitout, '("  Building Rmat.")')
      call timesec(t0)
    end if
    ! Build R-matrix from P-matrix
    ! \tilde{R}_{a,i} = 
    !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
    !     P^i_{o_a,u_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a})
    call setup_distributed_rmat(drmat, bi2d)
    if(mpiglobal%rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Building resonant oscillator strengths.")')
      call timesec(t0)
    end if

    !! Resonant oscillator strengths
    ! t^R_{\lambda,i} = < \tilde{R}^{i*} | X_\lambda> =
    !   ( \Sum_{a} \tilde{R}^T_{i, a} X_{a, \lambda} )^T =
    !     \Sum_{a} X^T_{\lambda, a} \tilde{R}_{a,i}
    call dzgemm(dbevecr, drmat, doszsr, transa='T', m=nexc)
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
    complex(8), allocatable :: pmouk(:,:,:,:), rmat(:,:)
    integer(4) :: io, iu, ioabs, iuabs, ik, iknr
    integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 
    integer(4) :: a1
    
    ! Allocate momentum matrix slice needed.
    allocate(pmouk(3, no_bse_max, nu_bse_max, nk_bse))
    allocate(rmat(hamsize, 3))

    write(unitout, '("  Reading Pmat.")')
    call timesec(t0)
    if(associated(input%gw)) then
      write(unitout, '("  Also renormalizing Pmat with GW eigenvalues.")')
    end if
    ! Read in all needed momentum matrix elements
    do ik = 1, nk_bse
      iknr = kmap_bse_rg(ik)
      iuabs1 = koulims(1,iknr)
      iuabs2 = koulims(2,iknr)
      ioabs1 = koulims(3,iknr)
      ioabs2 = koulims(4,iknr)
      inu = iuabs2 - iuabs1 + 1
      ino = ioabs2 - ioabs1 + 1
      call getpmat(iknr, vkl,&
        & ioabs1, ioabs2, iuabs1, iuabs2,&
        & .true., 'PMAT_XS.OUT', pmouk(:,1:ino,1:inu,ik))
      ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
      ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
      !   Where E are quasi-particle energies and e are KS energies.
      if(associated(input%gw)) then
        do io = 1, ino
          do iu = 1, inu
            pmouk(:,io,iu,ik) = pmouk(:,io,iu,ik)&
              &* (evalsv(iuabs,iknr)-evalsv(ioabs,iknr))/(eval0(iuabs,iknr)- eval0(ioabs,iknr))
          end do
        end do
      end if 
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Building Rmat.")')
    call timesec(t0)
    do a1 = 1, hamsize

      ! Absolute indices
      iuabs = smap(a1,1)
      ioabs = smap(a1,2)
      iknr  = smap(a1,3)

      ! Relative indices
      iu = smap_rel(a1,1)
      io = smap_rel(a1,2)
      ik = smap_rel(a1,3)
      
      ! Build R-matrix from P-matrix
      ! \tilde{R}_{a,i} = 
      !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
      !     P^i_{o_a,u_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a})
      rmat(a1, :) = ofac(a1) * pmouk(:, io, iu, ik)/(evalsv(iuabs, iknr) - evalsv(ioabs, iknr))

    end do
    ! Momentum matrix elements no longer needed
    deallocate(pmouk)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Building resonant oscillator strengths.")')
    call timesec(t0)
    !! Resonant oscillator strengths
    ! t^R_{\lambda,i} = < \tilde{R}^{i*} | X_\lambda> =
    !   ( \Sum_{a} \tilde{R}^T_{i, a} X_{a, \lambda} )^T =
    !     \Sum_{a} X^T_{\lambda, a} \tilde{R}_{a,i}
    call zgemm('t','n', nexc, 3, hamsize,&
      & zone, bevecr(1:hamsize,1:nexc), hamsize, rmat, hamsize, zzero, oszstrr, nexc)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    if(fcoup .and. present(oszstra)) then
      write(unitout, '("  Building anti-resonant oscillator strengths.")')
      call timesec(t0)
      !! Anti-resonant oscillator strengths
      ! t^A_{\lambda,i} = < \tilde{R}^{i} | Y_\lambda> =
      !   ( \Sum_{a} \tilde{R}^\dag_{i, a} Y_{a, \lambda} )^T =
      !     \Sum_{a} Y^T_{\lambda, a} \tilde{R}^*_{a,i}
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
    complex(8), allocatable :: tplus(:,:), tmat(:,:), enw(:,:), tminus(:,:), op(:,:)

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
    allocate(tminus(nexc,3), op(nexc,3))
    tminus = oszsr-oszsa
    op = conjg(oszsr+oszsa)

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
    
    write(unitout, '("  Making helper matrix tplus.")')
    ! Make helper matrices 
    ! t^+_\lambda,j = \Sum_\lambda' Q^-1_{\lambda,\lambda'} \tilde{O}_\lambda',j
    ! = \Sum_\lambda' Q^-1_{\lambda,\lambda'} ( t^R_{\lambda', j} + t^A_{\lambda',j})^*
    call timesec(t0)
    allocate(tplus(nexc,3))
    call zgemm('N','N', nexc, 3, nexc, zone, invoverlap, nexc,&
      & op, nexc, zzero, tplus, nexc)
    ! Inverse of overlap no longer need
    deallocate(invoverlap)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Making helper matrix tmat.")')
    ! tmat_{\lambda,ij} = O^*_{\lambda,i} t^+_{\lambda,j}
    !              = t^-_{\lambda, i} t^+_{\lambda,j}
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      allocate(tmat(nexc,9))
      do o1 = 1, 3
        do o2 = 1, 3
          j = o2 + (o1-1)*3
          tmat(:,j) = tminus(:,o1)*tplus(:,o2)
        end do
      end do
    else
      allocate(tmat(nexc,3))
      do o1 = 1, 3
        tmat(:,o1) = tminus(:,o1)*tplus(:,o1)
      end do
    end if
    ! tplus and tminus not needed anymore
    deallocate(tplus, tminus)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
         
    write(unitout, '("  Calculating spectrum.")')
    ! Make non-lattice-symmetrized spectrum
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      allocate(ns_spectr(nfreq,9))
      ! \epsilon^M_{w,ij} \prop \Sum_\lambda E^-1_{w,\lambda} tmat_{\lambda,ij}
      call zgemm('N','N', nfreq, 9, nexc, zone, enw, nfreq,&
        & tmat, nexc, zzero, ns_spectr, nfreq)
      ns_spectr = ns_spectr*8.d0*pi/omega/nkptnr
      ns_spectr(:,1) = zone + ns_spectr(:,1)
      ns_spectr(:,5) = zone + ns_spectr(:,5)
      ns_spectr(:,9) = zone + ns_spectr(:,9)
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
        & tmat, nexc, zzero, ns_spectr, nfreq)
      ns_spectr = zone + ns_spectr*8.d0*pi/omega/nkptnr
      ! Write to buffer for symmetry routine
      buf=zzero
      do o1=1,3
        buf(o1,o1,:) = ns_spectr(:,o1)
      end do
    end if
    ! Helper no longer needed
    deallocate(tmat,enw,ns_spectr)
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
    type(dzmat) :: denwr, denwa, dtmatr, dns_spectr

    if(mpiglobal%rank == 0) then
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)
    end if

    ! Shift energies back to absolute values.
    bevalre = bevalre + egap - bsed

    ! Make energy denominator for each frequency (resonant & anti-resonant) 
    ! enwr_{w, \lambda} = 1/(E_\lambda - w - i\delta)
    call new_dzmat(denwr, nfreq, nexc, bi2d)
    do j = 1, denwr%ncols_loc
#ifdef SCAL
      jg = denwr%c2g(j)
#endif
      do i = 1, denwr%nrows_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = denwr%r2g(i)
        denwr%za(i,j) = zone/(bevalre(jg)-freq(ig)-zi*input%xs%broad)
#else
        denwr%za(i,j) = zone/(bevalre(j)-freq(i)-zi*input%xs%broad)
#endif
      end do
    end do
    ! enwa_{w, \lambda} = 1/(E_\lambda + w + i\delta)
    call new_dzmat(denwa, nfreq, nexc, bi2d)
    do j = 1, denwa%ncols_loc
#ifdef SCAL
      jg = denwa%c2g(j)
#endif
      do i = 1, denwa%nrows_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = denwa%r2g(i)
        denwa%za(i,j) = zone/(bevalre(jg)+freq(ig)+zi*input%xs%broad)
#else
        denwa%za(i,j) = zone/(bevalre(j)+freq(i)+zi*input%xs%broad)
#endif
      end do
    end do

    if(mpiglobal%rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Making helper tmatr bmat.")')
      call timesec(t0)
    end if
    
    if(input%xs%dfoffdiag) then
      nopt = 9
    else
      nopt = 3
    end if

    ! tmatr_{\lambda, j} = t^R_{\lambda, o1_j} t^{R*}_{\lambda, o2_j} 
    ! where j combines the 2 cartesian directions
    call new_dzmat(dtmatr, nexc, nopt, bi2d,&
      & rblck=1, cblck=bi2d%nblck)
    do j = 1, dtmatr%ncols_loc
#ifdef SCAL
        jg = dtmatr%c2g(j)
#endif
      do i = 1, dtmatr%nrows_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = dtmatr%r2g(i)
#else
        ig = i
        jg = j
#endif
        ! Get individual opical indices
        if(input%xs%dfoffdiag) then
          o2 = (jg-1)/3 + 1
          o1 = jg-(o2-1)*3
        else
          o2 = jg
          o1 = jg
        end if
        dtmatr%za(i,j) = oszsr(ig,o1)*conjg(oszsr(ig,o2))
      end do
    end do

    if(mpiglobal%rank == 0) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)
    end if
         
    ! Make non-lattice-symmetrized spectrum
    call new_dzmat(dns_spectr, nfreq, nopt, bi2d,&
      & rblck=1, cblck=bi2d%nblck)

    ! nsspectr_{w,j} = \Sum_{\lambda} enwr_{w,\lambda} tmatr_{\lambda, j}
    !   i.e. nsspectr_{w,j} = 
    !     \Sum_{\lambda} 1/(E_\lambda - w - i\delta)
    !       t^R_{\lambda, o1_j} t^{R*}_{\lambda, o2_j} 
    call dzgemm(denwr, dtmatr, dns_spectr)
    ! nsspectr_{w,j} += \Sum_{\lambda} enwa_{w,\lambda} tmatr^*_{\lambda, j}
    !   i.e. nsspectr_{w,j} = nsspectr_{w,j}
    !     \Sum_{\lambda} 1/(E_\lambda + w + i\delta)
    !       t^{R*}_{\lambda, o1_j} t^R_{\lambda, o2_j} 
    dtmatr%za = conjg(dtmatr%za)
    call dzgemm(denwa, dtmatr, dns_spectr, beta=zone)

    ! Add 1 to diagonal elements o1=o2
    if(input%xs%dfoffdiag) then
      dns_spectr%za = dns_spectr%za*8.d0*pi/omega/nkptnr
      do j = 1, dns_spectr%ncols_loc
#ifdef SCAL
        jg = dns_spectr%c2g(j)
#endif
        do i = 1, dns_spectr%nrows_loc
#ifdef SCAL
          ig = dns_spectr%r2g(i)
#else
          ig = i
          jg = j
#endif
          if(jg == 1 .or. jg == 5 .or. jg == 9) then
            dns_spectr%za(i,j) = zone + dns_spectr%za(i,j)
          end if
        end do
      end do
    else
      dns_spectr%za = zone + dns_spectr%za*8.d0*pi/omega/nkptnr
    end if
      
    call del_dzmat(denwr)
    call del_dzmat(denwa)
    call del_dzmat(dtmatr)

    ! Send spectrum to a global matrix at rank 0 to 
    ! interface with non parallel post-processing routines.
    call dzmat_send2global_root(ns_spectr, dns_spectr, bi2d)
    call del_dzmat(dns_spectr)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Symmetrizing spectrum.")')
      call timesec(t0)

      ! Write to buffer for symmetry routine
      if(allocated(buf)) deallocate(buf)
      allocate(buf(3,3,nfreq))
      buf=zzero

      if(input%xs%dfoffdiag) then
        do o2=1,3
          do o1=1,3
            j = o1 + (o2-1)*3
            buf(o1,o2,:) = ns_spectr(:,j)
          end do
        end do
      else
        do o1=1,3
          buf(o1,o1,:) = ns_spectr(:,o1)
        end do
      end if
      deallocate(ns_spectr)

      symspectr = zzero
      do o2=1,3
        do o1=1,3
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
    complex(8), allocatable :: enwr(:,:), enwa(:,:), tmatr(:,:)

    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)
    ! Shift energies back to absolute values.
    bevalre = bevalre + egap - bsed

    ! Make energy denominator for each frequency (resonant & anti-resonant) 
    ! enwr_{w, \lambda} = 1/(E_\lambda - w - i\delta)
    allocate(enwr(nfreq, nexc))
    do j = 1, nexc
      do i = 1, nfreq
        enwr(i,j) = zone/(bevalre(j)-freq(i)-zi*input%xs%broad)
      end do
    end do
    ! enwa_{w, \lambda} = 1/(E_\lambda + w + i\delta)
    allocate(enwa(nfreq, nexc))
    do j = 1, nexc
      do i = 1, nfreq
        enwa(i,j) = zone/(bevalre(j)+freq(i)+zi*input%xs%broad)
      end do
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    
    write(unitout, '("  Making helper matrix tmatr.")')
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      ! tmatr_{\lambda, j} = t^R_{\lambda, o1_j} t^{R*}_{\lambda, o2_j} 
      ! where j combines the 2 cartesian directions
      allocate(tmatr(nexc,9))
      do o2 = 1, 3
        do o1 = 1, 3
          j = o1 + (o2-1)*3
          tmatr(:,j) = oszsr(:,o1)*conjg(oszsr(:,o2))
        end do
      end do
    else
      ! tmatr_{\lambda, o1} = t^R_{\lambda, o1} t^{R*}_{\lambda, o1} 
      allocate(tmatr(nexc,3))
      do o1 = 1, 3
        tmatr(:,o1) = oszsr(:,o1)*conjg(oszsr(:,o1))
      end do
    end if
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
         
    write(unitout, '("  Calculating spectrum.")')
    if(input%xs%dfoffdiag) then
      write(unitout, '("    Including off diagonals.")')
    end if
    ! Make non-lattice-symmetrized spectrum
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      ! nsspectr_{w,j} = \Sum_{\lambda} enwr_{w,\lambda} tmatr_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} 1/(E_\lambda - w - i\delta)
      !       t^R_{\lambda, o1_j} t^{R*}_{\lambda, o2_j} 
      allocate(ns_spectr(nfreq,9))
      call zgemm('N','N', nfreq, 9, nexc, zone, enwr, nfreq,&
        & tmatr, nexc, zzero, ns_spectr, nfreq)
      ! nsspectr_{w,j} += \Sum_{\lambda} enwa_{w,\lambda} tmatr^*_{\lambda, j}
      !   i.e. nsspectr_{w,j} = nsspectr_{w,j} +
      !     \Sum_{\lambda} 1/(E_\lambda + w + i\delta)
      !       t^{R*}_{\lambda, o1_j} t^R_{\lambda, o2_j} 
      call zgemm('N','N', nfreq, 9, nexc, zone, enwa, nfreq,&
        & conjg(tmatr), nexc, zone, ns_spectr, nfreq)
      ! Adjusting prfactor and add 1 to diagonal elements
      ns_spectr = ns_spectr*8.d0*pi/omega/nkptnr
      ns_spectr(:,1) = zone + ns_spectr(:,1)
      ns_spectr(:,5) = zone + ns_spectr(:,5)
      ns_spectr(:,9) = zone + ns_spectr(:,9)
      ! Write to buffer for symmetry routine
      buf=zzero
      do o2=1,3
        do o1=1,3
          j = o1 + (o2-1)*3
          buf(o1,o2,:) = ns_spectr(:,j)
        end do
      end do
    else
      ! Same as in other case, but j=o1=o2
      allocate(ns_spectr(nfreq,3))
      call zgemm('N','N', nfreq, 3, nexc, zone, enwr, nfreq,&
        & tmatr, nexc, zzero, ns_spectr, nfreq)
      call zgemm('N','N', nfreq, 3, nexc, zone, enwa, nfreq,&
        & conjg(tmatr), nexc, zone, ns_spectr, nfreq)
      ns_spectr = zone + ns_spectr*8.d0*pi/omega/nkptnr
      ! Write to buffer for symmetry routine
      buf=zzero
      do o1=1,3
        buf(o1,o1,:) = ns_spectr(:,o1)
      end do
    end if
    ! Helper no longer needed
    deallocate(tmatr, enwr, enwa, ns_spectr)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Symmetrizing spectrum.")')
    ! Symmetrize spectrum 
    call timesec(t0)
    spectrum = zzero
    do o2=1,3
      do o1=1,3
        ! Symmetrize the macroscopic dielectric tensor
        call symt2app(o1, o2, nfreq, symt2, buf, spectrum(o1,o2,:))
      end do 
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
  end subroutine makespectrum_tda

end subroutine b_bse
!EOC
