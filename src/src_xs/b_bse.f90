! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_bse
! !INTERFACE:
subroutine b_bse(iqmt)
! !USES:
  ! Basics
  use mod_misc, only: filext
  use modinput, only: input
  use mod_constants, only: zone, zi, zzero, pi, h2ev
  use mod_kpoint, only: nkptnr, vkl
  use mod_eigenvalue_occupancy, only: evalsv, nstsv
  ! MPI and BLACS/ScaLAPACK
  use modmpi
  use modscl
  ! XS
  use modxs, only: evalsv0, unitout, bcbs, ksgapval,&
                 & qgap, ikmapikq, iqmtgamma, vqlmt,&
                 & vkl0, ivgmt
  ! BSE
  use modbse
  ! Spectrum
  use invert
  use mod_lattice, only: omega
  use modxs, only: symt2
  use mod_symmetry, only: nsymcrys
  ! Interface modules
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
  use m_setup_pwmat
  use m_genexevec
  use m_putgetexcitons
  use modgw, only: eferqp
  use mod_eigenvalue_occupancy, only: efermi

use m_writecmplxparts
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

  ! I/O
  integer(4), intent(in) :: iqmt

  ! Local variables
  ! Variables
  integer(4) :: iknr, iq, ik, ikq, io, iu, a1
  integer(4) :: nexc
  real(8) :: ts0, ts1
  logical :: fcoup, fwp, fscal, fti

  ! Allocatable arrays
  real(8), allocatable, dimension(:) :: bevalim, bevalre, w
  complex(8), allocatable, dimension(:,:) :: ham, bevecr, bevecaux
  complex(8), allocatable, dimension(:,:) :: resvec, aresvec 
  complex(8), allocatable, dimension(:,:) :: oszsr, oszsa
  complex(8), allocatable, dimension(:,:) :: cmat, cpmat
  complex(8), allocatable, dimension(:,:,:) :: symspectr

  real(8) :: bsegap
  real(8) :: v1, v2, en1, en2
  integer(4) :: i1, i2, i,j, iex1, iex2, nreq
  integer(4) :: nsymcrys_save
  logical :: efind
  ! Distributed arrays
  integer(4) :: ip
  complex(8), allocatable :: buff(:,:)
  type(dzmat) :: dham, dbevecr, doszsr, dresvec, daresvec
  type(dzmat) :: dcmat, dcpmat

  ! Write Hamilton matrix parts to readable file
  fwp = input%xs%bse%writeparts

  ! Include coupling terms
  fcoup = input%xs%bse%coupling

  ! Use time inverted anti-resonant basis
  fti = input%xs%bse%ti

  ! Use ScaLAPACK
  fscal = input%xs%bse%distribute

  if(fscal .and. fcoup .and. .not. fti ) then 
    write(*,*) "Coupling and scalapack not supported in standard basis"
    call terminate
  end if
  if(iqmt/=1 .and. .not. fti) then 
    write(*,*) "Finite momentum transfer only supported in ti ar basis"
    call terminate
  end if

  ! Non-parallelized code.
  if(.not. fscal .and. mpiglobal%rank == 0) then 

    !write(*,*) "b_bse: Running non parallel version."

    ! General init
    call init0
    ! k-grid init
    call init1
    ! Save variables of the unshifted (apart from xs:vkloff) k grid 
    ! to modxs (vkl0, ngk0, ...)
    call xssave0
    ! q-point and qmt-point setup
    !   Init 2 sets up (task 445):
    !   * A list of momentum transfer vectors form the q-point list 
    !     (modxs::vqlmt and mod_qpoint::vql)
    !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
    !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
    !   * G+qmt quantities (modxs)
    !   * The square root of the Coulomb potential for the G+qmt points
    !   * Reads STATE.OUT
    !   * Generates radial functions (mod_APW_LO)
    call init2

    !write(*,*) "b_bse: iqmt=", iqmt
    !write(*,*) "b_bse: vqlmt=", vqlmt(1:3, iqmt)

    ! Read Fermi energy from file 
    ! (only needed in setranges_modxs::findocclims to check whether system has a gap)
    ! Use EFERMI_QMT001.OUT
    call genfilname(iqmt=iqmtgamma, setfilext=.true.)
    call readfermi

    ! Set ist* variables and ksgap in modxs using findocclims
    ! This also reads in 
    ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
    ! modxs:evalsv0, modxs:occsv0
    call setranges_modxs(iqmt, fcoup, fti)

    ! WARNING: ONTOP OF GW STILL IS INCONSISTENT, SINCE OCCUPATION SEARCH IS
    ! NOT DONE ON EVALQP.OUT !?!
    ! If on top of GW
    if(associated(input%gw) .and. iqmt==1) then
      ! Save KS eigenvalues to use them later for renormalizing PMAT
      allocate(eval0(nstsv, nkptnr))
      eval0=evalsv
      ! If scissor correction is presented, one should nullify it
      input%xs%scissor=0.0d0
      ! Read QP Fermi energies and eigenvalues from file
      ! NOTE: QP evals are shifted by -efermi-eferqp with respect to KS evals
      ! NOTE: getevalqp sets mod_symmetry::nsymcrys to 1
      ! NOTE: getevalqp needs the KS eigenvalues as input
      nsymcrys_save = nsymcrys
      call getevalqp(nkptnr,vkl0,evalsv)
      nsymcrys = nsymcrys_save
      !if(mpiglobal%rank == 0) then 
      !  write(*,'("efermi=", f10.7)') efermi
      !  write(*,'("eferqp=", f10.7)') eferqp
      !  do ik=1, nkptnr
      !    write(*,'("k-point #", i6,":", 3f10.7)') ik, vkl0(:,ik)
      !    write(*,'(" state    Eks    E+efermi+eferqp    E+eferqp    E+efermi    E")')
      !    do i=input%gw%ibgw, input%gw%nbgw
      !      write(*,'(i6, 5(4x, f10.7))') i, eval0(i,ik), evalsv(i, ik)+efermi+eferqp,&
      !        & evalsv(i, ik)+eferqp, evalsv(i, ik)+efermi, evalsv(i, ik) 
      !    end do
      !  end do
      !  write(unitout,'("Info(b_bse): Quasi particle energies are read from EVALQP.OUT")')
      !end if
      ! Set k and k'=k grid eigenvalues to QP energies
      evalsv0=evalsv
    else if(associated(input%gw) .and. iqmt /= 1) then 
      if(mpiglobal%rank==0) then 
        write(*,'("Error(b_bse): BSE+GW only supported for 0 momentum transfer.")')
      end if
      call terminate
    end if

    ! Select relevant transitions for the construction
    ! of the BSE hamiltonian
    ! Also sets nkkp_bse, nk_bse 
    ! NOTE: Resets vkl0 to k grid and vkl to k+qmt grid
    ! Resets eigenvalues and occupancies correspondingly
    ! NOTE: If GW and iqmt=1, then it reads in QP energies.
    ! NOTE: Sets transition energy array used on the diagonal of the BSE (modbse::de)
    call select_transitions(iqmt, serial=.true.)

    ! Determine "BSE"-gap, i.e. the lowest energy 
    ! KS transition which is considered.
    bsegap = 1.0d16
    !$OMP PARALLEL DO &
    !$OMP& DEFAULT(SHARED), PRIVATE(a1,iu,io,ik,ikq),&
    !$OMP& REDUCTION(min:bsegap)
    do a1 = 1, hamsize
      iu = smap(1,a1)
      io = smap(2,a1)
      ik = smap(3,a1)
      ikq = ikmapikq(ik, iqmt)
      bsegap = min(bsegap, evalsv(iu,ik)-evalsv0(io,ikq))
    end do
    !$OMP END PARALLEL DO

    write(unitout,*)
    write(unitout, '("Info(b_bse):&
      & bsegap, bsegap+scissor (eV):", E23.16,1x,E23.16)')&
      & bsegap*h2ev, (bsegap+sci)*h2ev

    ! Write Info
    write(unitout,*)
    write(unitout, '("Info(b_bse): Assembling BSE matrix")')
    write(unitout, '("  RR/RA blocks of global BSE-Hamiltonian:")')
    write(unitout, '("  Shape=",i8)') hamsize
    write(unitout, '("  nk_bse=", i8)') nk_bse
    if(fcoup) then
      write(unitout, '(" Including coupling terms ")')
      if(fti) then 
        write(unitout, '(" Using time inverted anti-resonant basis")')
        write(unitout, '(" Using squared EVP")')
      else
        write(unitout, '(" Full BSE-Hamiltonian:")')
        write(unitout, '("  Shape=",i8)') 2*hamsize
      end if
    end if
    if(fwp) then
      write(unitout, '("Info(b_bse):&
        & Writing real and imaginary parts of Hamiltonian to file ")')
    end if

    ! Assemble Hamiltonian matrix 
    if(fcoup) then 
      if(fti) then 
        call setup_ti_bse(ham, cmat, cpmat,  iqmt)
        if(fwp) then 
          call writecmplxparts('HamS', dble(ham), immat=aimag(ham))
          call writecmplxparts('cmat', dble(cmat), immat=aimag(cmat))
          call writecmplxparts('cpmat', dble(cpmat), immat=aimag(cpmat))
        end if
      else
        allocate(ham(2*hamsize,2*hamsize))
        call setup_full_bse(ham, iqmt)
        if(fwp) then 
          call writecmplxparts('Global_Ham', dble(ham), immat=aimag(ham))
        end if
      end if
    else
      allocate(ham(hamsize,hamsize))
      call setup_bse(ham, iqmt, .false., .false.)
    end if

    ! Write Info
    write(unitout,*)
    if(fcoup) then
      if(fti) then 
        write(unitout, '("Info(b_bse): Solving hermitian squared EVP")')
        write(unitout, '("Info(b_bse): Invoking lapack routine ZHEEVR")')
      else
        write(unitout, '("Info(b_bse): Diagonalizing full non symmetric Hamiltonian")')
        write(unitout, '("Info(b_bse): Invoking lapack routine ZGEEVX")')
      end if
    else
      write(unitout, '("Info(b_bse): Diagonalizing RR Hamiltonian (TDA)")')
      write(unitout, '("Info(b_bse): Invoking lapack routine ZHEEVR")')
    end if

    ! Allocate eigenvector and eigenvalue arrays
    if(fcoup) then 
      if(fti) then 
        allocate(bevalre(hamsize))
        allocate(bevecaux(hamsize, hamsize))
      else
        allocate(bevalre(2*hamsize))
        allocate(bevalim(2*hamsize))
        allocate(bevecr(2*hamsize, 2*hamsize))
      end if
    else
      allocate(bevalre(hamsize), bevecr(hamsize, hamsize))
    end if
    bevalre = 0.0d0

    ! Find only eigenvalues relevant for requested 
    ! spectrum?
    efind = input%xs%bse%efind

    ! Diagonalize Hamiltonian (destroys the content of ham)
    call timesec(ts0)
    if(fcoup) then
      if(.not. fti) then 
        call diagfull(2*hamsize, ham, bevalre,&
          & evalim=bevalim, evecr=bevecr, fbalance=.false., frcond=.false., fsort=.true.)
        nexc = 2*hamsize
      else
        if(efind) then
          v1=(max(wl, 0.0d0))**2
          v2=(wu)**2
          call hesolver(hemat=ham, evec=bevecaux, eval=bevalre,&
           & v1=v1, v2=v2, found=nexc)
        else
          if(input%xs%bse%nexc == -1) then 
            i2 = hamsize
          else
            i2 = input%xs%bse%nexc
          end if
          i1 = 1
          call hesolver(hemat=ham, evec=bevecaux, eval=bevalre,&
           & i1=i1, i2=i2, found=nexc)
        end if
      end if
    else
      if(efind) then
          v1=max(wl, 0.0d0)
          v2=wu
        call hesolver(hemat=ham, evec=bevecr, eval=bevalre,&
         & v1=v1, v2=v2, found=nexc)
      else
        if(input%xs%bse%nexc == -1) then 
          i2 = hamsize
        else
          i2 = input%xs%bse%nexc
        end if
        i1 = 1
        call hesolver(hemat=ham, evec=bevecr, eval=bevalre,&
         & i1=i1, i2=i2, found=nexc)
      end if
    end if
    call timesec(ts1)

    !write(*,*) "writing s evals"
    !call writecmplxparts('nd_s_evals', revec=bevalre, veclen=size(bevalre))

    ! Test write out right-eigenvectors
    if(fwp) then
      if(fcoup .and. fti) then 
        call writecmplxparts('bevecaux', dble(bevecaux), immat=aimag(bevecaux))
      else
        call writecmplxparts('bevecr', dble(bevecr), immat=aimag(bevecr))
        call writecmplxparts('evals',revec=dble(bevalre),veclen=size(bevalre))
      end if
    end if

    ! Deallocate BSE-Hamiltonian
    deallocate(ham)

    if(fcoup .and. fti) then 

      ! Take square root of auxiliary eigenvalues
      ! to retrieve actual eigenvalues
      if(any(bevalre < 0.0d0)) then 
        write(*,*) "Error(b_bse): Negative squared EVP evals occured"
        write(*,'(E23.16)') bevalre
        call terminate
      end if
      bevalre = sqrt(bevalre)
      
    end if

    if(fcoup) then 
      if(fti) then 
        if(efind) then
          write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
          write(unitout, '("  [",E10.3,",",E10.3,"]/H^2")') v1, v2
          write(unitout, '("  [",E10.3,",",E10.3,"]/eV^2")') v1*h2ev**2, v2*h2ev**2
        else
          write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
          write(unitout, '("  i1=",i8," i2=",i8)') i1, i2
        end if
      else
        write(unitout, '("  All eigen solutions found.")')
      end if
    else
      if(efind) then
        write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
        write(unitout, '("  [",E10.3,",",E10.3,"]/H")') v1, v2
        write(unitout, '("  [",E10.3,",",E10.3,"]/eV")') v1*h2ev, v2*h2ev
      else
        write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
        write(unitout, '("  i1=",i8," i2=",i8)') i1, i2
      end if
    end if
    write(unitout, '("  Timing (in seconds)	   :", f12.3)') ts1 - ts0
    write(unitout,*)

    ! Store excitonic energies and wave functions to file
    if(associated(input%xs%storeexcitons)) then

      if(input%xs%storeexcitons%selectenergy) then 
        en1=input%xs%storeexcitons%minenergyexcitons
        en2=input%xs%storeexcitons%maxenergyexcitons
        !write(*,*) "b_bse: en1, en2", en1, en2
        if(input%xs%storeexcitons%useev) then 
          en1=en1/h2ev
          en2=en2/h2ev
        end if
        !write(*,*) "b_bse: en1, en2", en1, en2
        if(fcoup .and. .not. fti) then 
          call energy2index(2*hamsize, nexc, bevalre, en1, en2, iex1, iex2)
        else
          call energy2index(hamsize, nexc, bevalre, en1, en2, iex1, iex2)
        endif 
        !write(*,*) "b_bse: iex1, iex2", iex1, iex2
      else
        iex1=input%xs%storeexcitons%minnumberexcitons
        iex2=input%xs%storeexcitons%maxnumberexcitons
      end if
      nreq=iex2-iex1+1
      !write(*,*) "b_bse: iex1=", iex1, " iex2=", iex2, " nreq=", nreq

      if(nreq < 1 .or. nreq > nexc .or. iex1<1 .or. iex2<1 ) then
        write(*,*) "Error(b_bse): storeexcitons index mismatch."
        write(*,*) "iex1, iex2, nreq, nex", iex1, iex2, nreq, nexc
        call terminate
      end if

      write(unitout, '("Info(b_bse): Writing excition eigenvectors for index range=",2i8)') iex1, iex2

      if(fcoup .and. fti) then 
        allocate(resvec(hamsize, nreq))
        allocate(aresvec(hamsize, nreq))

        write(unitout, '("Info(b_bse): Generating resonant&
          & and anti-resonant exciton coefficients from&
          & auxilliary squared EVP eigenvectors (pos. E).")')
        call genexevec(iex1, iex2, nexc, cmat, cpmat, bevecaux, bevalre,&
          & rvecp=resvec, avecp=aresvec)

        write(unitout, '("Info(b_bse): Writing exciton eigenvectors to file (pos. E).")')
        call put_excitons(bevalre(iex1:iex2), rvec=resvec, avec=aresvec,&
          & iqmt=iqmt, a1=iex1, a2=iex2)

        if(allocated(resvec)) deallocate(resvec)
        if(allocated(aresvec)) deallocate(aresvec)

      else if(fcoup .and. .not. fti) then 

        write(unitout, '("Info(b_bse): Writing exciton eigenvectors to file.")')
        call put_excitons(bevalre(iex1:iex2), bevecr(1:hamsize,iex1:iex2),&
          & avec=bevecr(hamsize+1:2*hamsize,iex1:iex2),&
          & iqmt=iqmt, a1=iex1, a2=iex2)

      else

        write(unitout, '("Info(b_bse): Writing exciton eigenvectors to file.")')
        call put_excitons(bevalre(iex1:iex2), bevecr(:,iex1:iex2),&
          & iqmt=iqmt, a1=iex1, a2=iex2)

      end if

    end if


    ! Calculate oscillator strengths.
    allocate(oszsr(nexc,3))
    if(fcoup .and. .not. fti) then
      allocate(oszsa(nexc,3))
    end if

    if(fcoup .and. .not. fti) then 
      call makeoszillatorstrength(oszsr, oszstra=oszsa)
    else
      call makeoszillatorstrength(oszsr)
    end if

    ! Write excition energies and oscillator strengths to 
    ! text file. 
    write(unitout, '("Info(b_bse):&
      & Writing excition energies and oszillator strengths to text file.")')
    if(fcoup .and. .not. fti) then
      call writeoscillator(2*hamsize, nexc, bsegap+sci, bevalre, oszsr,&
        & evalim=bevalim, oszstra=oszsa, sort=.true.)
    else
      call writeoscillator(hamsize, nexc, -(bsegap+sci), bevalre, oszsr, iqmt=iqmt)
    end if

    ! Allocate frequency array used in spectrum construction
    allocate(w(nw))
    ! Generate an evenly spaced frequency grid 
    call genwgrid(nw, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)

    ! Allocate arrays used in spectrum construction
    allocate(symspectr(3,3,nw))

    ! Calculate lattice symmetrized spectrum.
    if(fcoup) then 
      if(fti) then
        call makespectrum_ti(nw, w, symspectr)
      else
        !call makespectrum_full(nw, w, symspectr)
        call makespectrum_full_lr(nw, w, symspectr)
      end if
    else
      call makespectrum_tda(nw, w, symspectr)
    end if

    write(unitout, '("Info(b_bse): Writing derived quantities.")')
    ! Generate and write derived optical quantities
    call writederived(iqmt, symspectr, nw, w)
    write(unitout, '("  Derived quantities written.")')

    ! Clean up
    deallocate(bevalre, oszsr, w, symspectr, evalsv)
    if(allocated(bevecr)) deallocate(bevecr)
    if(fcoup) then 
      if(fti) then 
        deallocate(bevecaux, cmat, cpmat)
      else
        deallocate(bevalim, oszsa)
      end if
    end if
    if(associated(input%gw)) deallocate(eval0)

    write(unitout, '("BSE calculation finished")')

    ! Rank 0 says I am finished to all others
    call barrier

  ! Parallel version 
  else if (fscal) then
    
    !write(*,*) "b_bse: Running parallel version at rank:", mpiglobal%rank

    ! Set up process grids for BLACS 
    !   Make square'ish process grid (context 0)
    call setupblacs(mpiglobal, 'grid', bi2d)
    !   Also make 1d grid with the same number of processes (context 1)
    call setupblacs(mpiglobal, 'row', bi1d, np=bi2d%nprocs)
    !   Also make 0d grid containing only the current processes (context 2)
    call setupblacs(mpiglobal, '0d', bi0d, np=1)

    ! General init
    call init0
    ! k-grid init
    call init1
    ! Save variables of the unshifted (apart from xs:vkloff) k grid 
    ! to modxs (vkl0, ngk0, ...)
    call xssave0
    ! q-point and qmt-point setup
    !   Init 2 sets up (task 445):
    !   * A list of momentum transfer vectors form the q-point list 
    !     (modxs::vqmtl and mod_qpoint::vql)
    !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
    !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
    !   * G+qmt quantities (modxs)
    !   * The square root of the Coulomb potential for the G+qmt points
    !   * Reads STATE.OUT
    !   * Generates radial functions (mod_APW_LO)
    call init2

    ! Read Fermi energy from file
    ! Use EFERMI_QMT001.OUT
    call genfilname(iqmt=iqmtgamma, setfilext=.true.)
    call readfermi

    ! Set ist* variables and ksgap in modxs using findocclims
    ! This also reads in 
    ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
    ! modxs:evalsv0, modxs:occsv0
    call setranges_modxs(iqmt, fcoup, fti)

    ! WARNING: ONTOP OF GW STILL IS INCONSISTENT, SINCE OCCUPATION SEARCH IS
    ! NOT DONE ON EVALQP.OUT !?!
    ! If on top of GW
    if(associated(input%gw) .and. iqmt==1) then
      ! Save KS eigenvalues to use them later for renormalizing PMAT
      allocate(eval0(nstsv, nkptnr))
      eval0=evalsv
      ! If scissor correction is presented, one should nullify it
      input%xs%scissor=0.0d0
      ! Read QP Fermi energies and eigenvalues from file
      ! NOTE: QP evals are shifted by -efermi-eferqp with respect to KS evals
      ! NOTE: getevalqp sets mod_symmetry::nsymcrys to 1
      ! NOTE: getevalqp needs the KS eigenvalues as input
      nsymcrys_save = nsymcrys
      call getevalqp(nkptnr,vkl0,evalsv)
      nsymcrys = nsymcrys_save
      if(bi2d%isroot) then
        write(*,'("efermi=", f10.7)') efermi
        write(*,'("eferqp=", f10.7)') eferqp
        do ik=1, nkptnr
          write(*,'("k-point #", i6,":", 3f10.7)') ik, vkl0(:,ik)
          write(*,'(" state    Eks    E+efermi+eferqp    E+eferqp    E+efermi    E")')
          do i=input%gw%ibgw, input%gw%nbgw
            write(*,'(i6, 5(4x, f10.7))') i, eval0(i,ik), evalsv(i, ik)+efermi+eferqp,&
              & evalsv(i, ik)+eferqp, evalsv(i, ik)+efermi, evalsv(i, ik) 
          end do
        end do
      end if
      ! Set k and k'=k grid eigenvalues to QP energies
      evalsv0=evalsv
      if(bi2d%isroot) then
        write(unitout,'("Info(b_bse): Quasi particle energies are read from EVALQP.OUT")')
      end if
    else if(associated(input%gw) .and. iqmt /= 1) then 
      if(bi2d%isroot) then 
        write(*,'("Error(b_bse): BSE+GW only supported for 0 momentum transfer.")')
      end if
      call terminate
    end if
    call barrier

    ! Select relevant transitions for the construction
    ! of the BSE hamiltonian
    ! Also sets nkkp_bse, nk_bse 
    !   Note: Operates on mpiglobal
    call select_transitions(iqmt, serial=.false.)

    ! Only MPI ranks that have an associated 
    ! BLACS grid process proceed.
    if(bi2d%isactive) then 

      ! Determine "BSE"-gap, i.e. the lowest energy 
      ! KS transition which is considered.
      bsegap = 1.0d16
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(a1,iu,io,ik,ikq),&
      !$OMP& REDUCTION(min:bsegap)
      do a1 = 1, hamsize
        iu = smap(1,a1)
        io = smap(2,a1)
        ik = smap(3,a1)
        ikq = ikmapikq(ik, iqmt)
        bsegap = min(bsegap, evalsv(iu,ik)-evalsv0(io,ikq))
      end do
      !$OMP END PARALLEL DO

      if(bi2d%isroot) then
        write(unitout,*)
        write(unitout, '("Info(b_bse):&
          & bsegap, bsegap+scissor (eV):", E23.16,1x,E23.16)')&
          & bsegap*h2ev, (bsegap+sci)*h2ev
      end if

      ! Define global distributed Hamiltonian matrix.
      call new_dzmat(dham, hamsize, hamsize, bi2d)

      ! Write Info
      if(bi2d%isroot) then
        write(unitout,*)
        write(unitout, '("Info(b_bse): Assembling distributed BSE matrix")')
        if(fcoup .and. fti) then
          write(unitout, '(" Including coupling terms ")')
          write(unitout, '(" Using time inverted anti-resonant basis")')
          write(unitout, '(" Using squared EVP")')
        end if
        write(unitout, '("  RR/RA blocks of global BSE-Hamiltonian:")')
        write(unitout, '("  Shape=",i8)') hamsize
        write(unitout, '("  nk=", i8)') nk_bse
        write(unitout, '("  Distributing matrix to ",i3," processes")') bi2d%nprocs
        write(unitout, '("  Local matrix shape ",i6," x",i6)')&
          & dham%nrows_loc, dham%ncols_loc
      end if

      ! Assemble Hamiltonian matrix
      if(bi2d%isroot) call timesec(ts0)
      if(fcoup .and. fti) then 
        call setup_dis_ti_bse(dham, dcmat, dcpmat,  iqmt)
      else
        call setup_distributed_bse(dham, iqmt, .false., .false., bi2d)
      end if

      if(bi2d%isroot) then
        call timesec(ts1)
        write(unitout, '("All processes build their local matrix")')
        write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0
      end if
      
      if(fwp) then
        call dzmat_send2global_root(ham, dham, bi2d)
        if(bi2d%isroot) then
          if(fcoup .and. fti) then 
            call writecmplxparts("GlobalHamS", dble(ham), immat=aimag(ham))
          else
            call writecmplxparts("GlobalHam", dble(ham), immat=aimag(ham))
          end if
          deallocate(ham)
        end if
      end if

      ! Write Info
      if(bi2d%isroot) then
        if(fcoup .and. fti) then 
          write(unitout, '("Info(b_bse): Solving hermitian squared EVP")')
          write(unitout, '("Info(b_bse): Invoking scalapack routine PZHEEVX")')
        else 
          write(unitout, '("Info(b_bse): Diagonalizing RR Hamiltonian (TDA)")')
          write(unitout, '("Info(b_bse): Invoking scalapack routine PZHEEVX")')
        end if
      end if

      ! Eigenvectors are distributed
      ! must be NxN because Scalapack solver expects it so
      call new_dzmat(dbevecr, hamsize, hamsize, bi2d) 
      ! Eigenvalues are global
      allocate(bevalre(hamsize))
      bevalre = 0.0d0

      ! Find only eigenvalues relevant for requested 
      ! spectrum?
      efind = input%xs%bse%efind

      ! Diagonalize Hamiltonian (destroys the content of ham)
      if(bi2d%isroot) call timesec(ts0)

      if(efind) then

        if(fcoup .and. fti) then 
          v1=(max(wl, 0.0d0))**2
          v2=(wu)**2
        else
          v1=max(wl, 0.0d0)
          v2=wu
        end if

        call dhesolver(dham, bevalre, bi2d, dbevecr,&
         & v1=v1, v2=v2, found=nexc,&
         & eecs=input%xs%bse%eecs)

      else

        if(input%xs%bse%nexc == -1) then 
          i2 = hamsize
        else
          i2 = input%xs%bse%nexc
        end if
        i1 = 1

        call dhesolver(dham, bevalre, bi2d, dbevecr,&
         & i1=i1, i2=i2, found=nexc,&
         & eecs=input%xs%bse%eecs)

      end if

      if(fcoup .and. fti) then 
        ! Take square root of auxilliary eigenvalues
        ! to retrieve actual eigenvalues
        if(any(bevalre < 0.0d0)) then 
          write(*,*) "Error(b_bse): Negative squared EVP evals occured"
          write(*,'(E23.16)') bevalre
          call terminate
        end if
        bevalre = sqrt(bevalre)
      end if

      ! Deallocate BSE-Hamiltonian
      call del_dzmat(dham)

      if(bi2d%isroot) call timesec(ts1)

      if(bi2d%isroot) then
        if(efind) then
          if(fcoup .and. fti) then 
            write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
            write(unitout, '("  [",E10.3,",",E10.3,"]/H^2")') sqrt(v1), sqrt(v2)
            write(unitout, '("  [",E10.3,",",E10.3,"]/eV")') sqrt(v1)*h2ev, sqrt(v2)*h2ev
          else
            write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
            write(unitout, '("  [",E10.3,",",E10.3,"]/H")') v1, v2
            write(unitout, '("  [",E10.3,",",E10.3,"]/eV")') v1*h2ev, v2*h2ev
          end if
        else
          write(unitout, '("  ",i8," eigen solutions found in the interval:")') nexc
          write(unitout, '("  i1=",i8," i2=",i8)') i1, i2
        end if
        write(unitout, '("  Timing (in seconds)	   :", f12.3)') ts1 - ts0
        write(unitout,*)
      end if


      ! Store excitonic energies and wave functions to file
      if(associated(input%xs%storeexcitons)) then

        if(input%xs%storeexcitons%selectenergy) then 
          en1=input%xs%storeexcitons%minenergyexcitons
          en2=input%xs%storeexcitons%maxenergyexcitons
          if(bi2d%isroot) then 
            !write(*,*) "b_bse: en1, en2", en1, en2
          end if
          if(input%xs%storeexcitons%useev) then 
            en1=en1/h2ev
            en2=en2/h2ev
          end if
          call energy2index(hamsize, nexc, bevalre, en1, en2, iex1, iex2)
          if(bi2d%isroot) then 
            !write(*,*) "b_bse: iex1, iex2", iex1, iex2
          end if
        else
          iex1=input%xs%storeexcitons%minnumberexcitons
          iex2=input%xs%storeexcitons%maxnumberexcitons
        end if
        nreq=iex2-iex1+1
        if(bi2d%isroot) then 
          !write(*,*) "b_bse: iex1=", iex1, " iex2=", iex2, " nreq=", nreq
        end if

        if(nreq < 1 .or. nreq > nexc .or. iex1<1 .or. iex2<1 ) then
          if(bi2d%isroot) then 
            write(*,*) "Error(b_bse): storeexcitons index mismatch."
            write(*,*) "iex1, iex2, nreq, nex", iex1, iex2, nreq, nexc
          end if
          call terminate
        end if

        if(bi2d%isroot) then 
          write(unitout, '("Info(b_bse): Writing excition eigenvectors for index range=",2i8)') iex1, iex2
        end if

        if(fcoup .and. fti) then 

          call new_dzmat(dresvec, hamsize, nreq, bi2d)
          call new_dzmat(daresvec, hamsize, nreq, bi2d)

          if(bi2d%isroot) then 
            write(unitout, '("Info(b_bse): Generating resonant&
              & and anti-resonant exciton coefficients from&
              & auxilliary squared EVP eigenvectors (pos. E).")')
          end if
          call gendexevec(iex1, iex2, nexc, dcmat, dcpmat, dbevecr, bevalre,&
            & drvecp=dresvec, davecp=daresvec)

          write(unitout, '("Info(b_bse): Writing exciton eigenvectors to file.")')
          call putd_excitons(bevalre(iex1:iex2), drvec=dresvec, davec=daresvec,&
            & iqmt=iqmt, a1=iex1, a2=iex2)

          call del_dzmat(dresvec)
          call del_dzmat(daresvec)

        else 

          write(unitout, '("Info(b_bse): Writing exciton eigenvectors to file.")')
          call setview_dzmat(dbevecr, hamsize, nreq, 1, iex1)
          call putd_excitons(bevalre(iex1:iex2), dbevecr,&
            & iqmt=iqmt, a1=iex1, a2=iex2)
          call setview_dzmat(dbevecr, hamsize, nexc, 1, 1)

        end if

      end if

      ! Calculate oscillator strengths.
      ! Note: Deallocates eigenvectors
      call make_doszstren(doszsr)

      ! Every process gets a copy of the oscillator strength
      ! (actually only rank 0 writes them to file, but is is not much 
      !  memory and it make the setup for the spectrum calculation easier) 
      call dzmat_send2global_all(oszsr, doszsr, bi2d)

      if(bi2d%isroot) then
        ! Write excition energies and oscillator strengths to 
        ! text file. 
        write(unitout, '("Info(b_bse):&
          & Writing excition energies and oszillator strengths to text file.")')
        call writeoscillator(hamsize, nexc, -(bsegap+sci), bevalre, oszsr, iqmt=iqmt)
      end if

      ! Allocate arrays used in spectrum construction
      allocate(w(nw))
      ! Generate an evenly spaced frequency grid 
      call genwgrid(nw, input%xs%energywindow%intv,&
        & input%xs%tddft%acont, 0.d0, w_real=w)

      ! Allocate arrays used in spectrum construction
      if(bi2d%isroot) then
        allocate(symspectr(3,3,nw))
      end if

      ! Only process 0 gets an acctual output for symspectr
      if(fcoup .and. fti) then 
        call make_dist_spectrum_ti(nw, w, symspectr)
      else
        call make_dist_spectrum_tda(nw, w, symspectr)
      end if

      if(bi2d%isroot) then
        write(unitout, '("Writing derived quantities.")')
        ! Generate and write derived optical quantities
        call writederived(iqmt, symspectr, nw, w)
        write(unitout, '("  Derived quantities written.")')
      end if

      ! Clean up
      deallocate(bevalre, oszsr, w, evalsv)
      if(associated(input%gw)) deallocate(eval0)
      if(bi2d%isroot) then 
        deallocate(symspectr)
      end if

      ! Ranks that are on the BLACS grid signal that they are done
      call barrier

    else

      write(*, '("Warning(b_bse): Rank", i4, " is idle.")') mpiglobal%rank

      ! Ranks that are not on the BLACS grid wait 
      call barrier

    end if

  else ! not fscal and not rank 0

    write(*, '("Warning(b_bse): Rank", i4, " is idle.")') mpiglobal%rank

    ! Rank /= 0 wait for rank 0
    call barrier

  end if

contains

  !! SERIAL VERSIONS
  subroutine setup_ti_bse(smat, cmat, cpmat, iqmt)

    use m_sqrtzmat

    ! I/O
    integer(4), intent(in) :: iqmt
    complex(8), allocatable, intent(inout) :: smat(:,:)
    complex(8), allocatable, intent(inout) :: cmat(:,:)
    complex(8), allocatable, intent(inout) :: cpmat(:,:)

    ! Local 
    complex(8), allocatable :: rrmat(:, :)
    complex(8), allocatable :: ramat(:, :)
    complex(8), allocatable :: auxmat(:, :)

    integer(4) :: info, lwork
    integer(4), allocatable :: ipiv(:)
    complex(8), allocatable :: work(:)

    real(8) :: ts0, ts1, t1, t0
    integer(4) :: i, j

    real(8) :: evals(hamsize)
    logical :: fwp

    fwp = input%xs%bse%writeparts

    evals = 0.0d0

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Setting up matrices for squared EVP")')
      call timesec(ts0)
    end if

    ! Get RR part of BSE Hamiltonian

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Setting up RR Block of orignial BSE")')
      call timesec(t0)
    end if
    allocate(rrmat(hamsize,hamsize))
    call setup_bse(rrmat, iqmt, .false., .false.)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    ! Get RA^{ti} part of BSE Hamiltonian
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse):&
        & Setting up RA^{ti} Block of orignial BSE")')
      call timesec(t0)
    end if
    allocate(ramat(hamsize,hamsize))
    call setup_bse(ramat, iqmt, .true., .true.)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    ! Make combination matrices

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Setting up RR+RA and RR-RA matrices")')
      call timesec(t0)
    end if
    ! RR - RA^{it}
    allocate(cpmat(hamsize,hamsize))
    do j = 1, hamsize
      do i = 1, hamsize
        cpmat(i,j) = rrmat(i,j) - ramat(i,j)
      end do
    end do
    ! RR + RA^{it}
    deallocate(rrmat)
    allocate(cmat(hamsize, hamsize))
    do j = 1, hamsize
      do i = 1, hamsize
        cmat(i,j) = cpmat(i,j) + 2.0d0*ramat(i,j)
      end do
    end do
    !deallocate(ramat)

    ! Check positive definitness of (A+B)
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Checking positve definitness of RR+RA")')
      call timesec(t0)
    end if
    ramat = cmat
    if(fwp) then
      call writecmplxparts("nd_apb_mat", remat=dble(ramat), immat=aimag(ramat))
    end if
    call hesolver(ramat, evals)

    !write(*,*) "writing apm evals"
    if(fwp) then
      call writecmplxparts("nd_apb_evals", revec=evals, veclen=size(evals))
    end if

    if(any(evals < 0.0d0)) then 
      write(*,*) "Error(setup_ti_bse): A+B matrix is not positive definit"
      write(*,'(E10.3)') evals
      call terminate
    end if
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  RR+RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if
    deallocate(ramat)

    ! Take the square root of (A-B) (currently in cpmat)
    ! Note: It is assumed to be positive definit.
    if(fwp) then 
      call writecmplxparts("nd_amb_mat", remat=dble(cpmat), immat=aimag(cpmat))
    end if

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Taking square root of RR-RA matrix")')
      call timesec(t0)
    end if
    call sqrtzmat_hepd(cpmat)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  RR-RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    if(fwp) then
      call writecmplxparts("nd_sqrtamb_mat", remat=dble(cpmat), immat=aimag(cpmat))
    end if

    ! Construct S Matrix
    ! S = (A-B)^{1/2} (A+B) (A-B)^{1/2}

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Constructing S matrix")')
      call timesec(t0)
    end if

    allocate(auxmat(hamsize,hamsize))
    call zgemm('N','N', hamsize, hamsize, hamsize, zone, cpmat, hamsize,&
      & cmat, hamsize, zzero, auxmat, hamsize)

    allocate(smat(hamsize, hamsize))
    call zgemm('N','N', hamsize, hamsize, hamsize, zone, auxmat, hamsize,&
      & cpmat, hamsize, zzero, smat, hamsize)

    !write(*,*) "printing s mat"
    if(fwp) then 
      call writecmplxparts('nd_s_mat', remat=dble(smat), immat=aimag(smat))
    end if

    deallocate(auxmat)

    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    ! Cmat = (A-B)^1/2
    do j = 1, hamsize
      do i = 1, hamsize
        cmat(i,j) = cpmat(i,j)
      end do
    end do

    ! Cpmat = (A-B)^-1/2

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_ti_bse): Inverting (RR-RA)^1/2 matrix")')
      call timesec(t0)
    end if
    allocate(ipiv(hamsize))
    call zgetrf(hamsize, hamsize, cpmat, hamsize, ipiv, info)
    if(info /= 0) then
      write(*,'("Info(setup_ti_bse): ZGETRF has returned non-zero info.")')
      if(info < 0) then
        write(*,'("  ZGETRF: Incorrect input.")')
      else
        write(*,'("  ZGETRF: U is exactly singular.")')
      end if
      call terminate
    end if
    lwork = -1
    allocate(work(3))
    call zgetri(hamsize, cpmat, hamsize, ipiv, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call zgetri(hamsize, cpmat, hamsize, ipiv, work, lwork, info)
    if(info /= 0) then
      write(*,'("dzinvert (ERROR): ZGETRI has returned non-zero info.")')
      if(info < 0) then
        write(*,'("  ZGETRF: Incorrect input.")')
      else
        write(*,'("  ZGETRF: Matrix singular.")')
      end if
      call terminate
    end if
    deallocate(work)
    deallocate(ipiv)

    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    if(mpiglobal%rank == 0) then 
      call timesec(ts1)
      write(unitout, '("Info(setup_ti_bse): Total time needed",f12.3,"s")') ts1-ts0
    end if

  end subroutine setup_ti_bse

  ! NOTE: only for qmt = 0
  subroutine setup_full_bse(ham, iqmt)

    ! I/O
    complex(8), intent(inout) :: ham(:, :)
    integer(4), intent(in) :: iqmt

    ! Timings
    real(8) :: ts0, ts1

    integer(4) :: i, j

    call timesec(ts0)
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_full_bse): Setting up full hamiltonian")')
    end if

    ! RR
    call setup_bse(ham(1:hamsize,1:hamsize), iqmt, .false., .false.)

    ! RA
    call setup_bse(ham(1:hamsize,hamsize+1:hamsize*2), iqmt, .true., .false.)
! ham(1:hamsize,hamsize+1:hamsize*2) = zzero

    ! AR
    ! Note: AR part is the negative complex conjugate of RA even if qmt /= 0
    ham(hamsize+1:2*hamsize, 1:hamsize) = -conjg(ham(1:hamsize,hamsize+1:hamsize*2))
! ham(hamsize+1:2*hamsize, 1:hamsize) = zzero

    ! AA
    ! Note: AA part is the negative complex conjugate of RR ONLY if qmt /= 0
    ham(hamsize+1:2*hamsize, hamsize+1:2*hamsize) = -conjg(ham(1:hamsize,1:hamsize))

    call timesec(ts1)
    write(unitout, '(" Matrix build.")')
    write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0

  end subroutine setup_full_bse

  subroutine makeoszillatorstrength(oszstrr, oszstra)
    use m_getpmat
    use m_invertzmat
    
    implicit none

    ! I/O
    complex(8), intent(out) :: oszstrr(nexc,3)
    complex(8), intent(out), optional :: oszstra(nexc, 3)

    ! Local
    real(8) :: t1, t0, ts0, ts1
    complex(8), allocatable :: projmat(:,:), xpy(:,:), rbarmat(:,:)
    integer(4) :: io, iu, ioabs, iuabs, ik, iknr
    integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 
    integer(4) :: a1, lambda, igqmt, nopt
    real(8) :: potcl
    
    write(unitout, '("Info(b_bse:makeos): Making oszillator strengths.")')
    write(unitout, '("Info(b_bse:makeos): iqmt=", i4)') iqmt
    write(unitout, '("Info(b_bse:makeos): vqlmt=", 3E11.3)') vqlmt(:,iqmt)

    call timesec(ts0)

    if(iqmt == 1) then 

      write(unitout, '("Info(b_bse:makeos): Using zero momentum transfer formalismus.")')
      write(unitout, '("Info(b_bse:makeos):&
        & Making oszillator strengths using position operator matrix elements.")')
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Building position operator matrix elements using momentum matrix elements  !
      ! and transition energies. If on top of GW, renormalize the p mat elements.  !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      nopt=3
      allocate(projmat(hamsize, nopt))
      projmat=zzero
      call setup_rmat(projmat)

    else

      write(unitout, '("Info(b_bse:makeos): Making oszillator strengths using&
        & plane wave matrix elements.")')
      nopt=1
      allocate(projmat(hamsize,nopt))
      ! Generate plane wave matrix elements for G=Gqmt and q=qmt
      igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
      call setup_pwmat(projmat, iqmt, igqmt)

    end if

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! TDA case: Build resonant oscillator strengths                              !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if(.not. fcoup) then 
      write(unitout, '("  Building resonant oscillator strengths.")')
      call timesec(t0)
      oszstrr = zzero
      !! Resonant oscillator strengths
      ! qmt=0 case:
      ! t^R_{\lambda,i} = <X_\lambda|\tilde{R}^{i*}> =
      !   \Sum_{a} X^H_{\lambda, a} \tilde{R}^*_{a,i}
      ! qmt/=0 case:
      ! t^R_{\lambda}(G,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
      !   \Sum_{a} X^H_{\lambda, a} \tilde{M}^*_{a}(G,qmt)
      call zgemm('c','n', nexc, nopt, hamsize,&
        & zone, bevecr(1:hamsize,1:nexc), hamsize, conjg(projmat), hamsize,&
        & zzero, oszstrr(1:nexc,1:nopt), nexc)
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Full case: Build left and right oscillator strengths                       !
    !            New version                                                     !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if(fcoup .and. .not. fti) then 

      allocate(rbarmat(2*hamsize,3))
      rbarmat(1:hamsize,1:3) = conjg(projmat)
      rbarmat(hamsize+1:2*hamsize,1:3) = -projmat
      deallocate(projmat)

      write(unitout, '("  Building t_r oscillator strengths.")')
      call timesec(t0)

      !! Oscillator strengths from right eigenvectors
      ! t_r_{\lambda,i} = <P_\lambda|\bar{R}_i> =
      !   =  \Sum_{a} P^H_{\lambda,a} \bar{R}_{a,i}
      call zgemm('c','n', nexc, 3, 2*hamsize,&
        & zone, bevecr(1:2*hamsize,1:nexc), 2*hamsize, rbarmat, 2*hamsize,&
        & zzero, oszstrr, nexc)

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0


      write(unitout, '("  Inverting right EV matrix.")')
      call timesec(t0)
      call zinvert(bevecr)
      call timesec(t1)

      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Building t_l oscillator strengths.")')
      call timesec(t0)

      rbarmat(hamsize+1:2*hamsize,:) = -rbarmat(hamsize+1:2*hamsize,:)
      !! Oscillator strengths from left eigenvectors
      ! t_l_{\lambda,i} = < P^-1_\lambda | \bar{R}^s_i> =
      !   \Sum_{a} P^{-1}_{\lambda, a} \tilde{R}^s_{a,i}
      call zgemm('n','n', nexc, 3, 2*hamsize,&
        & zone, bevecr(1:nexc,1:2*hamsize), 2*hamsize,&
        & rbarmat, 2*hamsize, zzero, oszstra, nexc)

      deallocate(rbarmat)
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Full case: Build left and right oscillator strengths                       !
    !            Old version                                                     !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if(fcoup .and. .not. fti .and. .false.) then 
      write(unitout, '("  Building resonant oscillator strengths.")')
      call timesec(t0)
      !! Resonant oscillator strengths
      ! t^R_{\lambda,i} = < \tilde{R}^{i*} | X_\lambda> =
      !   ( \Sum_{a} \tilde{R}^T_{i, a} X_{a, \lambda} )^T =
      !     \Sum_{a} X^T_{\lambda, a} \tilde{R}_{a,i}
      call zgemm('t','n', nexc, 3, hamsize,&
        & zone, bevecr(1:hamsize,1:nexc), hamsize, projmat, hamsize, zzero, oszstrr, nexc)
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      write(unitout, '("  Building anti-resonant oscillator strengths.")')
      call timesec(t0)
      !! Anti-resonant oscillator strengths
      ! t^A_{\lambda,i} = < \tilde{R}^{i} | Y_\lambda> =
      !   ( \Sum_{a} \tilde{R}^\dag_{i, a} Y_{a, \lambda} )^T =
      !     \Sum_{a} Y^T_{\lambda, a} \tilde{R}^*_{a,i}
      call zgemm('t','n', nexc, 3, hamsize,&
        & zone, bevecr(hamsize+1:2*hamsize,1:nexc), hamsize,&
        & conjg(projmat), hamsize, zzero, oszstra, nexc)
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! TI case: Build oscillator strength                                         !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if(fcoup .and. fti) then 
      write(unitout, '("  Building (X+Y) from squared EVP EVs.")')
      call timesec(t0)

      ! Interested in X^+ + Y^+, so we rescale the 
      ! auxiliary eigenvectors Z by the square root of the eigenvalues
      ! so that 
      ! (X+Y)_{a,lambda} = \Sum_{a'} (A-B)^{1/2}_{a,a'} * E^{-1/2}_lambda * Z_{a',lambda}
      do lambda = 1, nexc
        do a1 = 1, hamsize
          bevecaux(a1, lambda) = bevecaux(a1, lambda) / sqrt(bevalre(lambda))
        end do
      end do
      allocate(xpy(hamsize, nexc))
      call zgemm('N','N',hamsize, nexc, hamsize, zone, cmat, hamsize,&
        & bevecaux(:,1:nexc), hamsize, zzero, xpy, hamsize)
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0

      write(unitout, '("  Building oscillator strengths for time inverted ar basis.")')
      call timesec(t0)
      !! Oscillator strengths
      ! qmt=0
      ! t_{\lambda,i} = < (| X_\lambda>+| Y_\lambda>)| \tilde{R}^{i*}> =
      !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{R}^*_{a,i}
      ! qmt/=0
      ! t_{\lambda}(G,qmt) = < (| X_\lambda>+| Y_\lambda>)| \tilde{M}^*(G,qmt)> =
      !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{M}^*_{a}
      call zgemm('c','n', nexc, nopt, hamsize,&
        & zone, xpy(1:hamsize,1:nexc), hamsize, conjg(projmat), hamsize,&
        & zzero, oszstrr(1:nexc,1:nopt), nexc)

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      deallocate(xpy)
    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    if(allocated(projmat)) deallocate(projmat)
    call timesec(ts1)
    write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0
  end subroutine makeoszillatorstrength

  subroutine makespectrum_ti(nfreq, freq, spectrum)
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(:,:,:)

    ! Local
    integer(4) :: i, j, o1, o2, nopt
    real(8) :: t1, t0, ts0, ts1
    complex(8) :: brd
    complex(8), allocatable :: ns_spectr(:,:)
    complex(8), allocatable :: enw(:,:), tmat(:,:)

    write(unitout, '("Info(b_bse:makesp_ti):&
     & Making spectrum using formula for coupling with time inverted ar basis.")')
    if(iqmt == 1) then 
      write(unitout, '("Info(b_bse:makesp_ti):&
       & Using formalism for vaninshing momentum transfer.")')
      if(input%xs%dfoffdiag) then
        write(unitout, '("Info(b_bse:makesp_ti): Including off diagonals.")')
      end if
    else
      write(unitout, '("Info(b_bse:makesp_ti):&
       & Using formalism for finite momentum transfer. iqmt=",i4)') iqmt
    end if

    ! Total time for spectrum construction
    call timesec(ts0)

    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make energy denominator for each frequency !
    ! (resonant & anti-resonant)                 !
    !++++++++++++++++++++++++++++++++++++++++++++!
    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)

    ! enw_{w, \lambda} = 1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta)
    allocate(enw(nfreq, nexc))

    ! Broadening factor
    brd = zi*input%xs%broad

    !$OMP PARALLEL DO &
    !$OMP& COLLAPSE(2),&
    !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
    do j = 1, nexc
      do i = 1, nfreq
        enw(i,j) = zone/(bevalre(j)-freq(i)-brd)&
                &+ zone/(bevalre(j)+freq(i)+brd)
      end do
    end do
    !$OMP END PARALLEL DO

    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    !++++++++++++++++++++++++++++++++++++++++++++!
    
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Helper matrix build form oscillator        !
    ! strengths                                  !
    !++++++++++++++++++++++++++++++++++++++++++++!
    write(unitout, '("  Making helper matrix tmat.")')
    call timesec(t0)

    if(iqmt == 1) then 

      if(input%xs%dfoffdiag) then
        nopt = 9
      else
        nopt = 3
      end if

      ! tmat_{\lambda, j} = t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! where j combines the Cartesian directions
      allocate(tmat(nexc,nopt))
      do j = 1, nopt
        if(input%xs%dfoffdiag) then
          o2 = (j-1)/3 + 1
          o1 = j-(o2-1)*3
        else
          o2 = j
          o1 = j
        end if
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nexc
          tmat(i,j) = conjg(oszsr(i,o1))*oszsr(i,o2)
        end do
        !$OMP END PARALLEL DO
      end do

    else

      nopt = 1

      ! tmat_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
      allocate(tmat(nexc,nopt))
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i)
      do i = 1, nexc
        tmat(i,1) = conjg(oszsr(i,1))*oszsr(i,1)
      end do
      !$OMP END PARALLEL DO

    end if

    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    !++++++++++++++++++++++++++++++++++++++++++++!
         
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make non-lattice-symmetrized spectrum      !
    !++++++++++++++++++++++++++++++++++++++++++++!
    write(unitout, '("  Calculating spectrum.")')
    call timesec(t0)

    allocate(ns_spectr(nfreq,nopt))

    ! qmt=0 case:
    ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
    !   i.e. nsspectr_{w,j} = 
    !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
    !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
    ! qmt/=0 case:
    ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
    !   i.e. nsspectr_{w,j} = 
    !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
    !       t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
    call zgemm('N','N', nfreq, nopt, nexc, zone, enw, nfreq,&
      & tmat, nexc, zzero, ns_spectr, nfreq)
    !++++++++++++++++++++++++++++++++++++++++++++!

    ! Helper no longer needed
    deallocate(tmat, enw)

    ! Post process non-lattice-symmetrized spectrum
    call finalizespectrum(ns_spectr, spectrum)
    deallocate(ns_spectr)

    ! Total time for spectrum construction
    call timesec(ts1)
    write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0
  end subroutine makespectrum_ti

  subroutine makespectrum_tda(nfreq, freq, spectrum)
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(3,3,nfreq)

    ! Local
    integer(4) :: i, j, o1, o2, nopt
    real(8) :: t1, t0, ts0, ts1
    complex(8), allocatable :: ns_spectr(:,:)
    complex(8), allocatable :: enw(:,:), enwr(:,:), enwa(:,:), tmatr(:,:)
    complex(8) :: brd

    write(unitout, '("Info(b_bse:makesp_tda):&
     & Making spectrum using TDA formula.")')
    if(iqmt == 1) then 
      write(unitout, '("Info(b_bse:makesp_tda):&
       & Using formalism for vaninshing momentum transfer.")')
      if(input%xs%dfoffdiag) then
        write(unitout, '("Info(b_bse:makesp_tda): Including off diagonals.")')
      end if
    else
      write(unitout, '("Info(b_bse:makesp_tda):&
       & Using formalism for finite momentum transfer. iqmt=",i4)') iqmt
    end if
    if(fti) then 
      write(unitout, '("Info(b_bse:makesp_tda): Using fomula for ti ar basis.")')
    end if

    if(iqmt /= 1 .and. .not. fti) then 
      write(unitout, '("Info(b_bse:makesp_tda):&
        & Finite momentum transfer only supported for ti ar basis.")')
      call terminate
    end if

    ! Total time for spectrum construction
    call timesec(ts0)

    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make energy denominator for each frequency !
    ! (resonant & anti-resonant)                 !
    !++++++++++++++++++++++++++++++++++++++++++++!
    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)

    ! Broadening factor
    brd = zi*input%xs%broad

    if(.not. fti) then 

      ! enwr_{w, \lambda} = 1/(E_\lambda - w - i\delta)
      allocate(enwr(nfreq, nexc))
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
      do j = 1, nexc
        do i = 1, nfreq
          enwr(i,j) = zone/(bevalre(j)-freq(i)-brd)
        end do
      end do
      !$OMP END PARALLEL DO

      ! enwa_{w, \lambda} = 1/(E_\lambda + w + i\delta)
      allocate(enwa(nfreq, nexc))
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
      do j = 1, nexc
        do i = 1, nfreq
          enwa(i,j) = zone/(bevalre(j)+freq(i)+brd)
        end do
      end do
      !$OMP END PARALLEL DO

    else

      ! enw_{w, \lambda} = 1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta)
      allocate(enw(nfreq, nexc))

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
      do j = 1, nexc
        do i = 1, nfreq
          enw(i,j) = zone/(bevalre(j)-freq(i)-brd)&
                  &+ zone/(bevalre(j)+freq(i)+brd)
        end do
      end do
      !$OMP END PARALLEL DO

    end if

    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    !++++++++++++++++++++++++++++++++++++++++++++!
    
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Helper matrix build form resonant          !
    ! oszillator strenghs                        !
    !++++++++++++++++++++++++++++++++++++++++++++!
    write(unitout, '("  Making helper matrix tmatr.")')
    call timesec(t0)

    if(iqmt == 1) then 

      if(input%xs%dfoffdiag) then
        nopt = 9
      else
        nopt = 3
      end if

      ! tmatr_{\lambda, j} = t^{R*}_{\lambda, o1_j} t^R_{\lambda, o2_j} 
      ! where j combines the 2 cartesian directions
      allocate(tmatr(nexc,nopt))
      do j = 1, nopt
        if(input%xs%dfoffdiag) then
          o2 = (j-1)/3 + 1
          o1 = j-(o2-1)*3
        else
          o2 = j
          o1 = j
        end if
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nexc
          tmatr(i,j) = conjg(oszsr(i,o1))*oszsr(i,o2)
        end do
        !$OMP END PARALLEL DO
      end do

    else
      
      nopt = 1
      ! tmat_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
      allocate(tmatr(nexc,nopt))
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i)
      do i = 1, nexc
        tmatr(i,nopt) = conjg(oszsr(i,nopt))*oszsr(i,nopt)
      end do
      !$OMP END PARALLEL DO

    end if

    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    !++++++++++++++++++++++++++++++++++++++++++++!
         
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make non-lattice-symmetrized spectrum      !
    !++++++++++++++++++++++++++++++++++++++++++++!
    write(unitout, '("  Calculating spectrum.")')
    call timesec(t0)

    allocate(ns_spectr(nfreq,nopt))

    if(.not. fti) then 
      ! nsspectr_{w,j} = \Sum_{\lambda} enwr_{w,\lambda} tmatr_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} 1/(E_\lambda - w - i\delta)
      !       (t^R)^H_{o1_j,\lambda} t^R_{\lambda, o2_j} 
      call zgemm('N','N', nfreq, nopt, nexc, zone, enwr, nfreq,&
        & tmatr, nexc, zzero, ns_spectr, nfreq)

      ! nsspectr_{w,j} += \Sum_{\lambda} enwa_{w,\lambda} tmatr^*_{\lambda, j}
      !   i.e. nsspectr_{w,j} = nsspectr_{w,j} +
      !     \Sum_{\lambda} 1/(E_\lambda + w + i\delta)
      !       ( (t^R)^H_{o1_j,\lambda} t^R_{\lambda, o2_j} )^*
      call zgemm('N','N', nfreq, nopt, nexc, zone, enwa, nfreq,&
        & conjg(tmatr), nexc, zone, ns_spectr, nfreq)
    else
      ! qmt=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! qmt/=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
      call zgemm('N','N', nfreq, nopt, nexc, zone, enw, nfreq,&
        & tmatr, nexc, zzero, ns_spectr, nfreq)
    end if

    !++++++++++++++++++++++++++++++++++++++++++++!

    ! Helper no longer needed
    if(allocated(enwr)) deallocate(enwr)
    if(allocated(enwa)) deallocate(enwa)
    if(allocated(enw)) deallocate(enw)
    if(allocated(tmatr)) deallocate(tmatr)

    ! Postprocess non-lattice-symmetrized spectrum
    call finalizespectrum(ns_spectr, spectrum)
    deallocate(ns_spectr)

    ! Total time for spectrum construction
    call timesec(ts1)
    write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0
  end subroutine makespectrum_tda

  subroutine makespectrum_full(nfreq, freq, spectrum)
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(3,3,nfreq)

    ! Local
    integer(4) :: i, j, o1, o2, nopt
    real(8) :: t1, t0, ts0, ts1
    complex(8), allocatable :: ns_spectr(:,:)
    complex(8), allocatable :: overlap(:,:), invoverlap(:,:)
    complex(8), allocatable :: tplus(:,:), tmat(:,:), enw(:,:), tminus(:,:), op(:,:)

    write(unitout, '("Info(b_bse:makesp): Making spectrum using general formula.")')
    if(input%xs%dfoffdiag) then
      write(unitout, '("  Including off diagonals.")')
    end if
    ! Total spectrum construction timer
    call timesec(ts0)

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
    write(unitout, '("  Making t- and t+.")')
    call timesec(t0)
    allocate(tminus(nexc,3), op(nexc,3))
    do j=1, 3
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i)
      do i=1, nexc
        tminus(i,j) = oszsr(i,j)-oszsa(i,j)
      end do
      !$OMP END PARALLEL DO
    end do
    do j=1, 3
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i)
      do i=1, nexc
        op(i,j) = conjg(oszsr(i,j)+oszsa(i,j))
      end do
      !$OMP END PARALLEL DO
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)
    ! Make energy denominator for each frequency 
    allocate(enw(nfreq, nexc))
    !$OMP PARALLEL DO &
    !$OMP& COLLAPSE(2),&
    !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
    do i = 1, nexc
      do j = 1, nfreq
        enw(j,i) = zone/(bevalre(i)-freq(j)-zi*input%xs%broad)
      end do
    end do
    !$OMP END PARALLEL DO
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
      nopt = 9
    else
      nopt = 3
    end if
    allocate(tmat(nexc,nopt))
    do j = 1, nopt
      if(input%xs%dfoffdiag) then
        o2 = (j-1)/3 + 1
        o1 = j-(o2-1)*3
      else
        o2 = j
        o1 = j
      end if
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i)
      do i= 1, nexc
        tmat(i,j) = tminus(i,o1)*tplus(i,o2)
      end do
      !$OMP END PARALLEL DO
    end do
    ! tplus and tminus not needed anymore
    deallocate(tplus, tminus)
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
         
    write(unitout, '("  Calculating spectrum.")')
    ! Make non-lattice-symmetrized spectrum
    call timesec(t0)
    allocate(ns_spectr(nfreq,nopt))
    ! \epsilon^M_{w,ij} \prop \Sum_\lambda E^-1_{w,\lambda} tmat_{\lambda,ij}
    call zgemm('N','N', nfreq, nopt, nexc, zone, enw, nfreq,&
      & tmat, nexc, zzero, ns_spectr, nfreq)
    ! Helper no longer needed
    deallocate(tmat,enw)

    ! Postprocess non-lattice-symmetrized spectrum
    call finalizespectrum(ns_spectr, spectrum)
    deallocate(ns_spectr)

    call timesec(ts1)
    write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0

  end subroutine makespectrum_full

  subroutine makespectrum_full_lr(nfreq, freq, spectrum)
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(3,3,nfreq)

    ! Local
    integer(4) :: i, j, o1, o2, nopt
    real(8) :: t1, t0, ts0, ts1
    complex(8), allocatable :: ns_spectr(:,:)
    complex(8), allocatable :: tmat(:,:), enw(:,:)

    write(unitout, '("Info(b_bse:makesp): Making spectrum using general formula&
      & with left and right EVs.")')
    if(input%xs%dfoffdiag) then
      write(unitout, '("  Including off diagonals.")')
    end if
    ! Total spectrum construction timer
    call timesec(ts0)

    write(unitout, '("  Making energy denominators ENW.")')
    call timesec(t0)
    ! Make energy denominator for each frequency 
    allocate(enw(nfreq, nexc))
    !$OMP PARALLEL DO &
    !$OMP& COLLAPSE(2),&
    !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
    do i = 1, nexc
      do j = 1, nfreq
        enw(j,i) = zone/(bevalre(i)-freq(j)-zi*input%xs%broad)
      end do
    end do
    !$OMP END PARALLEL DO
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0

    write(unitout, '("  Making helper matrix tmat.")')
    ! tmat_{\lambda,ij} = t_r^*_{\lambda,i}*t_l_{\lambda,j}
    call timesec(t0)
    if(input%xs%dfoffdiag) then
      nopt = 9
    else
      nopt = 3
    end if
    allocate(tmat(nexc,nopt))
    do j = 1, nopt
      if(input%xs%dfoffdiag) then
        o2 = (j-1)/3 + 1
        o1 = j-(o2-1)*3
      else
        o2 = j
        o1 = j
      end if
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i)
      do i= 1, nexc
        tmat(i,j) = conjg(oszsr(i,o1))*oszsa(i,o2)
      end do
      !$OMP END PARALLEL DO
    end do
    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
         
    write(unitout, '("  Calculating spectrum.")')

    ! Make non-lattice-symmetrized spectrum
    call timesec(t0)
    allocate(ns_spectr(nfreq,nopt))
    ! \epsilon^M_{w,ij} \prop \Sum_\lambda E^-1_{w,\lambda} tmat_{\lambda,ij}
    call zgemm('N','N', nfreq, nopt, nexc, zone, enw, nfreq,&
      & tmat, nexc, zzero, ns_spectr, nfreq)
    ! Helper no longer needed
    deallocate(tmat,enw)

    ! Postprocess non-lattice-symmetrized spectrum
    ! -\bar{P}_00(q->0) is needed, + was constructed
    call finalizespectrum(ns_spectr, spectrum)
    deallocate(ns_spectr)

    call timesec(ts1)
    write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0

  end subroutine makespectrum_full_lr

  subroutine finalizespectrum(nsp, sp)
    use modxs, only: sptclg

    complex(8), intent(inout) :: nsp(:,:)
    complex(8), intent(out) :: sp(:,:,:)

    complex(8), allocatable :: buf(:,:,:)
    integer(4) :: nfreq, nopt, i, j, o1, o2, igqmt
    real(8) :: pref, t0, t1
    logical :: foff

    ! Check input
    nfreq = size(nsp,1) 
    nopt = size(nsp,2)
    if(iqmt==1) then 
      if(nopt == 9) then 
        foff = .true.
      else if(nopt == 3) then
        foff = .false.
      else
        write(*,*) "Error (finalizespectrum): nopt invalid", nopt
        call terminate
      end if
    else
      if(nopt /= 1) then 
        write(*,*) "Error (finalizespectrum): nopt invalid", nopt, " iqmt=",iqmt
        call terminate
      end if
      foff=.false.
    end if

    write(unitout, '("Info(b_bse:finalizesp): Finalizing spectrum.")')
    call timesec(t0)

    ! Adjusting prfactor and add 1 to diagonal elements
    if(iqmt==1) then 
      pref = 2.d0*4.d0*pi/omega/nk_bse
    else
      igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
      pref = 2.0d0*sptclg(igqmt,iqmt)**2/omega/nk_bse
    end if
    do j = 1, nopt
      if(foff) then
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nfreq
          if(j == 1 .or. j == 5 .or. j == 9) then 
            nsp(i,j) = nsp(i,j)*pref+zone
          else
            nsp(i,j) = nsp(i,j)*pref
          end if
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nfreq
          nsp(i,j) = nsp(i,j)*pref+zone
        end do
        !$OMP END PARALLEL DO
      end if
    end do

    if(iqmt == 1) then 

      ! Write to buffer for symmetry routine
      allocate(buf(3,3,nfreq))
      buf=zzero
      do j = 1, nopt
        if(foff) then
          o2 = (j-1)/3 + 1
          o1 = j-(o2-1)*3
        else
          o2 = j
          o1 = j
        end if
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nfreq
          buf(o1,o2,i) = nsp(i,j)
        end do
        !$OMP END PARALLEL DO
      end do
    
      ! Symmetrize spectrum 
      sp = zzero
      do o1=1,3
        do o2=1,3
          ! Symmetrize the macroscopic dielectric tensor
          call symt2app(o1, o2, nfreq, symt2, buf, sp(o1,o2,:))
        end do 
      end do

    else

      sp = zzero
      sp(1,1,:) = nsp(:,1)

    end if

    call timesec(t1)
    write(unitout, '("    Time needed",f12.3,"s")') t1-t0
  end subroutine finalizespectrum

  !! DISTRIBUTED VERSIONS

  subroutine setup_dis_ti_bse(smat, cmat, cpmat, iqmt)

    use m_sqrtzmat
    use m_invertzmat

    ! I/O
    integer(4), intent(in) :: iqmt
    type(dzmat), intent(inout) :: smat
    type(dzmat), intent(inout) :: cmat
    type(dzmat), intent(inout) :: cpmat

    ! Local 
    type(dzmat) :: rrmat
    type(dzmat) :: ramat
    type(dzmat) :: auxmat

    real(8) :: ts0, ts1, t1, t0
    integer(4) :: i, j, ig, jg

    real(8) :: evals(hamsize)

    ! Test writeout
    complex(8), allocatable :: localmat(:,:)
    logical fwp 

    fwp = input%xs%bse%writeparts

    evals = 0.0d0

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse):&
        & Setting up distributed matrices for squared EVP")')
      call timesec(ts0)
    end if

    ! Get RR part of BSE Hamiltonian
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse): Setting up RR Block of orignial BSE")')
      call timesec(t0)
    end if
    call new_dzmat(rrmat, hamsize, hamsize, bi2d)
    call setup_distributed_bse(rrmat, iqmt, .false., .false., bi2d)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    ! Get RA^{ti} part of BSE Hamiltonian
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse):&
        & Setting up RA^{ti} Block of orignial BSE")')
      call timesec(t0)
    end if
    call new_dzmat(ramat, hamsize, hamsize, bi2d)
    call setup_distributed_bse(ramat, iqmt, .true., .true., bi2d)

    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    ! Make combination matrices
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse): Setting up RR+RA and RR-RA matrices")')
    !  call timesec(t0)
    end if

    ! RR - RA^{it}
    call new_dzmat(cpmat, hamsize, hamsize, bi2d)
    do j = 1, cpmat%ncols_loc
      do i = 1, cpmat%nrows_loc
        cpmat%za(i,j) = rrmat%za(i,j) - ramat%za(i,j)
      end do
    end do
    ! RR + RA^{it}
    call del_dzmat(rrmat)
    call new_dzmat(cmat, hamsize, hamsize, bi2d)
    do j = 1, cmat%ncols_loc
      do i = 1, cmat%nrows_loc
        cmat%za(i,j) = cpmat%za(i,j) + 2.0d0*ramat%za(i,j)
      end do
    end do
    !call del_dzmat(ramat)

    ! Check positive definitness of (A+B)

    if(fwp) then 
      call dzmat_send2global_root(localmat, cmat, bi2d)
      if(mpiglobal%rank == 0) then
        call writecmplxparts('apb_mat', remat=dble(localmat), immat=aimag(localmat))
      end if
    end if

    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse): Checking positve definitness of RR+RA")')
      call timesec(t0)
    end if
    ramat%za = cmat%za
    call dhesolver(ramat, evals, bi2d)
    if(fwp) then
      if(mpiglobal%rank == 0) then 
        write(*,*) "Writing evals for A+B"
        call writecmplxparts('apb_evals', revec=evals, veclen=size(evals))
      end if
    end if
    if(any(evals < 0.0d0)) then 
      write(*,*) "Error(setup_dis_ti_bse): RR+RA matrix is not positive definit"
      write(*,'(E10.3)') evals
      call terminate
    end if
    call del_dzmat(ramat)

    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  RR+RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    ! Take the square root of (A-B) (currently in cpmat)
    if(fwp) then 
      call dzmat_send2global_root(localmat, cpmat, bi2d)
      if(mpiglobal%rank == 0) then
        call writecmplxparts('amb_mat', remat=dble(localmat), immat=aimag(localmat))
      end if
    end if

    ! Note: It is assumed to be positive definit.
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse): Taking square root of RR-RA matrix")')
      call timesec(t0)
    end if
    call sqrtdzmat_hepd(cpmat, bi2d, eecs=input%xs%bse%eecs)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  RR-RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    if(fwp) then 
      call dzmat_send2global_root(localmat, cpmat, bi2d)
      if(mpiglobal%rank == 0) then
        call writecmplxparts('sqrtamb_mat', remat=dble(localmat), immat=aimag(localmat))
      end if
    end if

    ! Construct S Matrix
    ! S = (A-B)^{1/2} (A+B) (A-B)^{1/2}
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse): Constructing S matrix")')
      call timesec(t0)
    end if
    call new_dzmat(auxmat, hamsize, hamsize, bi2d)
    call dzgemm(cpmat, cmat, auxmat)
    call new_dzmat(smat, hamsize, hamsize, bi2d)
    call dzgemm(auxmat, cpmat, smat)
    call del_dzmat(auxmat)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    if(fwp) then 
      call dzmat_send2global_root(localmat, smat, bi2d)
      if(mpiglobal%rank == 0) then
        call writecmplxparts('s_mat', remat=dble(localmat), immat=aimag(localmat))
      end if
    end if

    ! Cmat = (A-B)^1/2
    do j = 1, cmat%ncols_loc
      do i = 1, cmat%nrows_loc
        cmat%za(i,j) = cpmat%za(i,j)
      end do
    end do

    ! Cpmat = (A-B)^-1/2
    if(mpiglobal%rank == 0) then 
      write(unitout, '("Info(setup_dis_ti_bse): Inverting (RR-RA)^1/2 matrix")')
      call timesec(t0)
    end if
    call dzinvert(cpmat)
    if(mpiglobal%rank == 0) then 
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
    end if

    if(mpiglobal%rank == 0) then 
      call timesec(ts1)
      write(unitout, '("Info(setup_dis_ti_bse): Total time needed",f12.3,"s")') ts1-ts0
    end if

  end subroutine setup_dis_ti_bse

  subroutine make_doszstren(doszsr)
    use m_getpmat
    
    implicit none

    ! I/O
    type(dzmat), intent(inout) :: doszsr

    ! Local
    real(8) :: ts0, ts1, t1, t0, sqrteval
    type(dzmat) :: dprojmat, dxpy
    integer(4) :: i,j, lambda, nopt, igqmt

    if(bi2d%isroot) then
      write(unitout, '("Info(b_bse:make_doszstren): Making oszillator strengths (distributed).")')
      write(unitout, '("Info(b_bse:make_doszstren): iqmt=", i4)') iqmt
      write(unitout, '("Info(b_bse:make_doszstren): vqlmt=", 3E11.3)') vqlmt(:,iqmt)
      call timesec(ts0)
    end if

    if(iqmt == 1) then 

      if(bi2d%isroot) then
        write(unitout, '("Info(b_bse:make_doszstren): Using zero momentum transfer formalismus.")')
        write(unitout, '("Info(b_bse:make_doszstren):&
          & Making oszillator strengths using position operator matrix elements.")')
      end if

      nopt=3

      ! Distributed oszillator strengths
      call new_dzmat(doszsr, nexc, nopt, bi2d, rblck=bi2d%mblck, cblck=1)

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Building position operator matrix elements using momentum matrix elements  !
      ! and transition energies. If on top of GW, renormalize the p mat elements.  !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Distributed position operator matrix
      call new_dzmat(dprojmat, hamsize, nopt, bi2d, rblck=bi2d%mblck, cblck=1)

      ! Build position operator matrix elements
      if(mpiglobal%rank == 0) then 
        write(unitout, '("  Building Rmat.")')
        call timesec(t0)
      end if

      ! Build R-matrix from P-matrix
      ! \tilde{R}_{a,i} = 
      !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
      !     P^i_{o_a,u_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a})
      call setup_distributed_rmat(dprojmat, bi2d)
      if(mpiglobal%rank == 0) then
        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      end if
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    else

      if(bi2d%isroot) then
        write(unitout, '("Info(b_bse:make_doszstren):&
          & Making oszillator strengths using plane wave matrix elements.")')
      end if

      nopt=1

      ! Distributed oszillator strengths
      call new_dzmat(doszsr, nexc, nopt, bi2d, rblck=bi2d%mblck, cblck=1)

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Building plane wave matrix elements                                        !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Distributed position operator matrix
      call new_dzmat(dprojmat, hamsize, nopt, bi2d, rblck=bi2d%mblck, cblck=1)

      ! Build plane wave matrix elements
      if(bi2d%isroot) then 
        write(unitout, '("  Building Pwmat.")')
        call timesec(t0)
      end if

      ! Generate plane wave matrix elements for G=Gqmt and q=qmt
      igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
      call setup_distributed_pwmat(dprojmat, iqmt, igqmt, bi2d)
    

      if(bi2d%isroot) then
        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      end if
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    end if

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! TDA case: Build resonant oscillator strengths                              !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if(.not. fcoup) then

      if(bi2d%isroot) then 
        write(unitout, '("  Building distributed resonant oscillator strengths.")')
        call timesec(t0)
      end if

      !! Resonant oscillator strengths
      ! qmt=0 case:
      ! t^R_{\lambda,i} = <X_\lambda|\tilde{R}^{i*}> =
      !   \Sum_{a} X^H_{\lambda, a} \tilde{R}^*_{a,i}
      ! qmt/=0 case:
      ! t^R_{\lambda}(G,qmt) = <X_\lambda|\tilde{M}^*(G,qmt)> =
      !   \Sum_{a} X^H_{\lambda, a} \tilde{M}^*_{a}(G,qmt)
      dprojmat%za=conjg(dprojmat%za)
      call dzgemm(dbevecr, dprojmat, doszsr, transa='C', m=nexc)

      ! Deallocating eigenvectors
      call del_dzmat(dbevecr)
    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! TI case: Build oscillator strength                                         !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if(fcoup .and. fti) then 

      if(mpiglobal%rank == 0) then 
        write(unitout, '("  Building (X+Y) from squared EVP EVs.")')
        call timesec(t0)
      end if

      ! Interested in X^+ + Y^+, so we rescale the 
      ! auxiliary eigenvectors Z by the square root of the eigenvalues
      ! so that 
      ! (X+Y)_{a,lambda} = \Sum_{a'} (A-B)^{1/2}_{a,a'} * E^{-1/2}_lambda * Z_{a',lambda}
      do j = 1, dbevecr%ncols_loc
        lambda = dbevecr%c2g(j)
        sqrteval = sqrt(bevalre(lambda))
        if(lambda <= nexc) then 
          do i = 1, dbevecr%nrows_loc
            dbevecr%za(i, j) = dbevecr%za(i, j) / sqrteval
          end do
        end if
      end do
      call new_dzmat(dxpy, hamsize, nexc, bi2d)
      call dzgemm(dcmat, dbevecr, dxpy, n=nexc)

      ! Deallocating eigenvectors
      call del_dzmat(dbevecr)

      if(mpiglobal%rank == 0) then 
        call timesec(t1)
        write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      end if

      if(mpiglobal%rank == 0) then 
        write(unitout, '("  Building distributed oscillator strengths&
          & for time inverted ar basis.")')
        call timesec(t0)
      end if
      !! Oscillator strengths
      ! qmt=0
      ! t_{\lambda,i} = < (| X_\lambda>+| Y_\lambda>)| \tilde{R}^{i*}> =
      !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{R}^*_{a,i}
      ! qmt/=0
      ! t_{\lambda}(G,qmt) = < (| X_\lambda>+| Y_\lambda>)| \tilde{M}^*(G,qmt)> =
      !     \Sum_{a} (X+Y)^H_{\lambda, a} \tilde{M}^*_{a}
      dprojmat%za=conjg(dprojmat%za)
      call dzgemm(dxpy, dprojmat, doszsr, transa='C', m=nexc) 
      call del_dzmat(dxpy)

    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    ! Projection matrix no longer needed
    call del_dzmat(dprojmat)

    if(bi2d%isroot) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if

    if(bi2d%isroot) then 
      call timesec(ts1)
      write(unitout, '("  Oszillator strengths made in:", f12.3,"s")') ts1-ts0
    end if

  end subroutine make_doszstren
  
  subroutine make_dist_spectrum_ti(nfreq, freq, symsp)
    use mod_lattice, only: omega
    use mod_constants, only: zone, zi, pi
    use modxs, only: symt2
    use invert
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(inout) :: symsp(:,:,:)

    ! Local
    integer(4) :: i, j, o1, o2, ig, jg, nopt
    real(8) :: t1, t0, ts0, ts1
    complex(8) :: brd
    complex(8), allocatable :: buf(:,:,:)
    complex(8), allocatable :: ns_spectr(:,:)
    type(dzmat) :: denw, dtmat, dns_spectr

    if(bi2d%isroot) then
      write(unitout, '("Info(b_bse:make_dist_sp_ti):&
        & Making spectrum using formula for coupling with time inverted ar basis.")')
      if(iqmt == 1) then 
        write(unitout, '("Info(b_bse:makesp_dist_sp_ti):&
         & Using formalism for vaninshing momentum transfer.")')
        if(input%xs%dfoffdiag) then
          write(unitout, '("Info(b_bse:makesp_dist_sp_ti): Including off diagonals.")')
        end if
      else
        write(unitout, '("Info(b_bse:makesp_dist_sp_ti):&
         & Using formalism for finite momentum transfer. iqmt=",i4)') iqmt
      end if
      call timesec(ts0)
    end if

    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make energy denominator for each frequency !
    ! (resonant & anti-resonant)                 !
    !++++++++++++++++++++++++++++++++++++++++++++!
    if(bi2d%isroot) then
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)
    end if

    ! enw_{w, \lambda} = 1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta)
    call new_dzmat(denw, nfreq, nexc, bi2d)

    ! Broadening factor
    brd = zi*input%xs%broad

    !$OMP PARALLEL DO &
    !$OMP& COLLAPSE(2),&
#ifdef SCAL
    !$OMP& DEFAULT(SHARED), PRIVATE(i,j,ig,jg)
#else
    !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
#endif
    do j = 1, denw%ncols_loc
      do i = 1, denw%nrows_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = denw%r2g(i)
        jg = denw%c2g(j)
        denw%za(i,j) = zone/(bevalre(jg)-freq(ig)-brd)&
                     &+ zone/(bevalre(jg)+freq(ig)+brd)
#else
        denw%za(i,j) = zone/(bevalre(j)-freq(i)-brd)&
                     &+ zone/(bevalre(j)+freq(i)+brd)
#endif
      end do
    end do
    !$OMP END PARALLEL DO

    if(bi2d%isroot) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if
    !++++++++++++++++++++++++++++++++++++++++++++!
    
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Helper matrix build form resonant          !
    ! oszillator strenghs                        !
    !++++++++++++++++++++++++++++++++++++++++++++!
    if(bi2d%isroot) then
      write(unitout, '("  Making helper matrix tmat.")')
      call timesec(t0)
    end if

    if(iqmt == 1) then 
      ! tmat_{\lambda, j} = t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! where j combines the Cartesian directions
      if(input%xs%dfoffdiag) then
        nopt = 9
      else
        nopt = 3
      end if
    else
      ! tmat_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
      nopt = 1 
    end if

    call new_dzmat(dtmat, nexc, nopt, bi2d,&
      & rblck=bi2d%mblck, cblck=1)

    do j = 1, dtmat%ncols_loc 
#ifdef SCAL
      ! Get corresponding global indices
      jg = dtmat%c2g(j)
#else
      jg = j
#endif
      ! Get individual opical indices
      if(iqmt == 1) then 
        if(input%xs%dfoffdiag) then
          o2 = (jg-1)/3 + 1
          o1 = jg-(o2-1)*3
        else
          o2 = jg
          o1 = jg
        end if
      else
        o2=1
        o1=1
      end if

      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i,ig)
      do i = 1, dtmat%nrows_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = dtmat%r2g(i)
#else
        ig = i
#endif
        dtmat%za(i,j) = conjg(oszsr(ig,o1))*oszsr(ig,o2)
      end do
      !$OMP END PARALLEL DO
    end do

    if(bi2d%isroot) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if
    !++++++++++++++++++++++++++++++++++++++++++++!
         
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make non-lattice-symmetrized spectrum      !
    !++++++++++++++++++++++++++++++++++++++++++++!
    if(bi2d%isroot) then
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)
    end if

    call new_dzmat(dns_spectr, nfreq, nopt, bi2d,&
      & rblck=bi2d%mblck, cblck=1)

    ! qmt=0 case:
    ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
    !   i.e. nsspectr_{w,j} = 
    !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
    !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
    ! qmt/=0 case:
    ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
    !   i.e. nsspectr_{w,j} = 
    !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
    !       t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
    call dzgemm(denw, dtmat, dns_spectr)
    !++++++++++++++++++++++++++++++++++++++++++++!

    ! Helper no longer needed
    call del_dzmat(dtmat)
    call del_dzmat(denw)

    ! Send spectrum to a global matrix at rank 0 to 
    ! interface with non parallel post-processing routines.
    call dzmat_send2global_root(ns_spectr, dns_spectr, bi2d)
    call del_dzmat(dns_spectr)

    ! Postprocess non-lattice-symmetrized spectrum
    if(bi2d%isroot) then
      call finalizespectrum(ns_spectr, symsp)
      deallocate(ns_spectr)
    end if

    if(bi2d%isroot) then 
      ! Total time for spectrum construction
      call timesec(ts1)
      write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0
    end if

  end subroutine make_dist_spectrum_ti

  subroutine make_dist_spectrum_tda(nfreq, freq, symsp)
    use mod_lattice, only: omega
    use mod_constants, only: zone, zi, pi
    use modxs, only: symt2
    use invert
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(inout) :: symsp(:,:,:)

    ! Local
    integer(4) :: i, j, o1, o2, ig, jg, nopt
    real(8) :: t1, t0, ts0, ts1
    complex(8) :: brd
    complex(8), allocatable :: buf(:,:,:)
    complex(8), allocatable :: ns_spectr(:,:)
    type(dzmat) :: denwr, denwa, dtmatr, dns_spectr

    if(bi2d%isroot) then
      write(unitout, '("Info(b_bse:make_dist_sp_tda):&
        & Making spectrum using TDA formula.")')
      if(iqmt == 1) then 
        write(unitout, '("Info(b_bse:make_dist_sp_tda):&
         & Using formalism for vaninshing momentum transfer.")')
        if(input%xs%dfoffdiag) then
          write(unitout, '("Info(b_bse:make_dist_sp_tda): Including off diagonals.")')
        end if
      else
        write(unitout, '("Info(b_bse:make_dist_sp_tda):&
         & Using formalism for finite momentum transfer. iqmt=",i4)') iqmt
      end if
      if(fti) then 
        write(unitout, '("Info(b_bse:make_dist_sp_tda): Using fomula for ti ar basis.")')
      end if
      call timesec(ts0)
    end if

    if(iqmt /=1 .and. .not. fti) then 
      write(unitout, '("Info(b_bse:make_dist_sp_tda):&
        & Finite momentum transfer only supported for ti ar basis.")')
      call terminate
    end if

    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make energy denominator for each frequency !
    ! (resonant & anti-resonant)                 !
    !++++++++++++++++++++++++++++++++++++++++++++!
    if(bi2d%isroot) then
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)
    end if

    ! Broadening factor
    brd = zi*input%xs%broad


    if(.not. fti) then 

      ! enwr_{w, \lambda} = 1/(E_\lambda - w - i\delta)
      call new_dzmat(denwr, nfreq, nexc, bi2d)
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
#ifdef SCAL
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j,ig,jg)
#else
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
#endif
      do j = 1, denwr%ncols_loc
        do i = 1, denwr%nrows_loc
#ifdef SCAL
          ! Get corresponding global indices
          ig = denwr%r2g(i)
          jg = denwr%c2g(j)
          denwr%za(i,j) = zone/(bevalre(jg)-freq(ig)-brd)
#else
          denwr%za(i,j) = zone/(bevalre(j)-freq(i)-brd)
#endif
        end do
      end do
      !$OMP END PARALLEL DO

      ! enwa_{w, \lambda} = 1/(E_\lambda + w + i\delta)
      call new_dzmat(denwa, nfreq, nexc, bi2d)
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
#ifdef SCAL
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j,ig,jg)
#else
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
#endif
      do j = 1, denwa%ncols_loc
        do i = 1, denwa%nrows_loc
#ifdef SCAL
          ! Get corresponding global indices
          ig = denwa%r2g(i)
          jg = denwa%c2g(j)
          denwa%za(i,j) = zone/(bevalre(jg)+freq(ig)+brd)
#else
          denwa%za(i,j) = zone/(bevalre(j)+freq(i)+brd)
#endif
        end do
      end do
      !$OMP END PARALLEL DO

    else

      ! enwr_{w, \lambda} = 1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta)
      call new_dzmat(denwr, nfreq, nexc, bi2d)
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
#ifdef SCAL
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j,ig,jg)
#else
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
#endif
      do j = 1, denwr%ncols_loc
        do i = 1, denwr%nrows_loc
#ifdef SCAL
          ! Get corresponding global indices
          ig = denwr%r2g(i)
          jg = denwr%c2g(j)
          denwr%za(i,j) = zone/(bevalre(jg)-freq(ig)-brd)&
                      &+ zone/(bevalre(jg)+freq(ig)+brd)
#else
          denwr%za(i,j) = zone/(bevalre(j)-freq(i)-brd)&
                      &+ zone/(bevalre(j)+freq(i)+brd)
#endif
        end do
      end do
      !$OMP END PARALLEL DO

    end if

    if(bi2d%isroot) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if
    !++++++++++++++++++++++++++++++++++++++++++++!
    
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Helper matrix build form resonant          !
    ! oszillator strenghs                        !
    !++++++++++++++++++++++++++++++++++++++++++++!
    if(bi2d%isroot) then
      write(unitout, '("  Making helper matrix tmatr.")')
      call timesec(t0)
    end if

    if(iqmt == 1) then 
      if(input%xs%dfoffdiag) then
        nopt = 9
      else
        nopt = 3
      end if
    else
      nopt = 1
    end if

    ! qmt=0
    ! tmatr_{\lambda, j} = t^{R^*}_{\lambda, o1_j} t^{R}_{\lambda, o2_j} 
    ! where j combines the 2 cartesian directions
    ! qmt/=0
    ! tmatr_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
    call new_dzmat(dtmatr, nexc, nopt, bi2d,&
      & rblck=bi2d%mblck, cblck=1)

    do j = 1, dtmatr%ncols_loc 
#ifdef SCAL
      ! Get corresponding global indices
      jg = dtmatr%c2g(j)
#else
      jg = j
#endif
      ! Get individual opical indices
      if(iqmt==1) then 
        if(input%xs%dfoffdiag) then
          o2 = (jg-1)/3 + 1
          o1 = jg-(o2-1)*3
        else
          o2 = jg
          o1 = jg
        end if
      else
        o2=1
        o1=1
      end if
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(i,ig)
      do i = 1, dtmatr%nrows_loc
#ifdef SCAL
        ! Get corresponding global indices
        ig = dtmatr%r2g(i)
#else
        ig = i
#endif
        dtmatr%za(i,j) = conjg(oszsr(ig,o1))*oszsr(ig,o2)
      end do
      !$OMP END PARALLEL DO
    end do

    if(bi2d%isroot) then
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end if
    !++++++++++++++++++++++++++++++++++++++++++++!
         
    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Make non-lattice-symmetrized spectrum      !
    !++++++++++++++++++++++++++++++++++++++++++++!
    if(bi2d%isroot) then
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)
    end if

    call new_dzmat(dns_spectr, nfreq, nopt, bi2d,&
      & rblck=bi2d%mblck, cblck=1)

    if(.not. fti) then 
      ! nsspectr_{w,j} = \Sum_{\lambda} enwr_{w,\lambda} tmatr_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} 1/(E_\lambda - w - i\delta)
      !       t^{R*}_{\lambda, o1_j} t^R_{\lambda, o2_j} 
      call dzgemm(denwr, dtmatr, dns_spectr)

      ! nsspectr_{w,j} += \Sum_{\lambda} enwa_{w,\lambda} tmatr^*_{\lambda, j}
      !   i.e. nsspectr_{w,j} = nsspectr_{w,j}
      !     \Sum_{\lambda} 1/(E_\lambda + w + i\delta)
      !       (t^{R*}_{\lambda, o1_j} t^R_{\lambda, o2_j})^*
      dtmatr%za = conjg(dtmatr%za)
      call dzgemm(denwa, dtmatr, dns_spectr, beta=zone)
    else
      ! qmt=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! qmt/=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
      call dzgemm(denwr, dtmatr, dns_spectr)
    end if
    !++++++++++++++++++++++++++++++++++++++++++++!

    ! Helper no longer needed
    call del_dzmat(dtmatr)
    call del_dzmat(denwr)
    call del_dzmat(denwa)

    ! Send spectrum to a global matrix at rank 0 to 
    ! interface with non parallel post-processing routines.
    call dzmat_send2global_root(ns_spectr, dns_spectr, bi2d)
    call del_dzmat(dns_spectr)

    ! Postprocess non-lattice-symmetrized spectrum
    if(bi2d%isroot) then
      call finalizespectrum(ns_spectr, symsp)
      deallocate(ns_spectr)
    end if

    if(bi2d%isroot) then 
      ! Total time for spectrum construction
      call timesec(ts1)
      write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0
    end if

  end subroutine make_dist_spectrum_tda

end subroutine b_bse
!EOC
