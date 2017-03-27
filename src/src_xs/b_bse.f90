! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_bse
! !INTERFACE:
subroutine b_bse(iqmt)
! !USES:
  ! Basics
  use modinput, only: input
  use mod_constants, only: zone, zi, zzero, pi, h2ev
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
  use m_setup_bse
  use m_genexevec
  use m_putgetexcitons
  use m_makeoscistr
  use m_makespectrum

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
  integer(4) :: i
  integer(4) :: ik, ikq, io, iu, a1
  integer(4) :: nexc
  real(8) :: ts0, ts1
  logical :: fcoup, fwp, fscal, fti, fip

  ! Allocatable arrays
  real(8), allocatable, dimension(:) :: bevalim, bevalre, w
  complex(8), allocatable, dimension(:,:) :: ham, bevecr
  complex(8), allocatable, dimension(:,:) :: resvec, aresvec 
  complex(8), allocatable, dimension(:,:) :: oscsr, oscsa
  complex(8), allocatable, dimension(:,:) :: cmat, cpmat
  complex(8), allocatable, dimension(:,:,:) :: symspectr

  real(8) :: bsegap
  real(8) :: v1, v2, en1, en2
  integer(4) :: i1, i2, iex1, iex2, nreq
  integer(4) :: nsymcrys_save
  logical :: efind
  ! Distributed arrays
  type(dzmat) :: dham, dbevecr, doscsr, dresvec, daresvec
  type(dzmat) :: dcmat, dcpmat

  ! Indepenedent particle formalism is used
  if(input%xs%bse%bsetype == "IP") then
    fip = .true.
  else
    fip = .false.
  end if

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
    call timesec(ts0)
    call init0
    call timesec(ts1)
    write(unitout, '("Info(b_bse):&
      & Init0 time:", f12.6)') ts1 - ts0
    ! k-grid init
    call timesec(ts0)
    call init1
    call timesec(ts1)
    write(unitout, '("Info(b_bse):&
      & Init1 time:", f12.6)') ts1 - ts0
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
    call timesec(ts0)
    call init2
    call timesec(ts1)
    write(unitout, '("Info(b_bse):&
      & Init2 time:", f12.6)') ts1 - ts0

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
      if(allocated(eval0)) deallocate(eval0)
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
    if(fip) then 
      write(unitout, '("Info(b_bse): IP requested, nothing much to do")') 
    else
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
    end if

    ! Assemble Hamiltonian matrix 
    if(.not. fip) then 

      if(fcoup) then 
        if(fti) then 
          call setup_bse_ti(ham, cmat, cpmat,  iqmt)
          if(fwp) then 
            call writecmplxparts('HamS', dble(ham), immat=aimag(ham))
            call writecmplxparts('cmat', dble(cmat), immat=aimag(cmat))
            call writecmplxparts('cpmat', dble(cpmat), immat=aimag(cpmat))
          end if
        else
          allocate(ham(2*hamsize,2*hamsize))
          call setup_bse_full(ham, iqmt)
          if(fwp) then 
            call writecmplxparts('Global_Ham', dble(ham), immat=aimag(ham))
          end if
        end if
      else
        allocate(ham(hamsize,hamsize))
        call setup_bse_block(ham, iqmt, .false., .false.)
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
      if(fcoup .and. .not. fti) then 
        allocate(bevalre(2*hamsize))
        allocate(bevalim(2*hamsize))
        allocate(bevecr(2*hamsize, 2*hamsize))
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
          call writecmplxparts('bevecaux', dble(bevecr), immat=aimag(bevecr))
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
          call genexevec(iex1, iex2, nexc, cmat, cpmat, bevecr, bevalre,&
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

    ! IP
    else

      if(fcoup .and. .not. fti) then
        nexc = 2*hamsize 
        allocate(bevalre(nexc))
        allocate(bevalim(nexc))
        bevalim = 0.0d0
        ! Eigenvalues are always saved in ascending order,
        ! doing the same for IP
        do i = 1, hamsize
          bevalre(i) = -de(ensortidx(hamsize-i+1))
        end do
        bevalre(hamsize+1:2*hamsize) = de(ensortidx) 
      else 
        nexc = hamsize
        allocate(bevalre(nexc))
        bevalre(1:hamsize) = de(ensortidx) 
      end if

    end if

    ! Calculate oscillator strengths.
    allocate(oscsr(nexc,3))
    if(fcoup .and. .not. fti) then
      allocate(oscsa(nexc,3))
    end if

    if(fcoup .and. .not. fti) then 
      call makeoscistr(iqmt, nexc, bevecr, oscsr, oscstra=oscsa)
    else
      if(fti) then
        call makeoscistr(iqmt, nexc, bevecr, oscsr, bevalre=bevalre, cmat=cmat)
      else
        call makeoscistr(iqmt, nexc, bevecr, oscsr)
      end if
    end if

    ! Write excition energies and oscillator strengths to 
    ! text file. 
    write(unitout, '("Info(b_bse):&
      & Writing excition energies and oscillator strengths to text file.")')
    if(fcoup .and. .not. fti) then
      call writeoscillator(2*hamsize, nexc, bsegap+sci, bevalre, oscsr,&
        & evalim=bevalim, oscstra=oscsa, sort=.true., iqmt=iqmt)
    else
      call writeoscillator(hamsize, nexc, -(bsegap+sci), bevalre, oscsr, iqmt=iqmt)
    end if

    ! Calculate lattice symmetrized spectrum.
    if(fcoup) then 
      if(fti) then
        call makespectrum_ti(iqmt, nexc, nk_bse, bevalre, oscsr, symspectr)
      else
        call makespectrum_full(iqmt, nexc, nk_bse, bevalre, oscsr, oscsa, symspectr)
      end if
    else
      call makespectrum_tda(iqmt, nexc, nk_bse, bevalre, oscsr, symspectr)
    end if

    ! Generate an evenly spaced frequency grid 
    allocate(w(nw))
    call genwgrid(nw, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)

    ! Generate and write derived optical quantities
    call writederived(iqmt, symspectr, nw, w)

    ! Clean up
    deallocate(bevalre, oscsr, w, symspectr, evalsv)
    if(allocated(bevecr)) deallocate(bevecr)
    if(fcoup) then 
      if(fti) then 
        if(allocated(cmat)) deallocate(cmat)
        if(allocated(cpmat)) deallocate(cpmat)
      else
        if(allocated(bevalim)) deallocate(bevalim)
        if(allocated(oscsa)) deallocate(oscsa)
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
    if(bi2d%isroot) then 
      call timesec(ts0)
    end if
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
    if(bi2d%isroot) then 
      call timesec(ts1)
      write(unitout, '("Info(b_bse):&
        & Init time:", f12.6)') ts1 - ts0
    end if

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
      if(allocated(eval0)) deallocate(eval0)
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
      !if(bi2d%isroot) then
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
      !end if
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

      if(.not. fip) then 

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
          call setup_bse_ti_dist(dham, dcmat, dcpmat, iqmt, bi2d)
        else
          call setup_bse_block_dist(dham, iqmt, .false., .false., bi2d)
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

        ! End root writeout
        end if

      ! IP
      else

        if(bi2d%isroot) then 
          write(unitout, '("Info(b_bse): IP requested, nothing much to do")') 
        end if

        nexc = hamsize
        allocate(bevalre(nexc))
        bevalre(1:hamsize) = de(ensortidx) 

      end if

    ! End on 2d process grid
    end if

    ! Calculate oscillator strengths.
    ! Note: Deallocates eigenvectors
    if(fti) then 
      call makeoscistr_dist(iqmt, nexc, dbevecr, doscsr, bi2d, bevalre, dcmat)
    else
      call makeoscistr_dist(iqmt, nexc, dbevecr, doscsr, bi2d)
    end if

    if(bi2d%isactive) then 

      ! Every process gets a copy of the oscillator strength
      ! (actually only rank 0 writes them to file, but is is not much 
      !  memory and it make the setup for the spectrum calculation easier) 
      call dzmat_send2global_all(oscsr, doscsr, bi2d)
      if(fip) then
        oscsr = oscsr(ensortidx,:)
      end if

      if(bi2d%isroot) then
        ! Write excition energies and oscillator strengths to 
        ! text file. 
        write(unitout, '("Info(b_bse):&
          & Writing excition energies and oscillator strengths to text file.")')
        call writeoscillator(hamsize, nexc, -(bsegap+sci), bevalre, oscsr, iqmt=iqmt)
      end if

      ! Allocate arrays used in spectrum construction
      if(bi2d%isroot) then
        allocate(symspectr(3,3,nw))
      end if

      ! Only process 0 gets an acctual output for symspectr
      if(fcoup .and. fti) then 
        call makespectrum_ti_dist(iqmt, nexc, nk_bse, bevalre, oscsr, symspectr, bi2d)
      else
        call makespectrum_tda_dist(iqmt, nexc, nk_bse, bevalre, oscsr, symspectr, bi2d)
      end if

      if(bi2d%isroot) then

        ! Allocate arrays used in spectrum construction
        allocate(w(nw))
        ! Generate an evenly spaced frequency grid 
        call genwgrid(nw, input%xs%energywindow%intv,&
          & input%xs%tddft%acont, 0.d0, w_real=w)

        ! Generate and write derived optical quantities
        call writederived(iqmt, symspectr, nw, w)

      end if

      ! Clean up
      deallocate(bevalre, oscsr, evalsv)
      if(associated(input%gw)) deallocate(eval0)
      if(bi2d%isroot) then 
        deallocate(w)
        deallocate(symspectr)
      end if

      ! Exit BLACS 
      call exitblacs(bi2d)
      call exitblacs(bi1d)
      call exitblacs(bi0d)

      ! MPI
      ! Ranks that are on the BLACS grid signal that they are done
      call barrier

    ! MPI rank not on 2D grid
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

end subroutine b_bse
!EOC
