! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bse
! !INTERFACE:
subroutine b_bse
! !USES:
  use modinput, only: input
  use mod_constants, only: zzero
  use mod_kpoint, only: nkptnr, vkl
  use mod_eigenvalue_occupancy, only: evalsv, nstsv
  use modmpi, only: rank, barrier
  use modxs, only: unitout, bcbs, bsed
  use modbse
  use m_genwgrid
  use m_genfilname
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
!EOP
!BOC

  implicit none

  ! Local variables
  ! Parameters
  integer(4), parameter :: iqmt = 0
  ! Variables
  integer(4) :: iknr, iq, nw
  integer(4) :: hamsize, nexc
  integer(4) :: no, nu, nou, nkkp
  real(8) :: egap, ts0, ts1

  ! Allocatable arrays
  integer(4), allocatable, dimension(:) :: kousize
  integer(4), allocatable, dimension(:,:) :: smap
  logical, allocatable, dimension(:,:) :: kouflag
  real(8), allocatable, dimension(:) :: ofac, beval, w
  complex(8), allocatable, dimension(:,:) :: ham, bevec, oszs
  complex(8), allocatable, dimension(:,:,:) :: symspectr
  ! GW 
  real(8), allocatable, dimension(:,:) :: eval0

  ! Routine not yet parallelized
  mpirank: if(rank .eq. 0) then

    write(*,*) "Hi, this is b_bse."

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
    egap = minval(evalsv(bcouabs%il2,1:nkptnr) - evalsv(bcouabs%iu1,1:nkptnr) + input%xs%scissor)
    write(unitout,*)
    write(unitout, '("Info(bse): gap:", E23.16)') egap
    ! Warn if the system has no gap even with scissor (or no scissor and on top of GW)
    if(egap .lt. input%groundstate%epspot) then
      write(unitout,*)
      write(unitout, '("Warning(bse): the system has no gap, setting it to 0")')
      write(unitout,*)
      egap = 0.0d0
    end if  

write(*,*) "About to check occupancies"
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
    call checkoccupancies(iq, nkptnr, bcouabs, hamsize, smap, kouflag, kousize, ofac)
    !!-->

    ! Allocate BSE-Hamiltonian (large matrix, up to several Gb)
    allocate(ham(hamsize, hamsize))
    ham(:, :) = zzero

    ! Number of kkp combinations with ikp >= ik
    nkkp = nkptnr*(nkptnr+1)/2
    ! Number of o u combinations
    no  = bcou%n1
    nu  = bcou%n2
    nou = no * nu

    ! Write Info
    write(unitout,*)
    write(unitout, '("Info(bse): Assembling BSE matrix")')
    write(unitout, '(" BSE-Hamiltonian: Shape=",i4, " Size=",i8)') hamsize, hamsize**2
    ! Assemble Hamiltonian matrix 
    call timesec(ts0)
    call setuphamiltonian(ham, .true.)
  write(*,*) "After setup"
    call timesec(ts1)
    write(unitout, '(" Timing (in seconds)	   :", f12.3)') ts1 - ts0


    ! Number of excitons to consider
    nexc = input%xs%bse%nexc
    if(nexc > hamsize .or. nexc < 1) then
      nexc = hamsize
    end if

    ! Write Info
    write(unitout,*)
    write(unitout, '("Info(bse): Invoking lapack routine ZHEEVR")')
    write(unitout, '(" Number of requested solutions : ", i8)') nexc

    ! Allocate eigenvector and eigenvalue arrays
    allocate(beval(hamsize), bevec(hamsize, nexc))

  write(*,*) "Before diag"
    ! Diagonalize Hamiltonian (destroys the content of ham)
    call timesec(ts0)
    call bsesoldiag(nexc, hamsize, ham, beval, bevec)
    call timesec(ts1)
  write(*,*) "after diag"

    ! Deallocate BSE-Hamiltonian
    deallocate(ham)

  write(*,*) "after ham deallocation"
    write(unitout, '(" Timing (in seconds)	   :", f12.3)') ts1 - ts0
    write(unitout,*)

    ! Calculate oscillator strengths.
    allocate(oszs(nexc,3))

    call makeoszillatorstrength(oszs)

    ! Write excition energies and oscillator strengths to 
    ! text file. 
    call writeoscillator(hamsize, nexc, egap, beval, oszs)

    ! Calculate macroscopic dielectric tensor with finite broadening, i.e. "the spectrum"
    nw = input%xs%energywindow%points

    ! Allocate arrays used in spectrum construction
    allocate(w(nw))
    allocate(symspectr(3,3,nw))

    ! Generate an evenly spaced frequency grid 
    call genwgrid(nw, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)

    ! Calculate lattice symmetrized spectrum.
    call makespectrum(nw, w, symspectr)

    ! Generate and write derived optical quantities
    iq = 1
    call writederived(iq, symspectr, nw, w)

    ! Store excitonic energies and wave functions to file
    if(associated(input%xs%storeexcitons)) then
      call storeexcitons(hamsize,nexc, nkptnr, iuref, bcou,smap,beval,bevec)
    end if

    ! Clean up
    deallocate(beval,bevec,oszs,w,symspectr, evalsv)
    if(associated(input%gw)) deallocate(eval0)

    call barrier

  else mpirank

    call barrier

  end if mpirank

contains

  subroutine setuphamiltonian(ham, writeparts)

    complex(8),  intent(inout) :: ham(hamsize, hamsize)
    logical, intent(in), optional :: writeparts

    complex(8), dimension(no,nu,no,nu) :: excli, sccli
    integer(4) :: ikkp
    integer(4) :: ik1, ik2
    integer(4) :: s1l, s1u, s2l, s2u
    logical :: wp

    real(8), allocatable :: wre(:,:), wim(:,:), vre(:,:), vim(:,:), kstransen(:,:)

    if(present(writeparts)) then
      wp = writeparts
    else
      wp = .false.
    end if
    if(wp) then
      allocate(wre(hamsize,hamsize),wim(hamsize,hamsize))
      allocate(vre(hamsize,hamsize),vim(hamsize,hamsize))
      allocate(kstransen(hamsize,hamsize))
      wre=0.0d0
      wim=0.0d0
      vre=0.0d0
      vim=0.0d0
      kstransen=0.0d0
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

      ! Read corresponding blocks of W and V from file
      select case(trim(input%xs%bse%bsetype))
        case('singlet', 'triplet')
          ! Read screened coulomb interaction W_{ouki,o'u'kj}
          call getbsemat('SCCLI.OUT', ikkp, no, nu, sccli)
      end select
      select case(trim(input%xs%bse%bsetype))
        case('RPA', 'singlet')
          ! Read exchange interaction v_{ouki,o'u'kj}
          call getbsemat('EXCLI.OUT', ikkp, no, nu, excli)
      end select

      ! Get ik1 and ik2 from ikkp
      call kkpmap(ikkp, nkptnr, ik1, ik2)

      ! Calculate where to write the BSE matrix
      s1l = sum(kousize(1:ik1-1))+1
      s1u = s1l + kousize(ik1) - 1 
      s2l = sum(kousize(1:ik2-1))+1
      s2u = s2l + kousize(ik2) - 1

!write(*,*) "Building for kkp=", ikkp, ":", ik1, ik2
!write(*,*) "Block", s1l, s1u, ":", s2l, s2u

      ! For blocks on the diagonal, add the KS transition
      ! energies to the diagonal of the block.
      if(ik1 .eq. ik2) then
        if(wp) then
          call writekstrans(ik1, ham(s1l:s1u,s2l:s2u), remat=kstransen(s1l:s1u,s2l:s2u))
        else
          call writekstrans(ik1, ham(s1l:s1u,s2l:s2u))
        end if
      end if
        
      ! Add exchange term
      ! + 2* v_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
      ! * sqrt(abs(f_{o_{s1} k_{s1}} - f_{u_{s1} k_{s1}}))
      ! * sqrt(abs(f_{o_{s2} k_{s2}} - f_{u_{s2} k_{s2}}))
      select case(trim(input%xs%bse%bsetype))
        case('RPA', 'singlet')
          if(wp) then
            call addexchange(ik1, ik2,&
              & ham(s1l:s1u,s2l:s2u), excli,&
              & ofac(s1l:s1u), ofac(s2l:s2u),&
              & remat=vre(s1l:s1u,s2l:s2u), immat=vim(s1l:s1u,s2l:s2u))
          else
            call addexchange(ik1, ik2,&
              & ham(s1l:s1u,s2l:s2u), excli, ofac(s1l:s1u), ofac(s2l:s2u))
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
              & ham(s1l:s1u,s2l:s2u), sccli,&
              & ofac(s1l:s1u), ofac(s2l:s2u),&
              & remat=wre(s1l:s1u,s2l:s2u), immat=wim(s1l:s1u,s2l:s2u))
          else
            call adddirect(ik1, ik2,&
              & ham(s1l:s1u,s2l:s2u), sccli, ofac(s1l:s1u), ofac(s2l:s2u))
          end if
      end select

    end do

    if(wp) then
      call writecmplxparts('W',wre,immat=wim)
      call writecmplxparts('V',vre,immat=vim)
      call writecmplxparts('E',kstransen)
    end if


    if(wp) then
      deallocate(wre,wim)
      deallocate(vre,vim)
      deallocate(kstransen)
    end if

  end subroutine setuphamiltonian

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
        ! Subtract gap if present, so that the exciton
        ! energies will be binding energies.
        tmp = devals(iou) - egap + bsed
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

    !! If excli would be in shape excli(nou, nou) the 
    !! following could be done more elegantly.
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


    !! If sccli would be in shape sccli(nou, nou) the 
    !! following could be done more elegantly.
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
  

  subroutine makeoszillatorstrength(oszstr)
    use mod_constants, only: zone
    use m_getpmat
    
    implicit none

    ! I/O
    complex(8), intent(out) :: oszstr(nexc,3)

    ! Local
    complex(8), allocatable :: pmou(:,:,:), rmat(:,:)
    integer(4) :: io, iu, iuabs, ik, ikprev
    integer(4) :: a1
    
    ! Allocate momentum matrix slice needed.
    allocate(pmou(3, no, nu))
    allocate(rmat(hamsize, 3))

    ! Loop over 3 directions
    ikprev = -1

    do a1 = 1, hamsize
      
      ! Get indices
      ik = smap(a1,1)
      io = smap(a1,2)
      iu = smap(a1,3)        ! Relative 
      iuabs = iu + iuref - 1  ! Absolute

      ! Read momentum matrix slice for given k-point
      ! if not already present (k index varies the slowest in the smap)
      if(ik /= ikprev) then 
        call getpmat(ik, vkl,& 
          & bcouabs%il1, bcouabs%iu1,&
          & bcouabs%il2, bcouabs%iu2,&
          & .true., 'PMAT_XS.OUT', pmou)
      end if
      ikprev = ik

!! Discuss
      ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
      ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
      !   Where E are quasi-particle energies and e are KS energies.
      if(associated(input%gw)) then
         pmou(1:3,io,iu) = pmou(1:3,io,iu)&
           &* (evalsv(iuabs,ik)-evalsv(io,ik))/(eval0(iuabs,ik)- eval0(io,ik))
      end if 


      ! Build complex conjugate R-matrix from p-matrix
      ! \tilde{R}^*_{u_{s1},o_{s1},k_{s1}},i = 
      ! (f_{o_{s1},k_{s1}}-f_{u_{s1},k_{s1}}) P_{o_{s1},u_{s1},k_{s1}},i /(e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}})
      rmat(a1, :) = ofac(a1) * pmou(:, io, iu)/(evalsv(iuabs, ik) - evalsv(io, ik))

    end do

    ! Momentum matrix elements no longer needed
    deallocate(pmou)

    oszstr = zzero

    ! t_\lambda,i = < \tilde{R}^i | \lambda> =
    ! \Sum_{s1} f_{s1} R^*_{{u_{s1} o_{s1} k_{s1}},i} A_{o_{s1} u_{s1} k_{s1}},\lambda= 
    ! \Sum_{s1} f_{s1} P^*_{{u_{s1} o_{s1} k_{s1}},i}/(e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) A_{o_{s1} u_{s1} k_{s1}},\lambda= 
    ! \Sum_{s1} f_{s1} P_{{o_{s1} u_{s1} k_{s1}},i}/(e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}}) A_{o_{s1} u_{s1} k_{s1}},\lambda
    call zgemm('t','n', nexc, 3, hamsize, zone, bevec, hamsize, rmat, hamsize, zzero, oszstr, nexc)
    deallocate(rmat)

  end subroutine makeoszillatorstrength


  subroutine makespectrum(nfreq, freq, spectrum)
    use mod_lattice, only: omega
    use mod_constants, only: zi, pi
    use modxs, only: symt2
    implicit none

    ! I/O
    integer(4), intent(in) :: nfreq
    real(8), intent(in) :: freq(nfreq)
    complex(8), intent(out) :: spectrum(3,3,nfreq)

    ! Local
    integer(4) :: o1, o2, ol, ou, iw
    integer(4) :: lambda
    integer(4) :: optcompt(3)
    real(8) :: eextot
    complex(8) :: buf(3,3,nw)
    complex(8) :: spectr(nw)

    ! External functions
    integer(4), external :: l2int

    ! Make non-lattice-symmetrized spectrum
    buf = zzero 

    do o1 = 1, 3

      ! Stk: compute off-diagonal optical components if requested
      if(input%xs%dfoffdiag) then
        ol = 1
        ou = 3
      else
        ol = o1
        ou = o1
      end if

      do o2 = ol, ou

        optcompt = (/ o1, o2, 0 /)
        spectr(:) = zzero

        do iw = 1, nfreq
          do lambda = 1, nexc
            
            ! Get total energy of excition, i.e. add gap energy
            eextot = beval(lambda) + egap
!! Discuss what is bsed and does anyone use it
            ! Substract bsed part (defautl is is zero)
            eextot = eextot - bsed

            ! Lorentzian line shape
            ! Resonant part
            spectr(iw) = spectr(iw) + oszs(lambda, o1) * conjg(oszs(lambda, o2))&
              &* (1.d0/(eextot-freq(iw)-zi*input%xs%broad))             
            ! Anti-resonant contribution (Default is true)
!! Discuss the conjugation was the other way around originally
            if(input%xs%bse%aresbse) spectr(iw) = spectr(iw)&
              &+ conjg(oszs(lambda, o1)) * oszs(lambda, o2)&
              &* (1.d0/(eextot+freq(iw)+zi*input%xs%broad))

          end do
        end do

        spectr(:) = l2int(o1 .eq. o2) * 1.d0 + spectr(:) * 8.d0*pi/omega/nkptnr

        buf(o1,o2,:)=spectr(:)

      end do
    end do


    ! Symmetrize spectrum 
    spectrum = zzero
    do o1=1,3
      do o2=1,3

        ! Symmetrize the macroscopic dielectric tensor
        call symt2app(o1, o2, nfreq, symt2, buf, spectrum(o1,o2,:))

      end do 
    end do
  end subroutine makespectrum

  subroutine writecmplxparts(fbasename, remat, immat, ik1, ik2)
    use m_getunit
    character(*), intent(in) :: fbasename
    real(8), intent(in) :: remat(:,:)
    real(8), intent(in), optional :: immat(:,:)
    integer(4), intent(in), optional :: ik1, ik2

    integer(4) :: un, a1, a2, n, m
    character(256) :: fname, tmp1, tmp2, tmp3

    n = size(remat,1)
    m = size(remat,2)

    tmp1 = ''
    if(present(ik1)) then
      write(tmp1, '(I8)') ik1
    end if
    tmp2 = ''
    if(present(ik2)) then
      write(tmp2, '(I8)') ik2
    end if
    tmp3= trim(adjustl(tmp1))//trim(adjustl(tmp2))
    if(present(ik1) .or. present(ik2)) then
      tmp3 = "_"//trim(adjustl(tmp3))
    end if

    fname =''
    write(fname,'("Re_",a,a,".OUT")') trim(adjustl(fbasename)),trim(adjustl(tmp3))

    call getunit(un)
    open(unit=un, file=fname, action='write', status='replace')
    do a1=1,n
      write(un,'(SP,E23.16)', advance='no') remat(a1,1)
      do a2=2,m
        write(un, '(SP,1x,E23.16)', advance='no') remat(a1,a2)
      end do
      write(un,*)
    end do
    close(un)

    if(present(immat)) then
      fname =''
      write(fname,'("Im_",a,a,".OUT")') trim(adjustl(fbasename)),trim(adjustl(tmp3))

      call getunit(un)
      open(unit=un, file=fname, action='write', status='replace')
      do a1=1,n
        write(un,'(SP,E23.16)', advance='no') immat(a1,1)
        do a2=2,m
          write(un, '(SP,1x,E23.16)', advance='no') immat(a1,a2)
        end do
        write(un,*)
      end do
      close(un)

    end if

  end subroutine writecmplxparts
    
end subroutine b_bse
!EOC
