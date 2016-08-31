! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bse
! !INTERFACE:
subroutine bse
! !USES:
  use modinput, only: input
  use mod_lattice, only: omega
  use mod_constants, only: zzero, zi, pi
  use mod_kpoint, only: nkptnr, vkl, ivknr
  use mod_qpoint, only: iqmap
  use mod_eigenvalue_occupancy, only: evalsv, nstsv
  use modmpi, only: rank, barrier
  use modxs, only: sta1, sta2, sto1, sto2,&
                 & istocc0, istocc, istunocc0, istunocc,&
                 & isto0, isto, istu0, istu,&
                 & unitout,&
                 & istl1, istl2, istl3, istl4,&
                 & istu1, istu2, istu3, istu4,&
                 & nst1, nst2, nst3, nst4,&
                 & bsed, iqmapr, dgrid, iksubpt,&
                 & escale, fneps, fnloss, fnsigma,&
                 & fnsumrules, symt2
  use m_genwgrid
  use m_getpmat
  use m_genfilname
  use m_getunit
  use m_genloss
  use m_gensigma
  use m_gensumrls
  use m_writeeps
  use m_writeloss
  use m_writesigma
  use m_writesumrls
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
  integer(4), parameter :: noptcmp = 3
  real(8), parameter :: epsortho = 1.d-12
  ! Variables
  character(256) :: fnexc, fnexcs, dotext
  integer(4) :: stat, ievec
  integer(4) :: iknr, jknr, iqr, iq, iw, iv2(3)
  integer(4) :: s1, s2, hamsiz, nexc, ne
  integer(4) :: unexc, ist1, ist2, ist3, ist4
  integer(4) :: ikkp, iv, ic, optcompt(3)
  integer(4) :: oct1, oct2, octu, octl
  integer(4) :: nrnst1, nrnst2, nrnst3, nrnst4
  real(8) :: de, egap, ts0, ts1, sumrls(3)
  ! Allocatable arrays
  integer(4), allocatable, dimension(:) :: sor
  real(8), allocatable, dimension(:) :: beval, w, oszsa, kdocc
  real(8), allocatable, dimension(:,:) :: docc, eval0
  real(8), allocatable, dimension(:,:,:) :: loss
  complex(8), allocatable, dimension(:) :: spectr, sigma
  complex(8), allocatable, dimension(:,:) :: ham, bevec, pmat, oszs
  complex(8), allocatable, dimension(:,:,:) :: pm, buf
  complex(8), allocatable, dimension(:,:,:,:) :: excli, sccli

  ! External functions
  integer(4), external :: l2int

  ! Routine not yet parallelized
  mpirank: if(rank .eq. 0) then

    call init0
    call init1
    call init2
    call xssave0

    ! Read Fermi energy from file
    call readfermi

    ! Initialize the selected ranges of valence/core and conduction states(note
    ! that core states are always understood to be local orbitals set by the user
    ! and actually treated as valence states)

    ! Lowest occupied state (absolute index)
    sta1 = input%xs%bse%nstlbse(1)
    ! Highest occupied state (absolute index)
    sto1 = input%xs%bse%nstlbse(2)
    ! Lowest unoccupied state (counted from first unoccupied state)
    sta2 = input%xs%bse%nstlbse(3)
    ! Highest unoccupied state (counted from first unoccupied state)
    sto2 = input%xs%bse%nstlbse(4)

    ! Number of occupied states
    nrnst1 = sto1-sta1+1
    nrnst2 = sto1-sta1+1
    ! Number of unoccupied states
    nrnst3 = sto2-sta2+1
    nrnst4 = sto2-sta2+1 

    ! Use eigenvector files from screening-calculation
    call genfilname(dotext='_SCR.OUT', setfilext=.true.)

    ! Find occupation limits for k and k+q, where q=0 ist set fixed by iqmt
    call findocclims(iqmt, istocc0, istocc, istunocc0, istunocc, isto0, isto, istu0, istu)

    ! emattype 2 selects o-o and u-u combinations
    input%xs%emattype = 2
    ! Set nstX, istlX, istuX variables for X: 1=o, 2=o, 3=u, 4=u
    call ematbdcmbs(input%xs%emattype)
    
    ! Write out state ranges
    write(unitout,*)
    write(unitout, '("Info(bse): information on number of states:")')
    write(unitout, '(" number of states below Fermi energy in Hamiltonian:", i6)') nrnst1
    write(unitout, '(" number of states above Fermi energy in Hamiltonian:", i6)') nrnst3
    write(unitout, '(" ranges of states according to BSE matrix:")')
    write(unitout, '("  range of first index and number  :", 2i6, 3x, i6)')&
      & istl1, istu1, nst1
    write(unitout, '("  range of second index and number :", 2i6, 3x, i6)')&
      & istl2, istu2, nst2
    write(unitout, '("  range of third index and number  :", 2i6, 3x, i6)')&
      & istl3, istu3, nst3
    write(unitout, '("  range of fourth index and number :", 2i6, 3x, i6)')&
      & istl4, istu4, nst4

    ! Size of BSE-Hamiltonian = #o * #u * #k
    hamsiz = nrnst1 * nrnst3 * nkptnr

    ! Allocations
    ! Allocate arrays for coulomb interactions
    ! W_{ouk,o'u'k}(q=0) 
    allocate(sccli(nrnst1, nrnst3, nrnst2, nrnst4)) ! #o,#u,#o,#u
    ! v_{ouk,o'u'k'}(q=0)
    allocate(excli(nrnst1, nrnst3, nrnst2, nrnst4)) ! #o,#u,#o,#u

    ! Allocate array for occupation number difference (future use)
    allocate(docc(nrnst1, nrnst3)) ! #o,#u
    allocate(kdocc(hamsiz))

    ! Allocate BSE-Hamiltonian (large matrix, up to several Gb)
    allocate(ham(hamsiz, hamsiz))
    ham(:, :) = zzero

    ! Read KS energies
    do iknr = 1, nkptnr
      call getevalsv(vkl(1, iknr), evalsv(1, iknr))
    end do

    ! If on top of GW
    if(associated(input%gw)) then
      ! To KS eigenvalues to use them later for renormalizing PMAT
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

      if(input%xs%bse%bsedirsing) then
        call getbsediag
        write(unitout, '("Info(bse): read diagonal of BSE kernel")')
        write(unitout, '(" mean value : ", 2g18.10)') bsed
      end if

    end if

    ! Determine the minimal optical gap
    egap = 1.d8
    do iknr = 1, nkptnr
      ! Occupied
      do ist1 = sta1, sto1
        ! Unoccupied
        do ist3 = sta2, sto2
          egap = min(egap, evalsv(ist3+istocc, iknr)-evalsv(ist1,iknr)+input%xs%scissor)
        end do
      end do
    end do

    write(unitout, '("Info(bse): gap:", g18.10)') egap

    ! Warn if the system has no gap
    if(egap .lt. input%groundstate%epspot) then
      write(unitout,*)
      write(unitout, '("Warning(bse): the system has no gap")')
      write(unitout,*)
    end if  

    ! Set up BSE-Hamiltonian
    ikkp = 0

    k: do iknr = 1, nkptnr
      kp: do jknr = iknr, nkptnr ! Uses hermiticity
      !! Note: If the Hamilton matrix 
      !! has the elements H_{i,j} and the indices enumerate the
      !! states according to
      !! i = o1u1k1, o1u1k2, ..., o1u1kN, o2u1k1, ..., o2u1kN, ...,
      !!     oMu1kN, oMu2k1, ..., oMuMkN
      !! then because of H_{j,i} = H^*_{i,j} only kj = ki,..,kN is 
      !! needed.

        ikkp = ikkp + 1

        ! Get kj-ki = q for extracting the correct W element
        iv2(:) = ivknr(:, jknr) - ivknr(:, iknr)
        iv2(:) = modulo(iv2(:), input%groundstate%ngridk(:))
        ! Q-point(reduced)
        iqr = iqmapr(iv2(1), iv2(2), iv2(3))
        ! Q-point(non-reduced)
        iq = iqmap(iv2(1), iv2(2), iv2(3))

        ! Read W_{i,j}
        select case(trim(input%xs%bse%bsetype))
          case('singlet', 'triplet')
            ! Read screened coulomb interaction W_{ouki,o'u'kj}
            call getbsemat('SCCLI.OUT', ikkp, nrnst1, nrnst3, sccli)
        end select

        ! Read exchange coulomb interaction
        select case(trim(input%xs%bse%bsetype))
          case('RPA', 'singlet')
            ! Read exchange interaction v_{ouki,o'u'kj}
            call getbsemat('EXCLI.OUT', ikkp, nrnst1, nrnst3, excli)
        end select

        ! Get occupation numbers (future use for partial occupancy) 			
        !   Calculates docc(l,m)= f_{o_l ki} - f_{u_m kj}
        !   Remark: Currently the iq argument does nothing at all
        call getdocc(iq, iknr, jknr, sta1, sto1, istl3, istl3+nrnst3-1, docc)
        
        ! Set up k/kp-part of Hamilton matrix
        ! Occupied
        do ist1 = sta1, sto1
          ! Unoccupied
          do ist3 = sta2, sto2
            ! Occupied
            do ist2 = sta1, sto1
              ! Unoccupied
              do ist4 = sta2, sto2

                ! The combined index labels the states as
                ! s1 = {o1u1k1, o1u2k1, ..., o1uMk1,
                !      o2uMk1, ..., oMuMk1, o1u1k2, ..., oMuMkN} <-> {1,...,M**2N}
                ! s2 does the same.
                s1 = hamidx(ist1-sta1+1, ist3-sta2+1, iknr, nrnst1, nrnst3)
                s2 = hamidx(ist2-sta1+1, ist4-sta2+1, jknr, nrnst2, nrnst4)

!!! Wrong occupation differences?
                ! Get occupation difference for current s1 combination
                ! but for different k's !!!!! those should be at the same k 
                ! at least in this optical case.
                ! kdocc(s1) = f_{o_{s1},k_i} - f_{u_{s1},k_j}
                kdocc(s1) = docc(ist1-sta1+1, ist3-sta2+1)

                ! Add diagonal term
                if(s1 .eq. s2) then
                  ! de = e_{u_{s1}, ki} - e_{o_{s1}, k_i} + scissor
                  de = evalsv(ist3+istl3-sta2, iknr) - evalsv(ist1, iknr)&
                    &+ input%xs%scissor
                  ! BSE-Diag method outputs BINDING energies.
                  !   i.e. subtract gap energy
                  ham(s1, s2) = ham(s1, s2) + de - egap + bsed
                end if

                ! Write partially occupied states in the output    
                if(kdocc(s1) .gt. input%groundstate%epsocc&
                  & .and. kdocc(s1) .lt. 2.d0-input%groundstate%epsocc) then
                  write(*,*) 'kdocc s1', kdocc(s1), s1
                end if

                ! Add exchange term
                ! + 2* v_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
!!! Wrong occupation differences?
                ! Also add occupation difference treatment
                select case(trim(input%xs%bse%bsetype))
                  case('RPA', 'singlet')
                    ham(s1,s2) = ham(s1,s2)&
                      &+ 2.d0*excli(ist1-sta1+1, ist3-sta2+1, ist2-sta1+1, ist4-sta2+1)&
                      &* kdocc(s1) * 0.5d0
                end select

                ! Add correlation term
                ! - W_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
!!! Wrong occupation differences?
                ! Also add occupation difference treatment
                select case(trim(input%xs%bse%bsetype))
                  case('singlet', 'triplet')
                    ham(s1,s2) = ham(s1,s2)&
                      &- sccli(ist1-sta1+1, ist3-sta2+1, ist2-sta1+1, ist4-sta2+1)&
                      &* kdocc(s1) * 0.5d0
                end select

              end do
            end do
          end do
        end do

      ! End loop over(k,kp)-pairs
      end do kp
    end do k

    deallocate(excli, sccli, docc)

    ! Write Info
    write(unitout,*)
    write(unitout, '("Info(bse): invoking lapack routine ZHEEVX")')
    write(unitout, '(" size of BSE-Hamiltonian	   : ", i8)') hamsiz
    write(unitout, '(" number of requested solutions : ", i8)') input%xs%bse%nexcitmax

    ! Allocate eigenvector and eigenvalue arrays
    allocate(beval(hamsiz), bevec(hamsiz, hamsiz))

    ! Set number of excitons (i.e. all of them)
    ne = hamsiz

    ! Diagonalize Hamiltonian
    call timesec(ts0)
    call bsesoldiag(hamsiz, ne, ham, beval, bevec)
    call timesec(ts1)

    ! Deallocate BSE-Hamiltonian
    deallocate(ham)

    write(unitout, '(" timing(in seconds)	   :", f12.3)') ts1 - ts0

    ! Number of excitons to consider
    nexc = hamsiz

    ! Allocate arrays used in spectrum construction
    allocate(oszs(nexc, 3), oszsa(nexc), sor(nexc), pmat(hamsiz, 3))
    allocate(w(input%xs%energywindow%points), spectr(input%xs%energywindow%points))
    allocate(buf(3,3,input%xs%energywindow%points))
    allocate(loss(3, 3, input%xs%energywindow%points),sigma(input%xs%energywindow%points))

    ! Generate an evenly spaced frequency grid 
    call genwgrid(input%xs%energywindow%points, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)

    buf(:,:,:)=zzero

    ! Loop over 3 directions
    optcmp: do oct1 = 1, noptcmp

      ! Allocate momentum matrix
      allocate(pm(3, nstsv, nstsv))

      do iknr = 1, nkptnr

        ! Read momentum matrix for given k-point
        call getpmat(iknr, vkl, 1, nstsv, 1, nstsv, .true., 'PMAT_XS.OUT', pm)

        ! Occupied
        do ist1 = sta1, sto1
          ! Unoccupied
          do ist2 = istl3 , istl3 + nrnst3 -1

            ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
            ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
            !   Where E are quasi-particle energies and e are KS energies.
            if(associated(input%gw)) then
               pm(1:3,ist1,ist2) = pm(1:3,ist1,ist2)&
                 &* (evalsv(ist2,iknr)-evalsv(ist1,iknr))&
                 &/ (eval0(ist2,iknr)- eval0(ist1,iknr))
            end if 

            ! Din
            ! Map the momentum matrix elements to same index convention as
            ! used in the BSE Hamiltonian.
            s1 = hamidx(ist1-sta1+1, ist2-istl3+1,iknr, nrnst1, nrnst3) 
            ! P^i_{o_{s1},u_{s1},k_{s1}}
            pmat(s1, oct1) = pm(oct1, ist1, ist2)

          end do
        end do

      end do

      deallocate(pm)

      ! Calculate oscillators for spectrum
      oszs(:, oct1) = zzero

      ! Exciton index
      do s1 = 1, nexc
        ! Non-reduced k-points
        do iknr = 1, nkptnr
          ! Occupied
          do iv = 1, nrnst1
            ! Unoccupied
            do ic = 1, nrnst3

              s2 = hamidx(iv, ic, iknr, nrnst1, nrnst3)

              ! t^i_s1 = t^i_s1 + A^s1  P^i_{o_{s2},u_{s2},k_{s2}}/(e_{o_{s2},k_{s2}} - e_{u_{s2},k_{s2}})
!!! Wrong occupation differences?
              oszs(s1, oct1) = oszs(s1, oct1) + bevec(s2, s1) * pmat(s2, oct1)&
                &* kdocc(s2)*0.5d0/(evalsv(ic+istl3-1, iknr) - evalsv(iv+sta1-1, iknr))

            end do
          end do
        end do
      end do

      ! Stk: Add case of double grid
      if(dgrid) then 

        write(dotext, '("_SG", i3.3, ".OUT")') iksubpt

        call genfilname(basename='EXCITON', tq0=.true., oc1=oct1, oc2=oct1,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, dotext=dotext, filnam=fnexc)

        call genfilname(basename='EXCITON_SORTED', tq0=.true., oc1=oct1, oc2=oct1,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, dotext=dotext, filnam=fnexcs)

      else

        call genfilname(basename='EXCITON', tq0=.true., oc1=oct1, oc2=oct1,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnexc)

        call genfilname(basename='EXCITON_SORTED', tq0=.true., oc1=oct1, oc2=oct1,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnexcs)

      endif

      ! Write out exciton energies and oscillator strengths
      call getunit(unexc)
      open(unexc, file=fnexc, form='formatted', action='write', status='replace')
      do s2 = 1, hamsiz
        write(unexc, '(i8, 6g18.10)') s2, (beval(s2)+egap-dble(bsed)) * escale,&
          & (beval(s2)+dble(bsed)) * escale, abs(oszs(s2, oct1)), dble(oszs(s2, oct1)),&
          & aimag(oszs(s2, oct1))
      end do
      write(unexc, '("# Nr.  E		      E - E_gap&
        &      |Osc.Str.|      Re      Im")')
      write(unexc, '("# E_gap : ", g18.10)') egap * escale
      if(input%xs%tevout) write(unexc, '("# energies are in electron volts")')
      close(unexc)

      ! Oscillator strengths sorted
      oszsa(:) = abs(oszs(:, oct1))

      call sortidx(hamsiz, oszsa, sor)

      sor = sor(hamsiz:1:-1)

      ! Write sorted oscillator strengths
      open(unexc, file=fnexcs, form='formatted', action='write', status='replace')
      do s1 = 1, hamsiz
        s2 = sor(s1)
        write(unexc, '(i8, 4g18.10)') s1, (beval(s2)+egap-dble(bsed)) * escale,&
          & (beval(s2)+dble(bsed)) * escale, abs(oszs(s2, oct1))
      end do
      write(unexc, '("#	  Nr.	E		  E - E_gap	   |Osc.Str.|")')
      write(unexc, '("# E_gap : ", g18.10)') egap * escale
      if(input%xs%tevout) write(unexc, '("# energies are in electron volts")')
      close(unexc)

    end do optcmp

    ! Stk: If run is within a double grid loop stop here
    notdgrid: if( .not. dgrid) then

      do oct1 = 1, noptcmp

        ! Stk: compute off-diagonal optical components if requested
        if(input%xs%dfoffdiag) then
          octl = 1
          octu = noptcmp
        else
          octl = oct1
          octu = oct1
        end if

        do oct2 = octl, octu

          optcompt = (/ oct1, oct2, 0 /)
          spectr(:) = zzero

          do iw = 1, input%xs%energywindow%points
            do s1 = 1, nexc

              ! Lorentzian lineshape
              spectr(iw) = spectr(iw) + oszs(s1, oct1) * conjg(oszs(s1, oct2))&
                &* (1.d0/(w(iw)-(beval(s1)+egap-bsed)+zi*input%xs%broad))             

              if(input%xs%bse%aresbse) spectr(iw) = spectr(iw)&
                &+ oszs(s1, oct1) * conjg(oszs(s1, oct2))&
                &* (1.d0/(-w(iw)-(beval(s1)+egap-bsed)-zi*input%xs%broad))

            end do
          end do

          spectr(:) = l2int(oct1 .eq. oct2) * 1.d0 - spectr(:) * 8.d0*pi/omega/nkptnr

          buf(oct1,oct2,:)=spectr(:)

        ! End loops over optical components
        end do
      end do

      ! B.A. Why is the loss function constructed with the non symmetrized spectra ?
      call genloss(buf, loss, noptcmp)

      do oct1 = 1, noptcmp

        if(input%xs%dfoffdiag) then
          octl = 1
          octu = noptcmp
        else
          octl = oct1
          octu = oct1
        end if

        do oct2 = octl, octu

          optcompt = (/ oct1, oct2, 0 /)

          ! Generate File names for resulting quantities
          call genfilname(basename='EPSILON', tq0=.true., oc1=oct1, oc2=oct2,&
            & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
            & nar= .not. input%xs%bse%aresbse, filnam=fneps)

          call genfilname(basename='LOSS', tq0=.true., oc1=oct1, oc2=oct2,&
            & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
            & nar= .not. input%xs%bse%aresbse, filnam=fnloss)

          call genfilname(basename='SIGMA', tq0=.true., oc1=oct1, oc2=oct2,&
            & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
            & nar= .not. input%xs%bse%aresbse, filnam=fnsigma)

          call genfilname(basename='SUMRULES', tq0=.true., oc1=oct1, oc2=oct2,&
            & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
            & nar= .not. input%xs%bse%aresbse, filnam=fnsumrules)

          ! Symmetrize the macroscopic dielectric function tensor
          call symt2app(oct1, oct2, input%xs%energywindow%points, symt2, buf, spectr)

          ! Generate optical functions
    ! Genloss generates error in debug mode !!!
    !      call genloss(spectr, loss)
          call gensigma(w, spectr, optcompt, sigma)
          call gensumrls(w, spectr, sumrls)

          ! Write optical functions to file
          call writeeps(iq, oct1, oct2, w, spectr, trim(fneps))
          call writeloss(iq, w, loss(oct1, oct2, :), trim(fnloss))
          call writesigma(iq, w, sigma, trim(fnsigma))
          call writesumrls(iq, sumrls, trim(fnsumrules))

        ! End loop over optical components
        end do
      end do

      ! Store exciton coefficients
      if(associated(input%xs%storeexcitons)) then

        !-----------------------------------------------------
        ! Upon request, store array with exciton coefficients
        !-----------------------------------------------------
        if( (input%xs%storeexcitons%minnumberexcitons .lt. 1) .or. &
          & (input%xs%storeexcitons%minnumberexcitons .gt. hamsiz) .or. &
          & (input%xs%storeexcitons%maxnumberexcitons .lt. 1) .or. &
          & (input%xs%storeexcitons%maxnumberexcitons .gt. hamsiz) .or. &
          & (input%xs%storeexcitons%minnumberexcitons .gt. &
          & input%xs%storeexcitons%maxnumberexcitons) ) then

          write(*,*)
          write(*,'("Error(bse): wrong range of exciton indices: ", 2i5)') &
            & input%xs%storeexcitons%minnumberexcitons,&
            & input%xs%storeexcitons%maxnumberexcitons
          write(*,*)

          stop
        end if  

        ! Write bin
        call getunit(unexc)
        open(unexc, file='EXCCOEFF.bin', action='write',form='unformatted', iostat=stat)

        if((stat/=0) .and. (rank==0)) then
          write(*,*) stat
          write(*,'("Error(bse): error creating EXCCOEFF.bin")')
          write(*,*)
          stop
        end if

        ! Write
        write(unexc) input%xs%storeexcitons%minnumberexcitons,&
          & input%xs%storeexcitons%maxnumberexcitons,&
          & nkptnr, istl3, sta1, sta2, nrnst1, nrnst3, hamsiz

        do ievec = input%xs%storeexcitons%minnumberexcitons,&
          & input%xs%storeexcitons%maxnumberexcitons

          write(unexc) beval(ievec), bevec(1:hamsiz,ievec)
        end do

        close(unexc)

      end if

      deallocate(beval,bevec,oszs,oszsa,sor,pmat,w,spectr,loss,sigma,buf)
      if(associated(input%gw)) deallocate(eval0)

    end if notdgrid
    
    call barrier

  else mpirank

    call barrier

  end if mpirank

contains

  ! hamidx calculates a combined index that labels the states as
  ! hamidx = {o1u1k1, o1u2k1, ..., o1uMk1,
  !      o2uMk1, ..., oMuMk1, o1u1k2, ..., oMuMkN} -> {1,...,M**2N}
  ! if i1=o i2=u
  integer(4) function hamidx(i1, i2, ik, n1, n2)
    implicit none
    integer(4), intent(in) :: i1, i2, ik, n1, n2
    hamidx = i2 + n2 * (i1-1) + n1 * n2 * (ik-1)
  end function hamidx

end subroutine bse
!EOC
