! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bse
! !INTERFACE:
subroutine bse(iqmt)
! !USES:
  ! Basics
  use modinput, only: input
  use mod_constants, only: zone, zi, zzero, pi, h2ev
  use mod_eigenvalue_occupancy, only: evalsv, nstsv
  use mod_kpoint, only: nkptnr
  ! MPI and BLACS/ScaLAPACK
  use modmpi
  use modscl
  ! XS
  use mod_xsgrids
  use modxs, only: vkl0, occsv0, evalsv0, unitout, iqmtgamma
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
  use m_dhesolver
  use m_hesolver
  use m_setup_bse
  use m_genexevec
  use m_putgetexcitons
  use m_makeoscistr
  use m_makespectrum
  use m_invertzmat
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
!   Forked from bse.F90 and adapted for non-TDA and Q-dependent BSE. (Aurich)
!EOP
!BOC

  implicit none

  ! I/O
  integer(4), intent(in) :: iqmt

  ! Local variables
  character(*), parameter :: thisname = "bse"

  integer(4) :: i, j
  integer(4) :: nexc
  real(8) :: ramscale
  real(8) :: ts0, ts1
  logical :: fcoup, fdist, fip

  ! Allocatable arrays
  real(8), allocatable, dimension(:) :: exeval, w
  complex(8), allocatable, dimension(:,:,:) :: symspectr
  complex(8), allocatable, dimension(:,:) :: oscsr

  real(8) :: bsegap
  real(8) :: v1, v2, en1, en2
  integer(4) :: i1, i2, iex1, iex2, nreq
  integer(4) :: nsymcrys_save
  logical :: efind

  ! Distributed arrays
  type(dzmat) :: dham, dexevec, doscsr, dresvec, daresvec
  type(dzmat) :: dcmat, dcpmat

  ! Write out the coupling measures ? 
  logical :: fmeasure
  fmeasure = input%xs%bse%measure

  !! Greeting
  !write(*, '("Info(",a,"): Running at rank", i3)')&
  !  & trim(thisname), mpiglobal%rank

  ! Independent particle formalism is used
  if(input%xs%bse%bsetype == "IP") then
    fip = .true.
  else
    fip = .false.
  end if

  ! Include coupling terms
  fcoup = input%xs%bse%coupling

  ! Distribute Hamilton matrix if ScaLapack is used
  ! and distribution is desired
  fdist = .false.
#ifdef SCAL
  fdist = input%xs%bse%distribute
#endif

  if(associated(input%gw)) then 
    ! If scissor correction is presented, one should nullify it
    input%xs%scissor=0.0d0
  end if

  !---------------------------------------------------------------------------!
  ! Find occupation limits using k, k-qmt/2 and k+qmt/2 grids
  !---------------------------------------------------------------------------!
  ! Read Fermi energy from file
  ! (only needed in setranges_modxs::findocclims to check whether system has a gap)
  ! Use EFERMI_QMT001.OUT
  call genfilname(iqmt=iqmtgamma, setfilext=.true.)
  call readfermi
  ! Set ist* variables and ksgap in modxs using findocclims
  ! Needs: init2 (is called by bselauncher)
  ! On exit: - k-grid quantities are stored in modxs: vkl0, evalsv0 etc
  !            k-qmt/2 grid quantities are stored in default locations: vkl, evalsv etc
  !            This is altered by select_transitions below
  !          - Saves mappings k -> k+qmt/2, k -> k-qmt/2 and k-qmt/2 -> k+qmt/2 
  !            to modbse
  call printline(unitout, '-')
  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Inspecting occupations...'
  call flushifc(unitout)
  
  call setranges_modxs(iqmt)
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  ! On top of GW
  !---------------------------------------------------------------------------!
  !
  ! WARNING: This still is inconsistent, since occupation search is
  !          not done on EVALQP.OUT
  ! WARNING: Only works for Qmt=0, in this case vkl0 = vkl
  !
  if(associated(input%gw) .and. iqmt==1) then

    ! Save KS eigenvalues of the k-grid to use them later for renormalizing PMAT
    if(allocated(eval0)) deallocate(eval0)
    allocate(eval0(nstsv, nkptnr))
    eval0=evalsv0

    ! Read QP Fermi energies and eigenvalues from file
    ! NOTE: QP evals are shifted by -efermi-eferqp with respect to KS evals
    ! NOTE: getevalqp sets mod_symmetry::nsymcrys to 1
    ! NOTE: getevalqp needs the KS eigenvalues as input
    nsymcrys_save = nsymcrys
    call getevalqp(nkptnr, vkl0, evalsv)
    nsymcrys = nsymcrys_save

    ! Set k and k'=k grid eigenvalues to QP energies
    evalsv0=evalsv

    write(unitout,'("Info(",a,"):&
      & Quasi particle energies are read from EVALQP.OUT")') trim(thisname)

  else if(associated(input%gw) .and. iqmt /= 1) then 

    if(bicurrent%isroot) then 
      write(*,'("Error(",a,"):&
       & BSE+GW only supported for 0 momentum transfer.")') trim(thisname)
    end if
    call terminate

  end if
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  ! Setup participating transitions and index maps
  !---------------------------------------------------------------------------!
  ! Select relevant transitions for the construction
  ! of the BSE Hamiltonian.
  ! Sets up transition energies in array "de", combined index maps and 
  ! other helpers in modbse
  !
  ! Note: Operates on mpiglobal
  ! Note: On exit k-grid/2 quantities are stored in vkl0, evalsv0 etc and
  !       k+qmt/2 grid quantities are stored in vkl, evalsv etc
  call printline(unitout, '-')
  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Selecting transitions...'
  call flushifc(unitout)

  if(fdist) then 
    ! Use all ranks on mpiglobal
    call select_transitions(iqmt, serial=.false.)
  else
    ! Use all current rank only
    call select_transitions(iqmt, serial=.true.)
  end if

  ! Determine "BSE"-gap, i.e. the lowest energy 
  ! KS transition with momentum transfer qmt which is considered.
  ! (includes scissor).
  bsegap = minval(de)
  ! Write some Info
  call printline(unitout, '-')
  write(unitout,*)
  write(unitout, '("Info(",a,"):&
    & bsegap-scissor, bsegap (eV):", E23.16,1x,E23.16)')&
    & trim(thisname), (bsegap-sci)*h2ev, bsegap*h2ev
  !---------------------------------------------------------------------------!

  ! Independent particle approximation needs no solution of any EVP
  if(fip) then 

    write(unitout, '("Info(",a,"):&
      & IP requested, nothing much to do")') trim(thisname)

    nexc = hamsize

    ! "Exciton" energies are just the KS transition energies
    allocate(exeval(nexc))
    exeval(1:hamsize) = de(ensortidx) 

  ! Do the actual BSE 
  else

    ! Process is on the process grid
    if(bicurrent%isactive) then 

      ! Define global distributed Hamiltonian matrix.
      ! Note: If the context of bicurrent is -1, i.e. the rank is not on the blacs
      !       process grid 
      ! Note: If not compiled with -DSCAL this just allocates a local array.
      call new_dzmat(dham, hamsize, hamsize, bicurrent)

      !------------------------------------------------------------------------!
      ! Write info
      !------------------------------------------------------------------------!
      write(unitout, '("Info(",a,"): Assembling BSE matrix")')&
        & trim(thisname)
      write(unitout, '("  RR/RA blocks of global BSE-Hamiltonian:")')
      write(unitout, '("  Size of the Hamiltonian:",i8)') hamsize
      write(unitout, '("  Number of k-points:", i8)') nk_bse
      ramscale=1.0d0

      if(fcoup) then
        write(unitout, '(" Including coupling terms ")')
        write(unitout, '(" Using squared EVP")')
        ramscale=3.0d0
      end if

      write(unitout, '("  Distributing matrix to ",i3," processes")') bicurrent%nprocs
      write(unitout, '("  Local matrix shape ",i6," x",i6)')&
        & dham%nrows_loc, dham%ncols_loc

      write(unitout, '("Info(",a,"):&
        & max local RAM needed for BSE matrices ~ ", f12.6, " GB" )')&
        & trim(thisname),&
        & dble(ramscale)*dham%nrows_loc*dham%ncols_loc*16.0d0/1024.0d0**3
      write(unitout,*)
      !------------------------------------------------------------------------!
      call printline(unitout, '-')
      !------------------------------------------------------------------------!
      ! Assemble Hamiltonian matrix
      !------------------------------------------------------------------------!
      call timesec(ts0)

      if(fcoup) then 
        ! Get aux. EVP matrices S=C(RR+RA)C and C=(RR-RA)^1/2 
        call setup_bse_tr_dist(iqmt, bicurrent, smat=dham, cmat=dcmat)
      else
        ! Get TDA Hamiltonian H=H^RR
        call setup_bse_block_dist(dham, iqmt, .false., bicurrent)
      end if

      call timesec(ts1)
      write(unitout, '("All processes build their local matrix")')
      if (input%xs%BSE%outputlevelnumber == 1) &
        & write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0
      !------------------------------------------------------------------------!

      !------------------------------------------------------------------------!
      ! Solve EVP or squared EVP
      !------------------------------------------------------------------------!

      ! Print some info
      if(fcoup) then 
        write(unitout, '("Info(",a,"):&
          & Solving hermitian squared EVP")') trim(thisname)
      else 
        write(unitout, '("Info(",a,"):&
          & Diagonalizing RR Hamiltonian (TDA)")') trim(thisname)
      end if
#ifdef SCAL
      if(fdist) then
        write(unitout, '("Info(",a,"):&
          & Invoking scalapack routine PZHEEVX")') trim(thisname)
      else
        write(unitout, '("Info(",a,"):&
          & Invoking Lapack routine ZHEEVR")') trim(thisname)
      end if
#else
      write(unitout, '("Info(",a,"):&
        & Invoking Lapack routine ZHEEVR")') trim(thisname)
#endif

      ! Eigenvectors are distributed
      ! must be NxN because ScaLapack solver expects it so
      call new_dzmat(dexevec, hamsize, hamsize, bicurrent) 

      ! Eigenvalues are global
      allocate(exeval(hamsize))
      exeval = 0.0d0

      ! Find only eigenvalues relevant for requested spectrum?
      efind = input%xs%bse%efind

      ! Diagonalize Hamiltonian (destroys the content of ham)
      call timesec(ts0)

      ! Use Lapack instead of ScaLapck if only current 
      ! rank diagonalizes, i.e. bi0d is used. (bit of a workaround)
      if(fdist .eqv. .false.) then 
        dham%isdistributed = .false.
        dexevec%isdistributed = .false.
      end if

      ! Find solutions in energy window only
      if(efind) then

        if(fcoup) then 
          v1=(max(wl, 0.0d0))**2
          v2=(wu)**2
        else
          v1=max(wl, 0.0d0)
          v2=wu
        end if

        call dhesolver(dham, exeval, bicurrent, dexevec,&
         & v1=v1, v2=v2, found=nexc,&
         & eecs=input%xs%bse%eecs)

      ! Find the nexc lowest solutions
      else

        if(input%xs%bse%nexc == -1) then 
          i2 = hamsize
        else
          i2 = input%xs%bse%nexc
        end if
        i1 = 1

        call dhesolver(dham, exeval, bicurrent, dexevec,&
         & i1=i1, i2=i2, found=nexc,&
         & eecs=input%xs%bse%eecs)

      end if

      ! Set it back (see above)
      if(fdist .eqv. .false.) then 
        dham%isdistributed = .true.
        dexevec%isdistributed = .true.
      end if

      ! Square root of EVs
      if(fcoup) then 
        ! Take square root of auxiliary eigenvalues
        ! to retrieve actual eigenvalues
        if(any(exeval < 0.0d0)) then 
          write(*, '("Error(",a,"): Negative squared EVP evals occured")')&
            & trim(thisname)
          write(*,'(E23.16)') exeval
          call terminate
        end if
        exeval = sqrt(exeval)
      end if

      ! Deallocate BSE-Hamiltonian
      call del_dzmat(dham)

      ! Stop timer for diagonalization 
      call timesec(ts1)

      call blacsbarrier(bicurrent)

      ! Print some info
      if(efind) then
        if(fcoup) then 
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
        write(unitout, '("         i1=",i8," i2=",i8)') i1, i2
      end if
      if (input%xs%BSE%outputlevelnumber == 1) &
        & write(unitout, '("  Timing (in seconds)	   :", f12.3)') ts1 - ts0
      write(unitout,*)

      if(nexc < 1) then 
        write(unitout, '("No eigenvalues found the energy window, check input")')
        call terminate
      end if  
      call printline(unitout, '-')
      !------------------------------------------------------------------------!

      ! Testing
      if(fmeasure) then 
        if(bicurrent%isroot) then 
          write(unitout, '("Info(",a,"): Writing coupling measures to file")')&
            & trim(thisname)
          call writemeasures(iqmt, nexc, exeval, fcoup)
        end if
      end if

      ! =================================================================!
      ! Store excitonic energies and wave functions to file if requested !
      ! =================================================================!
      store: if(associated(input%xs%storeexcitons)) then

        ! Select which solutions to store
        !   Select by energy
        if(input%xs%storeexcitons%selectenergy) then 
          en1=input%xs%storeexcitons%minenergyexcitons
          en2=input%xs%storeexcitons%maxenergyexcitons
          if(input%xs%storeexcitons%useev) then 
            en1=en1/h2ev
            en2=en2/h2ev
          end if
          call energy2index(hamsize, nexc, exeval, en1, en2, iex1, iex2)
        !   Select by number
        else
          iex1=input%xs%storeexcitons%minnumberexcitons
          iex2=input%xs%storeexcitons%maxnumberexcitons
        end if
        nreq=iex2-iex1+1

        ! Check requested range
        if(nreq < 1 .or. nreq > nexc .or. iex1<1 .or. iex2<1 ) then
          if(bicurrent%isroot) then 
            write(*,'("Error(",a,"):&
              & storeexcitons index mismatch.")') trim(thisname)
            write(*,*) "iex1, iex2, nreq, nex", iex1, iex2, nreq, nexc
          end if
          call terminate
        end if

        write(unitout, '("Info(",a,"):&
          & Writing excition eigenvectors for index range=",2i8)')&
          & trim(thisname), iex1, iex2

        ! When full BSE with time reversal symmetry is used 
        ! the original eigenvectors need to be reconstructed from the 
        ! auxiliary ones.
        if(fcoup) then 

          write(unitout, '("Info(",a,"): Generating resonant&
            & and anti-resonant exciton coefficients from&
            & auxilliary squared EVP eigenvectors (pos. E).")') trim(thisname)

          !===========================================================!
          ! Make C' Matrix                                            !
          ! Cpmat = (A-B)^{-1/2} = C^-1                               !
          ! Cpmat is needed additionally to build eigenvectors        !
          !===========================================================!
          write(unitout, '("Info(",a,"):&
            & Inverting C=(RR-RA)^1/2 matrix")') trim(thisname)
          call timesec(ts0)
          call new_dzmat(dcpmat, hamsize, hamsize, bicurrent)
          do j = 1, dcmat%ncols_loc
            do i = 1, dcmat%nrows_loc
              dcpmat%za(i,j) = dcmat%za(i,j)
            end do
          end do
          call dzinvert(dcpmat)
          call timesec(ts1)
          if (input%xs%BSE%outputlevelnumber == 1) &
            & write(unitout, '("  Time needed",f12.3,"s")') ts1-ts0
          !===========================================================!

          ! Make selected eigenvectors
          ! Note: allocates dresvec, daresvec and deallocates dcpmat
          call gendexevec(iex1, iex2, nexc, dcmat, dcpmat, dexevec, exeval,&
            & dresvec, daresvec)

          ! Write to binary file
          write(unitout, '("Info(",a,"):&
            & Writing exciton eigenvectors to file.")') trim(thisname)
          call putd_excitons(exeval(iex1:iex2), drvec=dresvec, davec=daresvec,&
            & iqmt=iqmt, a1=iex1, a2=iex2)

          ! Explicit resonant/anti-resonant eigenvectors no longer needed
          call del_dzmat(dresvec)
          call del_dzmat(daresvec)

        ! For TDA we can directly write out the desired eigenvectors
        else if(.not. fcoup) then  

          write(unitout, '("Info(",a,"):&
            & Writing exciton eigenvectors to file.")') trim(thisname)
          call setview_dzmat(dexevec, hamsize, nreq, 1, iex1)
          call putd_excitons(exeval(iex1:iex2), dexevec,&
            & iqmt=iqmt, a1=iex1, a2=iex2)
          call setview_dzmat(dexevec, hamsize, nexc, 1, 1)

        end if

      end if store
      ! =================================================================!

    ! Process is not on 2d process grid
    else
      write(*,*) "Info(",trim(thisname),"): rank,", rank," is not on the grid"
    end if

  ! IP or not if
  end if

  ! Note: The following allocation is only needed for processes that are not
  !       on the process grid and serves the only purpose that the debugger
  !       is not complaining here.
  if(.not. allocated(exeval)) allocate(exeval(hamsize))

  ! Calculate oscillator strengths.
  ! Note: Deallocates dexevec and dcmat
  if(fcoup) then 
    call makeoscistr_dist(iqmt, nexc, dexevec, doscsr, bicurrent, exeval, dcmat)
  else
    call makeoscistr_dist(iqmt, nexc, dexevec, doscsr, bicurrent)
  end if

  ! Back to the process grid.
  if(bicurrent%isactive) then 

    ! Every process gets a copy of the oscillator strength
    ! (actually only rank 0 writes them to file, but is is not much 
    !  memory and it make the setup for the spectrum calculation easier) 
    call dzmat_send2global_all(oscsr, doscsr, bicurrent)
    if(fip) then
      oscsr = oscsr(ensortidx,:)
    end if

    if(bicurrent%isroot) then
      ! Write excition energies and oscillator strengths to 
      ! text file. 
      write(unitout,*)
      call printline(unitout, '-')
      write(unitout, '("Info(",a,"):&
        & Writing excition energies and oscillator strengths to text file.")')&
        & trim(thisname)
      call writeoscillator(hamsize, nexc, nk_bse, -bsegap, exeval, oscsr, iqmt=iqmt)
    end if

    ! Allocate arrays used in spectrum construction
    if(bicurrent%isroot) then
      allocate(symspectr(3,3,nw))
    end if

    ! Only process root gets an acctual output for symspectr
    call makespectrum_dist(iqmt, nexc, nk_bse, exeval, oscsr, symspectr, bicurrent)

    ! Only root writes derived outputs
    if(bicurrent%isroot) then

      ! Allocate arrays used in spectrum construction
      allocate(w(nw))
      ! Generate an evenly spaced frequency grid 
      call genwgrid(nw, input%xs%energywindow%intv,&
        & input%xs%tddft%acont, 0.d0, w_real=w)

      ! Generate and write derived optical quantities
      call writederived(iqmt, symspectr, nw, w)

    end if

    ! Clean up
    deallocate(exeval, oscsr)
    if(associated(input%gw)) deallocate(eval0)

    if(bicurrent%isroot) then 
      deallocate(w)
      deallocate(symspectr)
    end if

    ! MPI
    ! Ranks that are on the BLACS grid signal that they are done
    if(fdist) then 
      call barrier(callername=thisname)
    end if

  ! Not on grid
  else

    write(*, '("Warning(",a,"): Rank", i4, " is idle.")')&
      & trim(thisname), mpiglobal%rank

    ! Ranks that are not on the BLACS grid wait 
    if(fdist) then 
      call barrier(callername=thisname)
    end if

  end if

end subroutine bse
!EOC
