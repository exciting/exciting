!BOP
! !MODULE: modbse
! !DESCRIPTION:
!   Supporting global variables and routines for the BSE scope.
!
! !REVISION HISTORY:
!   Created 2016 (Aurich)
!EOP   
!BOC
module modbse
  use modmpi
#ifdef USEOMP
  use omp_lib
#endif
  use modinput, only: input
  use mod_constants, only: h2ev
  use modxs, only: unitout, totalqlmt
  use modxs, only: evalsv0, occsv0
  use modxs, only: istocc0, istocc, istunocc0, istunocc,&
                 & isto0, isto, istu0, istu, ksgapval, qgap, iqmtgamma
  use mod_eigenvalue_occupancy, only: evalsv, occsv, nstsv
  use mod_kpoint, only: nkptnr, nkpt
  use m_getunit

  implicit none

  ! Eigenvaulues and eigenvectors for
  ! k-qmt, k, k+qmt
  real(8), pointer :: p_eval_plus(:), p_eval(:), p_eval_minus(:)
  complex(8), pointer :: p_evec_plus(:,:), p_evec(:,:), p_evec_minus(:,:)
  ! Reference states for occupied and unoccupied state indices
  integer(4) :: ioref, iuref 
  ! Size of the resonant-resonant block of the Hamiltonian
  integer(4) :: hamsize
  ! Number of all k-points and k-points/k-k'-combinations
  ! to be considered in BSE
  integer(4) :: nk_max, nk_bse, nkkp_bse
  ! Number of occupied and unoccupied bands present
  integer(4) :: no_max, nu_max, nou_max
  ! Maximum of needed number of o/u bands over all needed k
  integer(4) :: no_bse_max, nu_bse_max, nou_bse_max
  ! Energy range of spectrum
  real(8) :: wl, wu
  integer(4) :: nw
  ! Cutoff for occupation
  real(8) :: cutoffocc
  ! Auto-select KS transitions
  logical :: fensel
  ! Convergence energy
  real(8) :: econv(2)
  ! Lorentian width cutoff
  ! BSE gap
  real(8) :: sci
  ! Occupation factors for Hamiltonian construction
  real(8), allocatable :: ofac(:)
  ! KS/QP energy differences
  real(8), allocatable :: de(:)
  ! KS eigenvalue in case of GW QP energie usage
  real(8), allocatable, dimension(:,:) :: eval0
  ! Combined BSE index map
  integer(4), allocatable :: smap(:,:)
  integer(4), allocatable :: smap_rel(:,:)
  ! Energy sorting of smap
  integer(4), allocatable :: ensortidx(:)
  ! Contributing k points
  ! relative index ik -> global index iknr
  integer(4), allocatable :: kmap_bse_rg(:)
  ! global index iknr -> relative index ik
  integer(4), allocatable :: kmap_bse_gr(:)
  ! Number of allowed o-u transitions at each k-point
  integer(4), allocatable :: kousize(:)
  ! Lower and upper bound for occupied and unoccupied band indices for each k
  integer(4), allocatable :: koulims(:,:)

  ! Mapping between k and k' grids
  integer(4), dimension(:), allocatable :: ik2ikqmtp, ik2ikqmtm, ikqmtm2ikqmtp
  ! Offsets of k grids
  real(8), dimension(3) :: vkloff, vkqmtploff, vkqmtmloff

  ! Testing array for coupling measures
  real(8), allocatable :: vwdiffrr(:), vwdiffar(:)

  ! Filebasenames
  character(256) :: infofbasename = "BSEINFO"
  character(256) :: scclifbasename = "SCCLI"
  character(256) :: scclicfbasename = "SCCLIC"
  character(256) :: exclifbasename = "EXCLI"
  character(256) :: infofname
  character(256) :: scclifname
  character(256) :: exclifname


  contains

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Routines to setup the combinded index of  ! 
    ! the BSE hamiltonian.                      !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: setranges_modxs
    ! !INTERFACE:
    subroutine setranges_modxs(iqmt)
      use modinput
      use mod_xsgrids
      use mod_Gkvector, only: gkmax
      use modxs, only: totalqlmt, evalsv0, usefilext0, filext0,&
                     & ksgap, ksgapval, qmtpgap, qmtmgap, unitout
      use m_genfilname
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iqmt ! Considered q-point index (must be on the not shifted k-grid)
    !
    ! !DESCRIPTION:
    !   A small wrapper for the routine {\tt findocclims}. Used to initialize
    !   modxs module variables for occupation limits for $\vec{k}$ and $\vec{k}+\vec{q}$ 
    !   in some BSE related routines. 
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none
      integer(4), intent(in) :: iqmt

      integer(4) :: iomax, iumin

      integer(4) :: iomax_kkqmtp, iumin_kkqmtp
      integer(4) :: iomax_kkqmtm, iumin_kkqmtm

      integer(4), dimension(:), allocatable :: io_k, iu_k
      integer(4), dimension(:), allocatable :: io_kqmtp, iu_kqmtp
      integer(4), dimension(:), allocatable :: io_kqmtm, iu_kqmtm 

      real(8), parameter :: epslat = 1.0d-8

      logical :: fgap, fsamek
      real(8) :: gap
      real(8) :: t0,t1

      call timesec(t0)
      !---------------------------------------------------!
      ! Get offsets of and mapping between k and k' grids !
      !---------------------------------------------------!
      ! Note: requires init2 to set up totalqlmt
      call xsgrids_init(totalqlmt(1:3,iqmt), gkmax)

      !! Offsets for k and k' grids
      ! k
      vkloff = k_kqmtp%kset%vkloff
      ! k+qmt/2
      vkqmtploff = k_kqmtp%kqmtset%vkloff
      ! k-qmt/2
      vkqmtmloff = k_kqmtm%kqmtset%vkloff

      !! Mappings between k and k' grids
      ! k --> k+qmt/2
      if(allocated(ik2ikqmtp)) deallocate(ik2ikqmtp)
      allocate(ik2ikqmtp(nkpt))
      ik2ikqmtp(:) = k_kqmtp%ik2ikqmt(:)
      ! k --> k-qmt/2
      if(allocated(ik2ikqmtm)) deallocate(ik2ikqmtm)
      allocate(ik2ikqmtm(nkpt))
      ik2ikqmtm(:) = k_kqmtm%ik2ikqmt(:)
      ! k-qmt/2 --> k+qmt/2
      if(allocated(ikqmtm2ikqmtp)) deallocate(ikqmtm2ikqmtp)
      allocate(ikqmtm2ikqmtp(nkpt))
      ikqmtm2ikqmtp(:) = ikm2ikp(:)

      call xsgrids_finalize()
      !---------------------------------------------------!

      !---------------------------------------------------!
      ! Inspect occupancies for k, k'=k+qmt/2             !
      !---------------------------------------------------!
      ! Reset mod_kpoint / mod_Gkvector variables to the unshifted k-grid
      ! (apart from xs%vkloff)
      call init1offs(vkloff)

      ! Save the k-grid to modxs vkl0 etc 
      call xssave0

      ! Set k and G+k variables in the standard locations mod_kpoint and mod_Gkvector 
      ! to those of the k' grid.
      call init1offs(vkqmtploff)
      
      write(unitout, '("Info(setranges_modxs): Determining k+qmt/2 grid")') 
      if(all(abs(vkloff-vkqmtploff) < epslat)) then 
        if (iqmt .ne. 1) then
          write(unitout, '("Info(setranges_modxs): (+) Same k-grids for iqmt=1 and iqmt=",i3)') iqmt
        end if
        fsamek=.true.
      else
        write(unitout, '("Info(setranges_modxs): (+) Different k-grids for iqmt=1 and iqmt=",i3)') iqmt
        fsamek=.false.
      end if

      ! Allocate k eigenvalues
      ! Note: If evalsv0 is allocated before findocclims is called it gets read in 
      ! and stays allocated.
      if(.not. allocated(evalsv0)) allocate(evalsv0(nstsv, nkpt))

      ! Explicitly specify k file extension
      usefilext0 = .true.
      ! Set EVALSV_QMT001.OUT as first reference for the occupation search
      call genfilname(iqmt=iqmtgamma, fileext=filext0)

      ! Set k' file extension
      if(fsamek) then 
        ! Set EVALSV_QMT001.OUT as second reference for the occupation search
        call genfilname(iqmt=iqmtgamma, setfilext=.true.)
      else
        ! Set EVALSV_QMTXYZ.OUT as second reference for the occupation search
        call genfilname(iqmt=iqmt, setfilext=.true.)
      end if

      allocate(io_k(nkpt), iu_k(nkpt))
      allocate(io_kqmtp(nkpt), iu_kqmtp(nkpt))

      ! Inspect occupation limits for k and k+qmt/2 grids
      call findocclims(iqmt, ik2ikqmtp, iomax_kkqmtp, iumin_kkqmtp,&
        & io_k, io_kqmtp, iu_k, iu_kqmtp)

      ! Save gap status
      !   Was an (indirect) gap found?
      fgap = ksgap
      !   Size of found gap
      gap = ksgapval
      !   Size of qmt dependent gap (with qmt=0 this is the direct gap)
      qmtpgap = qgap

      deallocate(io_k, iu_k)
      deallocate(io_kqmtp, iu_kqmtp)
      !---------------------------------------------------!

      !---------------------------------------------------!
      ! Inspect occupancies for k, k'=k-qmt/2             !
      !---------------------------------------------------!
      ! Reset mod_kpoint / mod_Gkvector variables to the unshifted k-grid
      ! (apart from xs%vkloff)
      call init1offs(vkloff)

      ! Save the k-grid to modxs vkl0 etc 
      call xssave0

      ! Set k and G+k variables in the standard locations mod_kpoint and mod_Gkvector 
      ! to those of the k' grid.
      call init1offs(vkqmtmloff)

      write(unitout, '("Info(setranges_modxs): Determining k-qmt/2 grid")') 
      
      if(all(abs(vkloff-vkqmtmloff) < epslat)) then 
        if (iqmt .ne. 1) then
          write(unitout, '("Info(setranges_modxs): (-) Same k-grids for iqmt=1 and iqmt=",i3)') iqmt
        end if
        fsamek=.true.
      else
        write(unitout, '("Info(setranges_modxs): (-) Different k-grids for iqmt=1 and iqmt=",i3)') iqmt
        fsamek=.false.
      end if

      ! Explicitly specify k file extension
      usefilext0 = .true.
      ! Set EVALSV_QMT001.OUT as first reference for the occupation search
      call genfilname(iqmt=iqmtgamma, fileext=filext0)

      ! Set k' file extension
      if(fsamek) then 
        ! Set EVALSV_QMT001.OUT as second reference for the occupation limits search
        call genfilname(iqmt=iqmtgamma, setfilext=.true.)
      else
        ! Set EVALSV_QMTXYZ_m.OUT as second reference for the occupation limits search
        call genfilname(iqmt=iqmt, auxtype="m", setfilext=.true.)
      end if

      allocate(io_k(nkpt), iu_k(nkpt))
      allocate(io_kqmtm(nkpt), iu_kqmtm(nkpt))

      call findocclims(iqmt, ik2ikqmtm, iomax_kkqmtm, iumin_kkqmtm,&
        & io_k, io_kqmtm, iu_k, iu_kqmtm)
      !---------------------------------------------------!

      ! Refine gap status
      !   Does the system still have a gap with this different choice of the
      !   k mesh?
      fgap = ksgap .and. fgap
      !   Was a smaller indirect gap found?
      gap = min(ksgapval, gap)
      !   Gap with -qmt direction
      qmtmgap = qgap

      ! Use the highest partially occupied band over all k, k+qmt/2, k-qmt/2
      iomax = max(iomax_kkqmtp, iomax_kkqmtm)
      ! Use the lowest partially unoccupied band over all k, k+qmt/2, k-qmt/2
      iumin = min(iumin_kkqmtp, iumin_kkqmtm)

      deallocate(io_k, iu_k)
      deallocate(io_kqmtm, iu_kqmtm)

      ! Set gap status in modxs
      ksgap = fgap
      ksgapval = gap

      ! Setting modxs variables for occupation limits
      istocc0 = iomax
      istocc = iomax
      istunocc0 = iumin
      istunocc = iumin

      ! Set additional modbse variables
      iuref = iumin 
      ioref = 1
      nk_max = nkptnr
      no_max = iomax
      nu_max = nstsv-iumin+1
      nou_max = no_max*nu_max

      ! Scissor
      sci = input%xs%scissor

      write(unitout, '("Info(setranges_modxs):&
        & Number of non-reduced k-points:",i9)') nk_max
      write(unitout, '("Info(setranges_modxs):&
        & Number of states considered:",i9)') nstsv
      write(unitout, '("Info(setranges_modxs):&
        & Number of (partially) occupied state:", i9)') no_max
      write(unitout, '("Info(setranges_modxs):&
        & Highest (partially) occupied state:", i9)') iomax
      write(unitout, '("Info(setranges_modxs):&
        & Number of (partially) unoccupied state:", i9)') nu_max
      write(unitout, '("Info(setranges_modxs):&
        & Lowest (partially) unoccupied state:", i9)') iumin

      if(ksgapval == 0.0d0) then
        write(unitout, '("Warning(setranges_modxs): The system has no gap")')
        if(sci /= 0.0d0) then 
          write(unitout, '("Warning(setranges_modxs):&
            &   Scissor > 0 but no gap. Setting scissor to 0.")')
          sci = 0.0d0
        end if
      end if  

      call timesec(t1)
      if (input%xs%BSE%outputlevelnumber == 1) then
        write(unitout, '("Info(setranges_modxs):&
          & Time needed/s = ", f12.7)') t1-t0
      end if
    end subroutine setranges_modxs
    !EOC

    !BOP
    ! !ROUTINE: select_transitions
    ! !INTERFACE:
    subroutine select_transitions(iqmt, serial, dirname)
      use mod_kpoint, only: vkl
      use modxs, only: usefilext0, filext0, vkl0
      use modxas, only: xasstart, xasstop, ecore
      use m_genfilname
      use mod_symmetry, only: nsymcrys
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iqmt  ! q-point index (on unshifted k-mesh)
    ! Module out:
    ! integer(4) :: hamsize            ! Dimension of the BSE Hamiltonian matrix
    ! real(8)    :: ofac(hamsize)      ! Occupation factors need for 
    !                                  ! the construction of the BSE matrix
    ! integer(4) :: smap(hamsize, 3)   ! Map between BSE matrix index and u,o,k indices
    ! integer(4) :: kousize(nk)        ! How many u o combinations allowed for each k
    !
    ! !DESCRIPTION:
    !   Given a selected energy range for the spectrum, this routine will 
    !   select relevant transitions for each k point. Appart from the KS transition
    !   energies the routine checks whether the transition contain problematic 
    !   occupancy differences and sorts them out if need be.
    !   The simple treatment of fractional occupancy does not allow for transitions
    !   between states of the same partial occupancy. Also cases of occupancy inversion
    !   where the occupancy difference is negative are filtered out, since those break
    !   any kind of hermiticity of the BSE Hamiltonian.\\
    !   The routine crates the compined index map
    !   $\alpha \leftrightarrow \{u_\alpha, o_\alpha, \vec{k}_\alpha\}$, auxilliary maps
    !   and determins the size of the resulting hamiltonian.
    !   In all cases is u the fastest index followed by o and k.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none 

      integer(4), intent(in) :: iqmt
      logical, intent(in), optional :: serial
      character(*), intent(in), optional :: dirname

      logical :: fserial
      integer(4) :: ik, ikqp, ikqm, s, iknr
      integer(4) :: io, iu, kous
      integer(4) :: gwiomin, gwiumax
      integer(4) :: io1, io2, iu1, iu2
      integer(4) :: iomax, iumax
      integer(4) :: iomin, iumin
      real(8) :: maxocc
      real(8), parameter :: epslat = 1.0d-8

      integer(4) :: nk_loc, hamsize_loc
      integer(4) :: nsymcrys_save
      integer(4), allocatable :: smap_loc(:,:)
      real(8) :: detmp, doctmp
      real(8), allocatable :: ofac_loc(:)
      real(8), allocatable :: de_loc(:)
      real(8) :: t0, t1
      logical :: fsamek0, fsamek1, fsamek
      logical, allocatable :: sflag(:)
      integer(4) :: k1, k2
      integer(4) :: i1, i2
      logical :: fxas, posdiff
      character(*), parameter :: thisname = "select_transitions"

      call timesec(t0)

      ! Check whether mpi is used for the selection
      if(present(serial)) then 
        fserial = serial
      else
        fserial = .false.
      end if

      ! Check is XAS is used
      if(input%xs%bse%xas) then
        fxas = .true.
      else
        fxas = .false.
      end if

      ! Search for needed IP/QP transitions automatically
      ! depending on the chosen energy window.
      if((any(input%xs%bse%nstlbse == 0) .and. .not. fxas) &
       & .or. (any(input%xs%bse%nstlxas ==0) .and. fxas)) then
        fensel = .true.
      else
        fensel = .false.
      end if
      ! Set maxocc factor: 2.0d0 for spin-unpolarized, 1.0d0 for spin-polarized
      if (input%groundstate%tevecsv) then
        maxocc=1.0d0
      else
        maxocc=2.0d0
      end if

      ! What energy range is of interest?
      nw = input%xs%energywindow%points
      wl = input%xs%energywindow%intv(1)
      wu = input%xs%energywindow%intv(2)
      if(wl < 0.0d0 .or. wu < 0.0d0 .or. wu-wl < 0.0d0) then
        write(*,*) "Error(select_transitions): Inproper energy interval", wl, wu
        call terminate
      end if

      ! Excitons are build from KS states and the KS transition energies 
      ! dominatly determine exciton energies. 
      ! Select KS transitions withing the energy energy window for the 
      ! spectrum plus extra convergence energy. (referenced only if fensel)
      econv = input%xs%bse%econv

      if((wu+econv(2)-max(wl+econv(1),0.0d0)) < 0.0d0) then 
        write(*,*) "Error(select_transitions): Conflicting econv", econv(1), econv(2)
        call terminate
      end if

      ! What is considered to be occupied
      cutoffocc = input%groundstate%epsocc
      
      ! Bands to inspect
      if(fensel) then

        ! By default use all available XS GS states for the search
        if(fxas) then
          io1=xasstart
          io2=xasstop
          iu1=istunocc0
          iu2=nstsv
        else 
          io1 = 1
          io2 = istocc0 ! highest partially occupied "band" over k_+, k_- and k
          iu1 = istunocc0 ! lowest partially unoccupied band over k_+, k_- and k
          iu2 = nstsv
        end if

        ! If GW QP energies were computed, restrict the search
        ! to the computed GW range, also finite momentum transfer
        ! is not yet included.
        if(associated(input%gw) .and. iqmt==1) then

          gwiomin = input%gw%ibgw
          gwiumax = input%gw%nbgw
          io1 = gwiomin
          iu2 = gwiumax

          write(unitout,'("Info(",a,"):&
            & Energy selection ontop of GW.")') trim(thisname)

        end if

      ! Use the specified bands, and only inspect occupations.
      else

        if(fxas) then
          io1 = xasstart
          io2 = xasstop
          iu1 = input%xs%bse%nstlxas(1)+istunocc0-1
          iu2 = input%xs%bse%nstlxas(2)+istunocc0-1
        else
          io1= input%xs%bse%nstlbse(1)
          io2= input%xs%bse%nstlbse(2)
          iu1= input%xs%bse%nstlbse(3)+istunocc0-1
          iu2= input%xs%bse%nstlbse(4)+istunocc0-1
        end if

      end if

      if(istunocc0 < io1) then
        write(*, '("Waring(select_transitions):", a, 2i4)') &
          & "Lowest (partially) unoccupied band is below the first &
           considered (partially) occupied band.", istunocc0, io1
      end if

      write(unitout, '("Info(select_transitions): Inspecting transitions in the &
        & band interval:")')
      !write(unitout, '("  io1:", i4, " io2:", i4, " iu1:", i4, " iu2:", i4)')&
      !  & io1, io2, iu1, iu2
      write(unitout, '("  lowest occupied state:", i4)') io1
      write(unitout, '("  highest occupied state:", i4)') io2
      write(unitout, '("  lowest unoccupied state:", i4)') iu1
      write(unitout, '("  highest unoccupied state:", i4)') iu2
      if(fensel) then
        write(unitout, '("Info(select_transitions): Selecting KS/QP transitions in&
          & the energy interval:")')
        write(unitout, '("  [",E10.3,",",E10.3,"]/H")') max(wl+econv(1),0.0d0), wu+econv(2)
        write(unitout, '("  [",E10.3,",",E10.3,"]/eV")')&
          & max(wl+econv(1),0.0d0)*h2ev, (wu+econv(2))*h2ev
        write(unitout, '("  Using convergence energy of:")')
        write(unitout, '("    ",2E10.3," /H")')&
          & max(econv(1), -wl), econv(2)
        write(unitout, '("    ",2E10.3," /eV")')&
          & max(econv(1), -wl)*h2ev, econv(2)*h2ev
      end if
      write(unitout, '("  Opening gap with a scissor of: ",&
        & F10.3,"/H ", F10.3,"/eV")'), sci, sci*h2ev

      !! Read in eigenvalues and occupancies for k-qmt/2 and k+qmt/2

      ! Set mod_kpoint / mod_Gkvector variables to the k-qmt/2-grid
      call init1offs(vkqmtmloff)
      ! Save the k-qmt/2-grid to modxs::vkl0 etc 
      call xssave0
      ! Set k and G+k variables in the standard locations mod_kpoint and mod_Gkvector 
      ! to those of the k'=k+qmt/2 grid.
      call init1offs(vkqmtploff)

      ! Check for coinciding k-grids 
      if(all(abs(vkqmtmloff-vkloff) < epslat)) then
        if (iqmt .ne. 1) then
          write(unitout, '("Info(select_transitions): Same k-grids for - and ref. grids at iqmt=",i3)') iqmt
        end if
        fsamek0=.true.
      else
        write(unitout, '("Info(select_transitions): Different k-grids for - and ref. grids at iqmt=",i3)') iqmt
        fsamek0=.false.
      end if
      if(all(abs(vkqmtploff-vkloff) < epslat)) then
        if (iqmt .ne. 1) then
          write(unitout, '("Info(select_transitions): Same k-grids for + and ref. grids at iqmt=",i3)') iqmt
        end if
        fsamek1=.true.
      else
        write(unitout, '("Info(select_transitions): Different k-grids for + and ref. grids at iqmt=",i3)') iqmt
        fsamek1=.false.
      end if
      if(all(abs(vkqmtmloff-vkqmtploff) < epslat)) then
        if (iqmt .ne. 1) then
          write(unitout, '("Info(select_transitions): Same k-grids for + and - grids at iqmt=",i3)') iqmt
        end if
        fsamek=.true.
      else
        write(unitout, '("Info(select_transitions): Different k-grids for + and - grids at iqmt=",i3)') iqmt
        fsamek=.false.
      end if

      ! Normal case: based on KS energies
      if(.not. associated(input%gw)) then 

        !! Get energies and occupancies for the k+qmt/2 grid
        ! Set EVALSV_QMTXYZ.OUT as read file
        if(fsamek1) then 
          call genfilname(iqmt=iqmtgamma, setfilext=.true.)
        else
          call genfilname(iqmt=iqmt, setfilext=.true.)
        end if
        do ik = 1, nkpt
          call getoccsv(vkl(1:3, ik), occsv(1:nstsv, ik))
          call getevalsv(vkl(1:3, ik), evalsv(1:nstsv, ik))
        end do

        !! Get energies and occupancies for k-qmt/2
        !! and save them in the 0 named variables
        ! Use modxs:filext0 in getoccsv0 and getevalsv0
        usefilext0 = .true.
        if(fsamek0) then
          ! Set EVALSV_QMT001.OUT as read file
          call genfilname(iqmt=iqmtgamma, fileext=filext0)
        else
          ! If + - grids are the same use the + one
          if(fsamek) then 
            ! Set EVALSV_QMTXYZ.OUT as read file
            call genfilname(iqmt=iqmt, fileext=filext0)
          ! Normal case: use the - grid 
          else
            ! Set EVALSV_QMTXYZ_m.OUT as read file
            call genfilname(iqmt=iqmt, auxtype="m", fileext=filext0)
          end if
        end if
        do ik = 1, nkpt
          call getoccsv0(vkl0(1:3, ik), occsv0(1:nstsv, ik))
          call getevalsv0(vkl0(1:3, ik), evalsv0(1:nstsv, ik))
        end do

        gwiomin = 1
        gwiumax = nstsv

      ! On top of GW
      else if(associated(input%gw) .and. iqmt==1) then 

        ! Get KS occupations/eigenvalues
        do ik = 1, nkpt
          call getoccsv(vkl(1:3, ik), occsv(1:nstsv, ik))
          call getevalsv(vkl(1:3, ik), evalsv(1:nstsv, ik))
        end do
        ! Read QP Fermi energies and eigenvalues from file
        ! only the eigenvalues are used in the following
        ! NOTE: QP evals are shifted by -efermi-eferqp with respect to KS evals
        ! NOTE: getevalqp sets mod_symmetry::nsymcrys to 1
        ! NOTE: getevalqp needs the KS eigenvalues as input
        nsymcrys_save = nsymcrys
        call getevalqp(nkptnr,vkl0,evalsv)
        nsymcrys = nsymcrys_save
        ! Set k and k'=k grid eigenvalues to QP energies
        evalsv0=evalsv
        occsv0=occsv

      else if(associated(input%gw) .and. iqmt /= 1) then 

        write(*,'("Error(b_bse): BSE+GW only supported for 0 momentum transfer.")')
        call terminate

      end if

      ! Sizes local/maximal
      if(fserial) then 
        nk_loc = nk_max
      else
        nk_loc = ceiling(real(nk_max,8)/real(mpiglobal%procs,8))
      end if
      hamsize_loc = nk_loc*nou_max

      ! The index mapping we want to build 
      ! s(1) = iuabs, s(2) = ioabs, s(3) = iknr
      allocate(smap_loc(3, hamsize_loc))

      ! Energy differences (local)
      allocate(de_loc(hamsize_loc))

      ! Flag map whether to use k-o-u combination or not
      allocate(sflag(hamsize_loc))
      sflag = .false.

      ! Occupation factor
      allocate(ofac_loc(hamsize_loc))

      ! How many o-u combinations at ik 
      if(allocated(kousize)) deallocate(kousize)
      allocate(kousize(nk_max))
      ! and limits of band ranges
      if(allocated(koulims)) deallocate(koulims)
      allocate(koulims(4,nk_max))
           
      ! Loop over kpoints (non-reduced)
      !   Distribute k-loop over global MPI communicator.
      !   Each participating rank gets a continuous ik interval.
      !   Not participating ranks have k1=0 and k2=-1.
      if(fserial) then 
        k1 = 1
        k2 = nk_max
      else
        k1 = firstofset(mpiglobal%rank, nk_max, mpiglobal%procs)
        k2 = lastofset(mpiglobal%rank, nk_max, mpiglobal%procs)
      end if

      ! Loop over reference k-grid
      ikloop: do ik = k1, k2

        ! Get k+qmt/2 index form k index
        ikqp = ik
        if(iqmt .ne. 0) ikqp = ik2ikqmtp(ik)
        ! Get k-qmt/2 index form k index
        ikqm = ik
        if(iqmt .ne. 0) ikqm = ik2ikqmtm(ik)

        kous = 0
        iomax = 0
        iumax = 0
        iomin = istocc0+1
        iumin = nstsv+1

        ! Loop over KS transition energies 
        !$OMP PARALLEL DO &
        !$OMP& COLLAPSE(2),&
        !$OMP& DEFAULT(SHARED), PRIVATE(io,iu,s,detmp,doctmp,posdiff),&
        !$OMP& REDUCTION(+:kous),&
        !$OMP& REDUCTION(max:iomax),&
        !$OMP& REDUCTION(min:iomin),&
        !$OMP& REDUCTION(max:iumax),&
        !$OMP& REDUCTION(min:iumin)
        do io = io1, io2
          do iu = iu1, iu2 

            if(fensel) then

              ! \Delta E = E_{u,k-qmt/2} - E_{o,k+qmt/2} + shift
              if (fxas) then
                detmp= evalsv0(iu, ikqm) - ecore(io) + sci
              else
                detmp = evalsv0(iu, ikqm) - evalsv(io, ikqp) + sci
              end if 

              ! Only consider transitions which are in the energy window
              if(detmp <= wu+econv(2) .and. detmp >= max(wl+econv(1),0.0d0)) then

                ! Only consider transitions which have a positve non-zero 
                ! occupancy difference f_{o ki+qmt/2} - f_{u ki-qmt/2}
                !
                ! If it is a XAS comutation:
                if(fxas) then 
                  doctmp = 1.0d0 - occsv0(iu, ikqm)
                ! If it is a optics comutation:
                else
                  doctmp = occsv(io, ikqp) - occsv0(iu, ikqm)
                end if
                ! Check if positive
                if(doctmp > cutoffocc) then
                  posdiff=.true.
                else 
                  posdiff=.false.
                end if

                if(posdiff) then 

                  ! Combine u, o and k index
                  ! u is counted from lumo=1 upwards
                  ! o is counted from lowest state upwards
                  ! ik index is shifted, due to MPI parallelization
                  s = hamidx(iu-istunocc0+1, io, ik-k1+1, nu_max, no_max)

                  ! Use that u-o-k combination
                  sflag(s) = .true.
                  
                  ! Write to combinded index map
                  smap_loc(1,s) = iu
                  smap_loc(2,s) = io
                  smap_loc(3,s) = ik

                  ! Save energy difference
                  de_loc(s) = detmp 

                  ! Save occupation factor
                  ofac_loc(s) = sqrt((occsv(io, ikqp) - occsv0(iu, ikqm))/maxocc)

                  ! Keep track of how many valid transitions
                  ! are considered at current k point.
                  kous = kous + 1

                  ! Keep track of the minimal/maximal io/iu
                  iomax = max(iomax, io)
                  iomin = min(iomin, io)
                  iumax = max(iumax, iu)
                  iumin = min(iumin, iu)

                end if

              end if

            ! Use user selected bands, but filter out transitions
            ! with non-positive occupation difference.
            else

              ! Only consider transitions which have a positve non-zero 
              ! occupancy difference f_{o ki+qmt/2} - f_{u ki-qmt/2}
              !
              ! If it is a XAS comutation:
              if(fxas) then 
                doctmp = 1.0d0 - occsv0(iu, ikqm)
              ! If it is a optics comutation:
              else
                doctmp = occsv(io, ikqp) - occsv0(iu, ikqm)
              end if
              ! Check if positive
              if(doctmp > cutoffocc) then
                posdiff=.true.
              else 
                posdiff=.false.
              end if

              if(posdiff) then 

                ! Combine u, o and k index
                ! u is counted from lumo=1 upwards
                ! o is counted from lowest state upwards
                ! ik index is shifted, due to MPI parallelization
                s = hamidx(iu-istunocc0+1, io, ik-k1+1, nu_max, no_max)

                ! Use that u-o-k combination
                sflag(s) = .true.
                
                ! Write to combinded index map
                smap_loc(1,s) = iu
                smap_loc(2,s) = io
                smap_loc(3,s) = ik

                ! Save energy difference
                ! \Delta E = E_{u,k-qmt/2} - E_{o,k+qmt/2} + shift
                if (fxas) then
                  de_loc(s)=evalsv0(iu,ikqm)-ecore(io)+sci
                else
                  de_loc(s) = evalsv0(iu,ikqm)-evalsv(io,ikqp)+sci
                end if

                ! Save occupation factor
                if (fxas) then
                  ofac_loc(s) = sqrt(1.0d0 - occsv0(iu, ikqm)/maxocc)
                else
                  ofac_loc(s) = sqrt((occsv(io, ikqp) - occsv0(iu, ikqm))/maxocc)
                end if

                ! Keep track of how many valid transitions
                ! are considered at current k point.
                kous = kous + 1

                ! Keep track of the minimal/maximal io/iu
                iomax = max(iomax, io)
                iomin = min(iomin, io)
                iumax = max(iumax, iu)
                iumin = min(iumin, iu)

              end if

            end if

          ! io iu loops
          end do
        end do
        !$OMP END PARALLEL DO

        ! u limits refer to the corresponding k_- 
        ! o limits refer to the corresponding k_+ 
        koulims(1,ik) = iumin
        koulims(2,ik) = iumax
        koulims(3,ik) = iomin
        koulims(4,ik) = iomax

        kousize(ik) = kous
        
      end do ikloop

#ifdef MPI
      if( .not. fserial) then
        ! Collect kousize on all processes 
        call mpi_allgatherv_ifc(set=nk_max, rlen=1, ibuf=kousize,&
          & inplace=.true., comm=mpiglobal)

        ! Collect koulims on all processes 
        call mpi_allgatherv_ifc(set=nk_max, rlen=4, ibuf=koulims,&
          & inplace=.true., comm=mpiglobal)
      end if
#endif
      ! Calculate maximal no(k_+) and nu(k_-)
      no_bse_max=0
      nu_bse_max=0
      do ik = 1, nk_max
        no_bse_max = max(koulims(4, ik)-koulims(3, ik)+1, no_bse_max)
        nu_bse_max = max(koulims(2, ik)-koulims(1, ik)+1, nu_bse_max)
      end do
      ! Maximal number of transitions at one k
      nou_bse_max=maxval(kousize)
      ! Global results
      ! Number of contributing k points
      nk_bse = count(kousize /= 0)
      nkkp_bse = nk_bse*(nk_bse+1)/2

      ! List of participating k-points
      if(allocated(kmap_bse_rg)) deallocate(kmap_bse_rg)
      allocate(kmap_bse_rg(nk_bse))
      if(allocated(kmap_bse_gr)) deallocate(kmap_bse_gr)
      allocate(kmap_bse_gr(nk_max))
      ik = 0
      do iknr = 1, nk_max
        ! If there is a valid transition from ik_- to ik_+,
        ! add the k-point to the bse-index. Keep maps
        ! between relative indexing (i.e. bse-k-point) and
        ! the global index of the full k-grid. 
        if(kousize(iknr) /= 0) then
          ik = ik + 1
          kmap_bse_rg(ik) = iknr
          kmap_bse_gr(iknr) = ik
        else
          kmap_bse_gr(iknr) = 0
        end if
      end do

      ! Size of the resulting resonant BSE Hamiltonian, i.e.
      ! number of transitions.
      hamsize = sum(kousize)
      ! Combined index map
      if(allocated(smap)) deallocate(smap)
      allocate(smap(3,hamsize))
      ! Occupation factors
      if(allocated(ofac)) deallocate(ofac)
      allocate(ofac(hamsize))
      ! Energy differences
      if(allocated(de)) deallocate(de)
      allocate(de(hamsize))

      ! If rank was participating in k-loop.
      if(k2 > 0) then 
        ! Apply selection flag to local arrays
        ! and store result in global counterparts.
        i1 = sum(kousize(1:k1-1))+1
        i2 = sum(kousize(k1:k2))+i1-1
        smap(1,i1:i2) = pack(smap_loc(1,:),sflag)
        smap(2,i1:i2) = pack(smap_loc(2,:),sflag)
        smap(3,i1:i2) = pack(smap_loc(3,:),sflag)
        ofac(i1:i2) = pack(ofac_loc,sflag)
        de(i1:i2) = pack(de_loc,sflag)
      end if

      ! Local auxiliary local arrays not needed anymore
      deallocate(smap_loc)
      deallocate(ofac_loc)
      deallocate(de_loc)
      deallocate(sflag)

#ifdef MPI
      if( .not. fserial) then 
        ! Collect ofac on all processes 
        call mpi_allgatherv_ifc(set=nk_max, rlenv=kousize, rbuf=ofac,&
          & inplace=.true., comm=mpiglobal)
        ! Collect de on all processes 
        call mpi_allgatherv_ifc(set=nk_max, rlenv=kousize, rbuf=de,&
          & inplace=.true., comm=mpiglobal)
        ! Collect smap on all processes 
        call mpi_allgatherv_ifc(set=nk_max, rlenv=kousize*3, ibuf=smap,&
          & inplace=.true., comm=mpiglobal)
      end if
#endif
      ! Energy sorting 
      if(allocated(ensortidx)) deallocate(ensortidx)
      allocate(ensortidx(hamsize))
      call sortidx(hamsize, de, ensortidx)
      
      ! Make relative combinded index map
      if(allocated(smap_rel)) deallocate(smap_rel)
      allocate(smap_rel(3,hamsize))
      do i1 = 1, hamsize
        iknr = smap(3, i1)
        smap_rel(1,i1) = smap(1, i1) - koulims(1, iknr) + 1
        smap_rel(2,i1) = smap(2, i1) - koulims(3, iknr) + 1
        smap_rel(3,i1) = kmap_bse_gr(iknr)
      end do

      write(unitout, '("Info(select_transitions):&
        & Number of participating transitions:", I8)') sum(kousize) 

      ! Print mappings to human readable file
      if(fserial) then 
        if(present(dirname)) then 
          call printso(iqmt, dirname)
        endif 
      else
        if(mpiglobal%rank == 0) then 
          if(present(dirname)) then 
            call printso(iqmt, dirname)
          end if 
        end if
      end if

      call timesec(t1)
      if (input%xs%BSE%outputlevelnumber == 1) then
        write(unitout, '("Info(select_transitions):&
          & Time needed/s:", f12.7)') t1-t0
      end if
      if( .not. fserial) then 
        call barrier(callername=trim(thisname))
      end if

    end subroutine select_transitions
    !EOC

    subroutine printso(iqmt, dirname)
      implicit none 

      integer(4) :: i, un, iqmt
      character(256) :: fdir, syscommand, fext, fname, fiqmt
      character(*), intent(in), optional :: dirname

      ! Make a folder 
      fdir = 'TRANSINFO'
      if(present(dirname)) then 
        fdir = trim(dirname)//'/'//trim(fdir)
      end if
      if(mpiglobal%rank == 0) then 
        syscommand = 'test ! -d '//trim(adjustl(fdir))&
         &//' && mkdir -p '//trim(adjustl(fdir))
        call system(trim(adjustl(syscommand)))
      end if
      write(fiqmt,*) iqmt
      fext = '_QMT'//trim(adjustl(fiqmt))//'.OUT'

      ! Write bse index map
      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'BSE_SINDEX'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# Combined BSE index @ Q =", 3(E10.3,1x))')&
        &  input%xs%qpointset%qpoint(:, iqmt)
      write(un,'("# s iu io ik iu_rel io_rel ik_rel occ")')
      do i = 1, size(ofac)
        write(un, '(7(I8,1x),1x,E23.16)')&
          & i, smap(:,i), smap_rel(:,i), ofac(i)
      end do
      close(un)

      ! Write ranges of occupied/unoccupied states involved at ik
      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'KOU'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# k-o-u ranges used in combined BSE index @ iqmt =", i3)')&
        & iqmt
      write(un,'("# u ranges refer to k_-, o ranges refer to k_+")')
      write(un,'("# ik ikm ikp iu1 iu2 io1 io2 nou")')
      do i = 1, nk_max
        write(un, '(8(I8,1x))')&
          & i, ik2ikqmtm(i), ik2ikqmtp(i),&
          & koulims(1,i), koulims(2,i), koulims(3,i), koulims(4,i), kousize(i)
      end do
      close(un)

      ! Write IP/QP transition energies
      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'EKSTRANS'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# KS transition energies associated with each combined index&
        & @ iqmt =", i3)') iqmt
      write(un,'("# s de de+sci")')
      do i = 1, hamsize
        write(un, '(I8,1x,E23.16,1x,E23.16)')&
          & i, de(i)-sci, de(i)
      end do
      close(un)
      ! Same but energy-sorted
      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'EKSTRANS_sorted'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# KS transition energies (eV) associated with each combined index&
        & @ iqmt =", i3)') iqmt
      write(un,'("# s de de+sci")')
      do i = 1, hamsize
        write(un, '(I8,1x,E23.16,1x,E23.16)')&
          & ensortidx(i), (de(ensortidx(i))-sci)*h2ev, de(ensortidx(i))*h2ev
      end do
      close(un)

    end subroutine printso

    ! Write out the coupling measures for each calculated exciton
    subroutine writemeasures(iqmt, nexc, evals, fcoup, dirname)
      use m_genfilname
      implicit none 

      integer(4), intent(in) :: iqmt, nexc
      logical, intent(in) :: fcoup
      real(8), intent(in) :: evals(:)
      character(*), intent(in), optional :: dirname

      character(256) :: fdir, syscommand, fext, fname
      integer(4) :: un

      real(8), allocatable :: measuresrr(:), measuresar(:)
      integer(4) :: alphamaxrr, alphamaxar, i, j
      integer(4), allocatable :: sorti(:)
      character(256) :: tdastring, bsetypestring, scrtypestring

      ! Make a folder 
      fdir = 'MEASURES'
      if(present(dirname)) then 
        fdir = trim(dirname)//'/'//trim(fdir)
      end if
      if(mpiglobal%rank == 0) then 
        syscommand = 'test ! -d '//trim(adjustl(fdir))//' && mkdir -p '//trim(adjustl(fdir))
        call system(trim(adjustl(syscommand)))
      end if

      ! Generate file name
      if(input%xs%bse%coupling) then
        tdastring=''
      else
        if(input%xs%bse%chibarq) then 
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if
      if(input%xs%bse%bsetype == "IP") then
        tdastring=''
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      ! Make measures
      allocate(measuresrr(hamsize))
      allocate(measuresar(hamsize))
      allocate(sorti(hamsize))

      measuresrr = 0.0d0
      alphamaxrr = 1
      measuresrr = vwdiffrr(1:hamsize)/de(1:hamsize)
      alphamaxrr = maxloc(measuresrr,1)

      measuresar = 0.0d0
      alphamaxar = 1
      if(fcoup) then
        measuresar = vwdiffar(1:hamsize)/de(1:hamsize)
        alphamaxar = maxloc(measuresar,1)
      end if

      call getunit(un)

      call genfilname(basename='Coupling_Measures', iqmt=iqmt,&
        & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
        & nar= .not. input%xs%bse%aresbse, filnam=fname, dirname=trim(fdir))

      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# Measures for excitions @ Q =", 3(E10.3,1x))')  input%xs%qpointset%qpoint(:, iqmt)
      write(un,'("#")')
      write(un,'("# RR: Max_{a,b} |V_ab - W^rr_ab|/dE^ip_a")')
      write(un,'("# a=",i8," dE=",E13.4)') alphamaxrr, de(alphamaxrr)*h2ev
      if(fcoup) then
        write(un,'("# AR: Max_{a,b} |V_ab - W^ar_ab|/dE^ip_a")')
        write(un,'("# a=",i8," dE=",E13.4)') alphamaxar, de(alphamaxar)*h2ev
      end if
      write(un,'("# lambda, exenrgy, MaxMeasure_rr, MaxMeasure_ar")')

      do i = 1, nexc

        measuresrr = 0.0d0
        measuresrr = vwdiffrr(1:hamsize)/(de(1:hamsize)+evals(i))

        measuresar = 0.0d0
        if(fcoup) then
          measuresar = vwdiffar(1:hamsize)/(de(1:hamsize)+evals(i))
        end if

        write(un, '(I8,1x,E13.4,1x,2(E13.4,1x))')&
          & i, evals(i)*h2ev,&
          & measuresrr(alphamaxrr),&
          & measuresar(alphamaxar)

      end do

      close(un)

      call getunit(un)

      call genfilname(basename='VW_diff', iqmt=iqmt,&
        & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
        & nar= .not. input%xs%bse%aresbse, filnam=fname, dirname=trim(fdir))

      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# max per row of |V-W| @ Q =", 3(E10.3,1x))')  input%xs%qpointset%qpoint(:, iqmt)
      write(un,'("# alpha, ipen, VWdiff_rr, VWdiff_ra")')

      call sortidx(hamsize, de, sorti)

      do i = 1, hamsize

        j = sorti(i)
        if(fcoup) then 
          write(un, '(I8,1x,E13.4,1x,2(E13.4,1x))')&
            & j, de(j)*h2ev,&
            & vwdiffrr(j), vwdiffar(j)
        else
          write(un, '(I8,1x,E13.4,1x,2(E13.4,1x))')&
            & j, de(j)*h2ev,&
            & vwdiffrr(j), 0.0d0
        end if

      end do

      close(un)

    end subroutine writemeasures

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Simple continuous combined index mappers  !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: hamidx
    ! !INTERFACE:
    integer(4) function hamidx(i1, i2, ik, n1, n2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: i1, i2, ik ! Indices counting from 1 continuously 
    ! integer(4) :: n1, n2     ! Maximum values of i1, i2 respectively
    ! Out:
    ! integer(4) :: hamidx     ! Combined index
    !
    ! !DESCRIPTION:
    !   The function returns a combined index given two band indices
    !   and a $\vec{k}$ index. It is use in the construction of
    !   the BSE Hamiltonian.\\
    !   Map:\\
    !   $\text{hamdix} = i_1 + (i_2 - 1) n_1 + n_1 n_2 (i_k-1)$\\
    !   Notes:\\
    !     $i_1$ is the fastest varying index, followed in order by $i_2$ and $i_k$.\\
    !     All indices are assumed to be counted from 1 onwards continuously.
    ! 
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Changed fastes index to i1. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: i1, i2, ik, n1, n2
      hamidx = i1 + n1 * (i2-1) + n1 * n2 * (ik-1)
    end function hamidx
    !EOC

    !BOP
    ! !ROUTINE: hamidx_back
    ! !INTERFACE:
    subroutine hamidx_back(s, i1, i2, ik, n1, n2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: s          ! Combined index created with {\tt hamidx}
    ! integer(4) :: n1, n2     ! Maximum values of i1, i2 respectively
    ! Out:
    ! integer(4) :: i1, i2, ik ! Individual indices 
    !
    ! !DESCRIPTION:
    !   The subroutine does the inverse operation of the function {\tt hamidx}.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: s, n1, n2
      integer(4), intent(out) :: i1, i2, ik
      integer(4) :: n12, tmp
      n12 = n1*n2
      ik = (s-1)/n12 + 1
      tmp = s - (ik-1)*n12
      i2 = (tmp-1)/n1 + 1
      i1 = tmp - (i2-1)*n1
    end subroutine hamidx_back
    !EOC

    !BOP
    ! !ROUTINE: subhamidx
    ! !INTERFACE:
    integer(4) function subhamidx(i1, i2, n1)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: i1, i2     ! Indices counting from 1 continuously 
    ! integer(4) :: n1         ! Maximum value of i1 
    ! Out:
    ! integer(4) :: subhamidx  ! Combined index
    !
    ! !DESCRIPTION:
    !   The function return a combined index given two indices.
    !   It is use in the construction of the BSE Hamiltonian.\\
    !   Map:\\
    !   $\text{hamdix} = i_1 + (i_2 - 1) n_1$\\
    !   Notes:\\
    !     $i_1$ is the fastest varying index, followed by $i_2$.\\
    !     All indices are assumed to be counted from 1 onwards continuously.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !   Changed fastest index to i1. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: i1, i2, n1
      subhamidx = i1 + n1 * (i2-1)
    end function subhamidx
    !EOC

    !BOP
    ! !ROUTINE: subhamidx_back
    ! !INTERFACE:
    subroutine subhamidx_back(s, i1, i2, n1)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: s    ! Combined index
    ! integer(4) :: n1   ! Maximum value of i1 
    ! Out:
    ! integer(4) :: i1, i2   ! Individual indices
    !
    ! !DESCRIPTION:
    !   The routine performs the inverse operation to 
    !   {\tt subhamidx}.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !   Changed fastest index to i1. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: s, n1
      integer(4), intent(out) :: i1, i2
      i2 = (s-1)/n1 + 1
      i1 = s - (i2-1)*n1
    end subroutine subhamidx_back
    !EOC

end module modbse
!EOC
