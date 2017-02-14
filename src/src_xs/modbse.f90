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
  use modinput, only: input
  use mod_constants, only: h2ev
  use modxs, only: unitout
  use modxs, only: evalsv0, occsv0, ikmapikq
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
  ! KS energy differences
  real(8), allocatable :: de(:)
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
  integer(4), dimension(:), allocatable :: ik2ikqmtp, ik2ikqmtm, imk2imkqmtp
  ! Offsets of k grids
  real(8), dimension(3) :: vkloff, vkqmtploff, vkqmtmloff

  ! Filebasenames
  character(256) :: infofbasename = "BSEINFO"
  character(256) :: scclifbasename = "SCCLI"
  character(256) :: scclicfbasename = "SCCLIC"
  character(256) :: scclictifbasename = "SCCLICTI"
  character(256) :: exclifbasename = "EXCLI"
  character(256) :: exclicfbasename = "EXCLIC"
  character(256) :: infofname
  character(256) :: scclifname
  character(256) :: exclifname

  ! Legacy
  ! GW eigenvalue backup 
  real(8), allocatable, dimension(:,:) :: eval0

  contains

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Routines to setup the combinded index of  ! 
    ! the BSE hamiltonian.                      !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: setranges_modxs
    ! !INTERFACE:
    subroutine setranges_modxs(iqmt, fcoup, fti)
      use modinput
      use mod_misc, only: filext
      use mod_kpoint, only: nkpt, vkl
      use mod_xsgrids
      use mod_Gkvector, only: gkmax
      use modxs, only: vqlmt, evalsv0, usefilext0, filext0,&
                     & ksgap, ksgapval, qmtpgap, qmtmgap, unitout, vkl0
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
      logical, intent(in) :: fcoup, fti

      integer(4) :: iomax, iumin

      integer(4) :: iomax_kkqmtp, iumin_kkqmtp
      integer(4) :: iomax_kkqmtm, iumin_kkqmtm
      integer(4) :: iomax_mkmkqmtp, iumin_mkmkqmtp

      integer(4), dimension(:), allocatable :: io_k, iu_k
      integer(4), dimension(:), allocatable :: io_kqmtp, iu_kqmtp
      integer(4), dimension(:), allocatable :: io_kqmtm, iu_kqmtm 
      integer(4), dimension(:), allocatable :: io_mk, iu_mk
      integer(4), dimension(:), allocatable :: io_mkqmtp, iu_mkqmtp

      integer(4) :: ik

      logical :: fgap
      real(8) :: gap

      !write(*,*) "setranges_modxs here!"
      !write(*,*) "iqmt=", iqmt
      !write(*,'(a,3E10.3)') "vqlmt=", vqlmt(1:3,iqmt)

      !---------------------------------------------------!
      ! Get offsets of and mapping between k and k' grids !
      !---------------------------------------------------!
      ! Note: requires init2 to set up vqlmt
      call xsgrids_init(vqlmt(1:3,iqmt), gkmax)

      !! Offsets for k and k' grids
      ! k
      vkloff = k_kqmtp%kset%vkloff
      ! k+qmt
      vkqmtploff = k_kqmtp%kqmtset%vkloff
      ! k-qmt
      vkqmtmloff = k_kqmtm%kqmtset%vkloff

      !! Mappings between k and k' grids
      ! k --> k+qmt
      if(allocated(ik2ikqmtp)) deallocate(ik2ikqmtp)
      allocate(ik2ikqmtp(nkpt))
      ik2ikqmtp(:) = k_kqmtp%ik2ikqmt(:)
      ! k --> k-qmt
      if(allocated(ik2ikqmtm)) deallocate(ik2ikqmtm)
      allocate(ik2ikqmtm(nkpt))
      ik2ikqmtm(:) = k_kqmtm%ik2ikqmt(:)

      call xsgrids_finalize()
      !---------------------------------------------------!

      !write(*,*)
      !write(*,'(a,3E10.3)') "vkloff = ", vkloff
      !write(*,'(a,3E10.3)') "vkqmtploff = ", vkqmtploff
      !write(*,'(a,3E10.3)') "vkqmtmloff = ", vkqmtmloff
      !write(*,*)
      !write(*,*) "ik2ikqmtp mapping"
      do ik = 1, nkpt
        !write(*,'(2i3)') ik, ik2ikqmtp(ik)
      end do
      !write(*,*)
      !write(*,*) "ik2ikqmtm mapping"
      do ik = 1, nkpt
        !write(*,'(2i3)') ik, ik2ikqmtm(ik)
      end do

      !---------------------------------------------------!
      ! Inspect occupancies for k, k'=k+qmt               !
      !---------------------------------------------------!
      ! Reset mod_kpoint / mod_Gkvector variables to the unshifted k-grid
      ! (apart from xs%vkloff)
      call init1offs(vkloff)

      !write(*,*)
      !write(*,*) "ik vkl"
      do ik = 1, nkpt
        !write(*,'(i3, 3E10.3)') ik, vkl(1:3,ik)
      end do

      ! Save the k-grid to modxs vkl0 etc 
      call xssave0

      !write(*,*)
      !write(*,*) "ik vkl0"
      do ik = 1, nkpt
        !write(*,'(i3, 3E10.3)') ik, vkl0(1:3,ik)
      end do

      ! Set k and G+k variables in the standard locations mod_kpoint and mod_Gkvector 
      ! to those of the k' grid.
      call init1offs(vkqmtploff)

      !write(*,*)
      !write(*,*) "ik vkl"
      do ik = 1, nkpt
        !write(*,'(i3, 3E10.3)') ik, vkl(1:3,ik)
      end do

      ! Allocate k eigenvalues
      ! Note: If evalsv0 is allocated before findocclims is called it gets read in 
      ! and stays allocated.
      if(.not. allocated(evalsv0)) allocate(evalsv0(nstsv, nkpt))

      ! Explicitly specify k file extension
      usefilext0 = .true.
      ! Set EVALSV_QMT001.OUT as first reference for the occupation search
      call genfilname(iqmt=iqmtgamma, setfilext=.true.)
      filext0 = filext

      !write(*,*) "filext0 =", trim(filext0)

      ! Set k' file extension
      ! Set EVALSV_QMTXYZ.OUT as second reference for the occupation search
      call genfilname(iqmt=iqmt, setfilext=.true.)

      !write(*,*) "filext =", trim(filext)

      allocate(io_k(nkpt), iu_k(nkpt))
      allocate(io_kqmtp(nkpt), iu_kqmtp(nkpt))

      ! Inspect occupation limits for k and k+qmt grids
      call findocclims(iqmt, ik2ikqmtp, iomax_kkqmtp, iumin_kkqmtp,&
        & io_k, io_kqmtp, iu_k, iu_kqmtp)

      ! Save gap status
      !   Was an indirect gap found?
      fgap = ksgap
      !   Size of found gap
      gap = ksgapval
      !   Size of qmt dependent gap (with qmt=0 this is the direct gap)
      qmtpgap = qgap

      !!write(*,*) "fgap =", fgap
      !write(*,*) "ksgap =", gap
      !write(*,*) "qmtpgap =", qgap

      !write(*,*) "ksgapval", ksgapval
      !write(*,*) "max(io_k)", maxval(io_k)
      !write(*,*) "max(io_kqmtp)", maxval(io_kqmtp)
      !write(*,*) "min(iu_k)", minval(iu_k)
      !write(*,*) "min(iu_kqmtp)", maxval(iu_kqmtp)

      deallocate(io_k, iu_k)
      deallocate(io_kqmtp, iu_kqmtp)
      !---------------------------------------------------!

      !---------------------------------------------------!
      ! When using couping terms in the BSE               !
      ! (in the std ar basis) also k, k'=k-qmt occupation !
      ! limits need to be inspected.                      !
      !---------------------------------------------------!
      if(fcoup .and. .not. fti) then 

        !write(*,*) "COUPLING k k-qmt"

        !---------------------------------------------------!
        ! Inspect occupancies for k, k'=k-qmt               !
        !---------------------------------------------------!
        ! Reset mod_kpoint / mod_Gkvector variables to the unshifted k-grid
        ! (apart from xs%vkloff)
        call init1offs(vkloff)

        !write(*,*)
        !write(*,*) "ik vkl"
        do ik = 1, nkpt
          !write(*,'(i3, 3E10.3)') ik, vkl(1:3,ik)
        end do

        ! Save the k-grid to modxs vkl0 etc 
        call xssave0

        !write(*,*)
        !write(*,*) "ik vkl0"
        do ik = 1, nkpt
          !write(*,'(i3, 3E10.3)') ik, vkl0(1:3,ik)
        end do

        ! Set k and G+k variables in the standard locations mod_kpoint and mod_Gkvector 
        ! to those of the k' grid.
        call init1offs(vkqmtmloff)

        !write(*,*)
        !write(*,*) "ik vkl"
        do ik = 1, nkpt
          !write(*,'(i3, 3E10.3)') ik, vkl(1:3,ik)
        end do

        ! Explicitly specify k file extension
        usefilext0 = .true.
        ! Set EVALSV_QMT001.OUT as first reference for the occupation search
        call genfilname(iqmt=iqmtgamma, setfilext=.true.)
        filext0 = filext

        !write(*,*) "filext0 =", trim(filext0)

        ! Set k' file extension
        ! Set EVALSV_QMTXYZ_mqmt.OUT as second reference for the occupation limits search
        if(iqmt /= 1) then 
          call genfilname(iqmt=iqmt, auxtype="mqmt", setfilext=.true.)
        else
          call genfilname(iqmt=iqmt, setfilext=.true.)
        end if

        !write(*,*) "filext =", trim(filext)

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

        !!write(*,*) "fgap =", fgap
        !write(*,*) "ksgap =", gap
        !write(*,*) "qmtmgap =", qgap

        !write(*,*) "ksgapval", ksgapval
        !write(*,*) "max(io_k)", maxval(io_k)
        !write(*,*) "max(io_kqmtm)", maxval(io_kqmtm)
        !write(*,*) "min(iu_k)", minval(iu_k)
        !write(*,*) "min(iu_kqmtm)", maxval(iu_kqmtm)

        ! Use the highest partially occupied band over all k, k+qmt, k-qmt
        iomax = max(iomax_kkqmtp, iomax_kkqmtm)
        ! Use the lowest partially unoccupied band over all k, k+qmt, k-qmt
        iumin = min(iumin_kkqmtp, iumin_kkqmtm)

        deallocate(io_k, iu_k)
        deallocate(io_kqmtm, iu_kqmtm)

      else

        ! Use the highest partially occupied band over all k, k+qmt
        iomax = iomax_kkqmtp
        ! Use the lowest partially unoccupied band over all k, k+qmt
        iumin = iumin_kkqmtp
        
      end if

      !write(*,*) "iomax=", iomax
      !write(*,*) "iumin=", iumin

      ! Set gap status in modxs
      ksgap = fgap
      ksgapval = gap

      ! Do I need to set modxs: 
      !   isto, isto0, istu0, istu, istunocc0, istunocc, istocc0, istocc
      ! ?
      !call findocclims(iqmt, istocc0, istocc, istunocc0,&
      !  & istunocc, isto0, isto, istu0, istu)
      
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

      if(ksgapval == 0.0d0) then
        write(unitout, '("Warning (setranges_modxs): The system has no gap")')
        if(sci /= 0.0d0) then 
          write(unitout, '("Warning (setranges_modxs):&
            &   Scissor > 0 but no gap. Setting scissor to 0.")')
          sci = 0.0d0
        end if
      end if  

    end subroutine setranges_modxs
    !EOC

    !BOP
    ! !ROUTINE: select_transitions
    ! !INTERFACE:
    subroutine select_transitions(iqmt, serial)
      use mod_kpoint, only: vkl
      use mod_misc, only: filext
      use modxs, only: usefilext0, filext0, vkl0
      use m_genfilname
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

      logical :: fserial
      integer(4) :: ik, ikq, s, iknr
      integer(4) :: io, iu, kous
      integer(4) :: io1, io2, iu1, iu2
      integer(4) :: iomax, iumax
      integer(4) :: iomin, iumin
      real(8), parameter :: maxocc = 2.0d0

      integer(4) :: nk_loc, hamsize_loc
      integer(4), allocatable :: smap_loc(:,:)
      real(8) :: detmp
      real(8), allocatable :: ofac_loc(:)
      real(8), allocatable :: de_loc(:)
      logical, allocatable :: sflag(:)
      integer(4) :: buflen, buflen_r
      integer(4) :: k1, k2, k1_r, k2_r
      integer(4) :: iproc, i1, i2, il_r, iu_r

      !write(*,*) "select_transitions here with iqmt=", iqmt

      if(present(serial)) then 
        fserial = serial
      else
        fserial = .false.
      end if

      ! Search for needed KS transitions automatically
      ! depending on the chosen energy window?
      if(any(input%xs%bse%nstlbse == 0)) then
        fensel = .true.
      else
        fensel = .false.
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
        io1 = 1
        io2 = istocc0
        iu1 = istunocc0
        iu2 = nstsv
      else
        io1 = input%xs%bse%nstlbse(1)
        io2 = input%xs%bse%nstlbse(2)
        iu1 = input%xs%bse%nstlbse(3)+istunocc0-1
        iu2 = input%xs%bse%nstlbse(4)+istunocc0-1
      end if

      if(mpiglobal%rank == 0) then 
        if(fensel) then
          write(unitout, '("Info(select_transitions): Searching for KS transitions in&
            & the energy interval:")')
          write(unitout, '("  [",E10.3,",",E10.3,"]/H")') max(wl+econv(1),0.0d0), wu+econv(2)
          write(unitout, '("  [",E10.3,",",E10.3,"]/eV")')&
            & max(wl+econv(1),0.0d0)*h2ev, (wu+econv(2))*h2ev
          write(unitout, '("  Using convergence energy of:")')
          write(unitout, '("    ",2E10.3," /H")')&
            & max(econv(1), -wl), econv(2)
          write(unitout, '("    ",2E10.3," /eV")')&
            & max(econv(1), -wl)*h2ev, econv(2)*h2ev
        else
          write(unitout, '("Info(select_transitions): Searching for KS transitions in&
            & the band interval:")')
          write(unitout, '("  io1:", i4, " io2:", i4, " iu1:", i4, " iu2:", i4)')&
            & io1, io2, iu1, iu2
        end if
        write(unitout, '("  Opening gap with a scissor of:",&
          & E10.3,"/H", E10.3,"/eV")'), sci, sci*h2ev
      end if

      !! Read in eigenvalues and occupancies for k and k+qmt
      !write(*,*) "reading in eigenvalues and occupancies"

      ! Reset mod_kpoint / mod_Gkvector variables to the unshifted k-grid
      ! (apart from xs%vkloff)
      call init1offs(vkloff)
      ! Save the k-grid to modxs vkl0 etc 
      call xssave0
      ! Set k and G+k variables in the standard locations mod_kpoint and mod_Gkvector 
      ! to those of the k'=k+qmt grid.
      call init1offs(vkqmtploff)
      ! Explicitly specify k file extension
      usefilext0 = .true.
      ! Set EVALSV_QMT001.OUT as first reference for the occupation search
      call genfilname(iqmt=iqmtgamma, setfilext=.true.)
      filext0 = filext
      !write(*,*) "filext0 =", trim(filext0)
      do ik = 1, nkpt
        call getoccsv0(vkl0(1:3, ik), occsv0(1:nstsv, ik))
        call getevalsv0(vkl0(1:3, ik), evalsv0(1:nstsv, ik))
      end do
      ! Set k' file extension
      ! Set EVALSV_QMTXYZ.OUT as second reference for the occupation search
      call genfilname(iqmt=iqmt, setfilext=.true.)
      !write(*,*) "filext =", trim(filext)
      do ik = 1, nkpt
        call getoccsv(vkl(1:3, ik), occsv(1:nstsv, ik))
        call getevalsv(vkl(1:3, ik), evalsv(1:nstsv, ik))
      end do

      ! Sizes local/maximal
      if(fserial) then 
        nk_loc = nk_max
        hamsize_loc = nk_loc*nou_max
      else
        nk_loc = ceiling(real(nk_max,8)/real(mpiglobal%procs,8))
        hamsize_loc = nk_loc*nou_max
      end if

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

      ikloop: do ik = k1, k2

        ! Get k+q index form k index
        ikq = ik
        if(iqmt .ne. 0) ikq = ik2ikqmtp(ik)

        kous = 0
        iomax = 0
        iumax = 0
        iomin = istocc0+1
        iumin = nstsv+1

        ! Loop over KS transition energies 
        !$OMP PARALLEL DO &
        !$OMP& COLLAPSE(2),&
        !$OMP& DEFAULT(SHARED), PRIVATE(io,iu,s,detmp),&
        !$OMP& REDUCTION(+:kous),&
        !$OMP& REDUCTION(max:iomax),&
        !$OMP& REDUCTION(min:iomin),&
        !$OMP& REDUCTION(max:iumax),&
        !$OMP& REDUCTION(min:iumin)
        do io = io1, io2
          do iu = iu1, iu2 

            if(fensel) then

              detmp = evalsv(iu, ikq) - evalsv0(io, ik) + sci 

              ! Only consider transitions which are in the energy window
              ! \Delta E = \epsilon_{u ki+q} - \epsilon_{o ki}
              if(detmp <= wu+econv(2) .and. detmp >= max(wl+econv(1),0.0d0)) then

                ! Only consider transitions which have a positve non-zero 
                ! occupancy difference f_{o ki} - f_{u ki+q}
                if( occsv0(io, ik) - occsv(iu, ikq) > cutoffocc) then 

                  ! Combine u, o and k index
                  ! u is counted from lumo=1 upwards
                  ! o is counted from lowest state upwards
                  ! ik index is shifted, due to MPI parallelization
                  s = hamidx(iu-istunocc0+1, io, ik-k1+1, nu_max, no_max)

                  ! Use that u-o-k combination
                  sflag(s) = .true.
                  
                  ! Write to combinded index map
                  smap_loc(:,s) = [iu, io, ik]

                  ! Save energy difference
                  de_loc(s) = detmp 

                  ! Save occupation factor
                  ofac_loc(s) = sqrt((occsv0(io, ik) - occsv(iu, ikq))/maxocc)

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
              ! occupancy difference f_{o ki} - f_{u ki+q}
              if( occsv0(io, ik) - occsv(iu, ikq) > cutoffocc) then 

                ! Combine u, o and k index
                ! u is counted from lumo=1 upwards
                ! o is counted from lowest state upwards
                ! ik index is shifted, due to MPI parallelization
                s = hamidx(iu-istunocc0+1, io, ik-k1+1, nu_max, no_max)

                ! Use that u-o-k combination
                sflag(s) = .true.
                
                ! Write to combinded index map
                smap_loc(:,s) = [iu, io, ik]

                ! Save energy difference
                de_loc(s) = evalsv(iu,ikq)-evalsv0(io,ik)+sci

                ! Save occupation factor
                ofac_loc(s) = sqrt((occsv0(io, ik) - occsv(iu, ikq))/maxocc)

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

          end do
        end do
        !$OMP END PARALLEL DO

        koulims(:,ik) = [ iumin, iumax, iomin, iomax ] 

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
      ! Calculate maximal no(k) and nu(k)
      no_bse_max=0
      nu_bse_max=0
      do ik = 1, nk_max
        no_bse_max = max(koulims(4, ik)-koulims(3, ik)+1, no_bse_max)
        nu_bse_max = max(koulims(2, ik)-koulims(1, ik)+1, nu_bse_max)
      end do
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
        if(kousize(iknr) /= 0) then
          ik = ik + 1
          kmap_bse_rg(ik) = iknr
          kmap_bse_gr(iknr) = ik
        else
          kmap_bse_gr(iknr) = 0
        end if
      end do

      ! Size of the resulting resonant BSE Hamiltonian
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

      if(mpiglobal%rank == 0) then 
        write(unitout, '("Info(select_transitions):&
          & Number of participating transitions:", I8)') sum(kousize) 
      end if

      if(mpiglobal%rank == 0) then 
        call printso(iqmt)
      end if

    end subroutine select_transitions
    !EOC

    subroutine printso(iqmt)
      implicit none 

      integer(4) :: i, un, iqmt
      character(256) :: fdir, syscommand, fext, fname, fiqmt

      ! Make a folder 
      fdir = 'TRANSINFO'
      if(rank == 0) then 
        syscommand = '[[ ! -e '//trim(adjustl(fdir))//' ]] && mkdir '//trim(adjustl(fdir))
        call system(trim(adjustl(syscommand)))
      end if
      write(fiqmt,*) iqmt
      fext = '_QMT'//trim(adjustl(fiqmt))//'.OUT'

      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'BSE_SINDEX'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# Combined BSE index @ q =", 3(E10.3,1x))')  0.0d0, 0.0d0, 0.0d0
      write(un,'("# s iu io ik iu_rel io_rel ik_rel occ")')
      do i = 1, size(ofac)
        write(un, '(7(I8,1x),1x,E23.16)')&
          & i, smap(:,i), smap_rel(:,i), ofac(i)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'KOU'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# k-o-u ranges used in combined BSE index @ iqmt =", i3)')&
        & iqmt
      write(un,'("# ik iu1 iu2 io1 io2 nou")')
      do i = 1, nk_max
        write(un, '(6(I8,1x))')&
          & i, koulims(1,i), koulims(2,i), koulims(3,i), koulims(4,i), kousize(i)
      end do
      close(un)

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

      call getunit(un)
      fname = trim(adjustl(fdir))//'/'//'EKSTRANS_sorted'//fext
      open(un, file=trim(adjustl(fname)), action='write', status='replace')
      write(un,'("# KS transition energies associated with each combined index&
        & @ iqmt =", i3)') iqmt
      write(un,'("# s de de+sci")')
      do i = 1, hamsize
        write(un, '(I8,1x,E23.16,1x,E23.16)')&
          & ensortidx(i), de(ensortidx(i))-sci, de(ensortidx(i))
      end do
      close(un)

    end subroutine printso

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
