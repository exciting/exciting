!BOP
! !ROUTINE: bselauncher
! !INTERFACE:
subroutine bselauncher
! !USES:
  use modmpi
  use modscl
  use modxs, only: unitout
  use modinput, only: input, input_type
  use modbse, only: select_bse_solver, setup_bse_type_list
  use bsemain
  use os_utils, only: make_directory

! !DESCRIPTION:
!   Launches the construction and solving of the Bethe-Salpeter Hamiltonian
!   for the specified $\vec{Q}_\text{mt}$ momentum transfer and approximation
!   (TDA or non-TDA).
!
! !REVISION HISTORY:
!   Created. 2016 (Aurich)
!EOP
!BOC

  implicit none

  ! Local vars
  character(*), parameter :: thisname = "bselauncher"
  integer, parameter :: num_bse_dirs = 3
  integer, parameter :: max_dir_len = 12
  !
  logical :: fdist, fcoup, fchibarq
  character(len=max_dir_len) :: bse_dir_list(num_bse_dirs)
  character(256) :: syscommand, epsilondir, lossdir, sigmadir
  integer(4) :: idir, iqmt, iqmti, iqmtf, nqmt, nqmtselected, iq1, iq2
  real(8) :: ts0, ts1
  real(8) :: vqmt(3)
  integer :: bse_type_index
  character(len=256), allocatable :: bsetypelist(:)
  ! Process grid shape
  integer(4) :: nprows, npcols
  ! Number of used processes to be used for blacs
  integer(4) :: nprocs, nprocs2d
  ! mpi instance only for active ranks, used for ELPA
  type(mpiinfo) :: mpicom_active
  ! MPI communicator of the active BLACS grid
  integer(4) :: comm_active
  ! distinguish idle ranks and error status
  integer(4) :: color, ierror
  
  !---------------------------------------------------------------------------!
  ! Init0,1,2 General inits
  !---------------------------------------------------------------------------!
  ! Start timer for init calls
  call timesec(ts0)
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
  !     (modxs::totalqmtl etc. and mod_qpoint::vql)
  !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
  !   * non-reduced mapping between (ik,qmt) and ik' grids (modxs::ikmapikq)
  !   * G+qmt quantities (modxs)
  !   * The square root of the Coulomb potential for the G+qmt points
  !   * Reads STATE.OUT
  !   * Generates radial functions (mod_APW_LO)
  call init2
  
  ! xas and xes specific init (has to come after init0 and init1)
  if(input%xs%bse%xas .or. input%xs%BSE%xes) call xasinit
  
  ! End timer for init calls
  call timesec(ts1)
  write(unitout, '("Info(",a,"):&
    & Init time:", f12.6)') trim(thisname), ts1 - ts0
  !---------------------------------------------------------------------------!
  ! Create Output directories                                                 !
  !---------------------------------------------------------------------------!
  bse_dir_list= (/ 'EPSILON','LOSS   ','SIGMA  ' /)
  do idir = 1, num_bse_dirs
    call make_directory(bse_dir_list(idir), mpiglobal)
  end do
  !---------------------------------------------------------------------------!
  ! Check Q-point sublist range
  !---------------------------------------------------------------------------!
  ! Use all
  nqmt = size(input%xs%qpointset%qpoint, 2)
  iqmti = 1
  iqmtf = nqmt
  !   or use range (iqmtrange defaults to "1 1")
  if(input%xs%bse%iqmtrange(1) /= -1) then 
    iqmti=input%xs%bse%iqmtrange(1)
    iqmtf=input%xs%bse%iqmtrange(2)
  end if
  nqmtselected = iqmtf-iqmti+1
  ! Check requested range is compatible with qpointlist
  if(iqmtf > nqmt .or. iqmti < -1 .or. iqmti > iqmtf) then 
    write(unitout, '("Error(",a,"):", a)') trim(thisname),&
      & " iqmtrange incompatible with qpointset list"
    call terminate
  end if

  ! Info out
  call printline(unitout, "+")
  write(unitout, '("Info(",a,"):", a)') trim(thisname),&
    & " Setting up and diagonalizing BSE Hamiltonian."
  write(unitout, '("Info(",a,"):", a, i3, a, i3)') trim(thisname),&
    & " Using momentum transfer vectors from list : ", iqmti, " to", iqmtf
  call printline(unitout, "+")
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  ! Set up process grids for BLACS (if compiled with DSCAL, dummies otherwise)
  ! DSCAL implies DMPI
  !---------------------------------------------------------------------------!
  ! overwrite the input parameter, distribute if ScaLAPACK or ELPA are requested
  call select_bse_solver(input)
  fdist = input%xs%bse%bsesolver /= 'lapack'

  !   Make square'ish process grid (context 0)
#ifdef _ELPA_
  ! Make square'ish porcess grid out of 'nprocs' processes to guarantee that only active ranks participate
  nprocs = mpiglobal%procs
  npcols = int(sqrt(dble(nprocs)))
  nprows = nprocs / npcols
  nprocs2d = npcols*nprows
  if(nprocs2d /= nprocs) then
    write(unitout,'("Info(setupblacs):&
      & Warning - Processes do not fit 2d grid")')
    write(unitout,'("Info(setup2dblacs):&
      & Warning - ",i2," processes not used")') nprocs-nprocs2d
  end if
  if(mpiglobal%rank < nprocs2d) then
    color = 1
  else
    color = MPI_UNDEFINED
  end if
  call MPI_Comm_split(mpiglobal%comm, color, 0, comm_active, ierror)
  if (color /= MPI_UNDEFINED) then
    call mpicom_active%init(comm_active)
    call setupblacs(mpicom_active, 'grid', bi2d, np=nprocs2d)
  else
    call setupblacs(mpiglobal, 'xxx', bi2d, np=nprocs2d)
  end if
#else
  ! legacy call for ScaLAPACK
  call setupblacs(mpiglobal, 'grid', bi2d)
#endif
  !   Also make 0d grid containing only the current processes (context 2) 
  !   (context 0 for a rank that is not on bi2d)
  call setupblacs(mpiglobal, '0d', bi0d, np=1)

  ! If fdist, then distribute Hamiltonian on 2d blacs grid, otherwise
  ! use only current rank.
  if(fdist) then 
    bicurrent => bi2d
  else
    bicurrent => bi0d
  end if
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  ! Check if we use MPI, otherwise purely serial computations.
  ! If we distribute the BSE matrix, then loop over Q points serially.
  !---------------------------------------------------------------------------!
#ifdef MPI
  ! Distribute over qmt points and do diagonalization serial 
  if(.not. fdist) then 
    write(unitout, '("Info(",a,"):", a, i0, a)') trim(thisname),&
      & " Distributing qmt-points over ", mpiglobal%procs, " processes."
    call printline(unitout, "+")
    iq1 = firstofset(mpiglobal%rank, nqmtselected, mpiglobal%procs)
    iq2 = lastofset(mpiglobal%rank, nqmtselected, mpiglobal%procs)
  ! Distribute diagonalization and loop over qmt points serially
  else
    write(unitout, '("Info(",a,"):", a)') trim(thisname),&
      & " Distributing BSE matrix, not qmt-points"
    call printline(unitout, "+")
    iq1 = 1
    iq2 = nqmtselected
  end if
#else
  write(unitout, '("Info(",a,"):", a)') trim(thisname),&
    & " Serial execution"
  call printline(unitout, "+")
  iq1 = 1
  iq2 = nqmtselected
#endif
  !---------------------------------------------------------------------------!
  
  !---------------------------------------------------------------------------!
  ! General BSE checks
  !---------------------------------------------------------------------------!
  fcoup = input%xs%bse%coupling
  fchibarq = input%xs%bse%chibarq
  ! If TDA and full Coulomb potential, print warning
  if(.not. fchibarq .and. .not. fcoup) then 
    call printline(unitout, "!")
    write(unitout, '("Warning(",a,"):", a)') trim(thisname),&
      & " TDA using full Chi produces bad results for finite Q, use \bar{Chi}!"
    write(unitout, '("Warning(",a,"):", a)') trim(thisname),&
      & " set input%xs%bse%chibarq = .true."
    call printline(unitout, "!")
  end if
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  ! Assemble and solve BSE for each Q-point in range
  !---------------------------------------------------------------------------
  call setup_bse_type_list(input, bsetypelist, input%xs%BseTypeSet%skipDoneBSEType)
  do bse_type_index = 1, size(bsetypelist)
    input%xs%bse%bsetype = trim(adjustl(bsetypelist(bse_type_index)))

    write(unitout, '("Info(",a,"):", a, a)') trim(thisname),&
      & " BSE type: ", trim(adjustl(input%xs%bse%bsetype))
    call printline(unitout, "+")
    do iqmt = iqmti+iq1-1, iqmti+iq2-1
  
      ! Get full Q vector for info out
      vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)
  
      ! Info out
      call printline(unitout, "-")
      write(unitout, '("Info(",a,"):", a, i3)') trim(thisname),&
        & " Momentum transfer list index: iqmt=", iqmt
      write(unitout, '("Info(",a,"):", a, 3f8.3)') trim(thisname),&
        & " Momentum transfer: vqmtl=", vqmt(1:3)
      call printline(unitout, "-")
  
      
      call bse(iqmt)
  
      ! Info out
      call printline(unitout, "-")
      write(unitout, '("Info(",a,"): BSE finished for iqmt=", i3)')&
        &trim(thisname), iqmt
      call printline(unitout, "-")
  
     end do
     write(unitout, '("Info(",a,"): All done for BSE type ", a)')&
       &trim(thisname), trim(adjustl(input%xs%bse%bsetype))
     if(bse_type_index < size(bsetypelist)) then
       call printline(unitout, " ")
       call printline(unitout, "+")
       call printline(unitout, "+")
       call printline(unitout, " ")
     end if

   end do


  !---------------------------------------------------------------------------!

  if(iq2<0) then
    write(*, '("Info(",a,"): Rank= ", i3, " is idle.")')&
      & trim(thisname), mpiglobal%rank
  end if

  ! Some xas specific finalizations
  if(input%xs%bse%xas .or. input%xs%bse%xes) call xasfinit
  
  call barrier(callername=thisname)

  ! Exit BLACS 
  call exitblacs(bi2d)
  call exitblacs(bi0d)

end subroutine bselauncher
!EOC
