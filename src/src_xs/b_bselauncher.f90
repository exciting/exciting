!BOP
! !ROUTINE: b_bselauncher
! !INTERFACE:
subroutine b_bselauncher
! !USES:
  use modmpi
  use modscl
  use modxs, only: unitout
  use modinput, only: input
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
  character(*), parameter :: thisname = "b_bselauncher"
  integer(4) :: iqmt, iqmti, iqmtf, nqmt, nqmtselected, iq1, iq2
  real(8) :: ts0, ts1
  real(8) :: vqmt(3)
  logical :: fdist, fcoup, fchibarq

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
  ! End timer for init calls
  call timesec(ts1)
  write(unitout, '("Info(",a,"):&
    & Init time:", f12.6)') trim(thisname), ts1 - ts0
  !---------------------------------------------------------------------------!

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
    nqmtselected = iqmtf-iqmti+1
  end if
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
  ! Check if code was compiled with -DSCAL, only then is distribute=true valid
  fdist = input%xs%bse%distribute 
#ifndef SCAL
  if(fdist) then 
    write(*,'("Warning(",a,"):",a)') trim(thisname),&
      & "Setting distribute=false, since code was not&
      & compiled with -DSCAL"
  end if
  input%xs%bse%distribute = .false.
  fdist = .false.
#endif

  !   Make square'ish process grid (context 0)
  call setupblacs(mpiglobal, 'grid', bi2d)
  !   Also make 1d grid with the same number of processes (context 1)
  call setupblacs(mpiglobal, 'row', bi1d, np=bi2d%nprocs)
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
    write(unitout, '("Info(",a,"):", a, i3, a)') trim(thisname),&
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
  !---------------------------------------------------------------------------!
  do iqmt = iqmti+iq1-1, iqmti+iq2-1

    ! Get full Q vector for info out
    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    ! Info out
    call printline(unitout, "-")
    write(unitout, '("Info(",a,"):", a, i3)') trim(thisname),&
      & " Momentum tranfer list index: iqmt=", iqmt
    write(unitout, '("Info(",a,"):", a, 3f8.3)') trim(thisname),&
      & " Momentum tranfer: vqmtl=", vqmt(1:3)
    call printline(unitout, "-")

    ! Assemble and solve BSE
    call b_bse(iqmt)

    ! Info out
    call printline(unitout, "-")
    write(unitout, '("Info(",a,"): Spectrum finished for iqmt=", i3)')&
      &trim(thisname), iqmt
    call printline(unitout, "-")

    if(mpiglobal%rank == 0) then
      write(6, '(a,"BSE(q) Progress:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(iqmt-iqmti+1)/dble(iq2-iq1+1)
      flush(6)
    end if

  end do
  !---------------------------------------------------------------------------!

  if(mpiglobal%rank == 0) then
    write(6, *)
  end if

  if(iq2<0) then
    write(*, '("Info(",a,"): Rank= ", i3, " is idle.")')&
      & trim(thisname), mpiglobal%rank
  end if

  call barrier(callername=thisname)

  ! Exit BLACS 
  call exitblacs(bi2d)
  call exitblacs(bi1d)
  call exitblacs(bi0d)

end subroutine b_bselauncher
!EOC
