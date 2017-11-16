!BOP
! !ROUTINE: b_xsgeneigveclauncher
! !INTERFACE:
subroutine b_xsgeneigveclauncher()
! !USES:
  use modmpi
  use modinput, only: input
  use mod_qpoint, only: nqpt, vql
  use modxs, only: unitout, totalqlmt, tscreen, qvkloff
  use m_genfilname, only: genfilname
  use m_writegqpts, only: writegqpts
  use mod_misc, only: filext
  use mod_xsgrids
  use mod_Gkvector, only: gkmax
! !DESCRIPTION:
!   Wrapper routine for \texttt{b\_xsgeneigvec}. Launches one-shot ground state 
!   calculations needed for Q-dependent BSE. Note: First Q-point in the Q-point
!   list needs to be the Gamma point.
!
! !REVISION HISTORY:
!   Created. 2017 (Aurich)
!EOP
!BOC

  implicit none

  ! Local variables
  character(*), parameter :: thisnam = 'b_xsgeneigveclauncer'
  integer(4) :: iq, qf
  logical :: tmqmt
  logical :: firstisgamma
  real(8), allocatable :: vkloff_kqmtm(:,:)
  real(8), allocatable :: vkloff_kqmtp(:,:)
  real(8), parameter :: epslat=1.d-6
  logical :: fwg

  ! Initialize universal variables
  call init0
  ! k-point setup
  ! Also allocated the radial functions (mod_APW_LO)
  call init1
  ! q-point and qmt-list setup
  !   Init 2 sets up (task 301/401):
  !   * A list of momentum transfer vectors form the q-point list 
  !     (modxs::vqmtl and mod_qpoint::vql)
  !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
  !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
  !   * G+qmt quantities (modxs)
  !   * The square root of the Coulomb potential for the G+qmt points
  !   * Reads STATE.OUT
  !   * Generates radial functions (mod_APW_LO)
  call init2

  ! Check if first entry in Q-point list is the Gamma point.
  if(NORM2(totalqlmt(1:3,1)) > epslat) then 
    firstisgamma = .false.
  else
    firstisgamma = .true.
  end if
  if(.not. firstisgamma) then 
    if(rank == 0) then 
      write(*,*) "First Q-point needs to be the gamma point."
    end if
    call terminate
  end if

  ! For the screening only the reference grid is needed.
  ! Note: tscreen is true, for the tasks scrgeneigvec.
  if(tscreen) then
    qf = 1
  else
    qf = nqpt
  end if

  ! Write Q-points
  if(rank .eq. 0) then 
    call writeqmtpts
    call writeqpts
    do iq = 1, qf
      if(tscreen) then 
        call genfilname(scrtype='', setfilext=.true.)
      else
        call genfilname(iqmt=iq, setfilext=.true.)
      end if
      call writegqpts(iq, filext)
    end do
  end if

  ! BSE uses +/- Q/2 grids (symmetrically shifted k grids
  ! for each momentum transfer Q), while TDDFT uses +Q grids (asymmetrically shifted).
  if(input%xs%xstype == "BSE") then 

    call printline(unitout, "=")
    write(unitout, '("One-shot GS runs for BSE calculations")')
    call printline(unitout, "=")

    ! Offsets for k+qmt/2 grids
    allocate(vkloff_kqmtp(3,nqpt))
    ! Offsets for the k-qmt/2 grid
    allocate(vkloff_kqmtm(3,nqpt))

    ! For each Q-point in the Q-point list generate grids and
    ! save offsets which are needed to call the gound state routine.

    ! Write out detailed information about the used grids.
    if(input%xs%writexsgrids) then 
      fwg = .true.
    else
      fwg = .false.
    end if

    do iq = 1, nqpt

      call xsgrids_init(totalqlmt(1:3, iq), gkmax, makegk_=fwg, makegq_=fwg)
      if(mpiglobal%rank == 0 .and. fwg) then 
        call xsgrids_write_grids(iq)
      end if

      !! Only save offsets
      ! Offset for (k+qmt/2) grid
      vkloff_kqmtp(1:3,iq) = k_kqmtp%kqmtset%vkloff
      ! Offset for (k-qmt/2) grid
      vkloff_kqmtm(1:3,iq) = k_kqmtm%kqmtset%vkloff

      ! Clear grids again
      call xsgrids_finalize()

    end do

    ! Depending on the BSE Hamiltonian to be constructed 
    ! the eigensolutions on differing grids are needed.

    !! (k+qmt/2) grid
    tmqmt=.false.
    call printline(unitout, "+")
    write(unitout, '("One-shot GS runs for k+qmt/2 grids")')
    call printline(unitout, "+")
    ! Do one-shot GS calculations for qmt-points number 1 to qf
    call b_xsgeneigvec(1, qf, nqpt, vql(1:3,1:nqpt), vkloff_kqmtp(1:3,1:nqpt),&
      & tscreen, tmqmt)

    !! (k-qmt/2) grid
    ! Skip Gamma (assumed to be the first entry)
    if(qf > 1) then 
      tmqmt=.true.
      call printline(unitout, "+")
      write(unitout, '("One-shot GS runs for k-qmt/2 grids")')
      call printline(unitout, "+")
      ! Do one-shot GS calculations for qmt-points number 2 to qf
      call b_xsgeneigvec(2, qf, nqpt, vql(1:3,1:nqpt), vkloff_kqmtm(1:3,1:nqpt),&
        & tscreen, tmqmt)
    end if

    deallocate(vkloff_kqmtp)
    deallocate(vkloff_kqmtm)

  else if(input%xs%xstype == "TDDFT") then

    call printline(unitout, "=")
    write(unitout, '("One-shot GS runs for TDDFT calculations")')
    call printline(unitout, "=")

    !! (k+qmt) grid
    tmqmt=.false.
    call printline(unitout, "+")
    write(unitout, '("One-shot GS runs for k+qmt grids")')
    call printline(unitout, "+")
    ! Do one-shot GS calculations for qmt-points number 1 to qf
    call b_xsgeneigvec(1, qf, nqpt, vql(1:3,1:nqpt), qvkloff(1:3,1:nqpt),&
      & tscreen, tmqmt)

  else

    write(*,*) "What? We should not be here."

  end if

end subroutine b_xsgeneigveclauncher
!EOC
