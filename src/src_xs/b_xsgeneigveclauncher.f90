! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine b_xsgeneigveclauncher
  use modmpi
  use modinput, only: input
  use mod_qpoint, only: nqpt, vql
  use modxs, only: nqmt, vqlmt, ivgmt, tscreen, qvkloff, unitout
  use m_genfilname, only: genfilname
  use m_writegqpts, only: writegqpts
  use mod_misc, only: filext
  use mod_xsgrids
  use mod_Gkvector, only: gkmax

  implicit none

  ! Local variables
  character(*), parameter :: thisnam = 'b_xsgeneigveclauncer'
  integer(4) :: iq, qi, qf
  logical :: tminus, tmqmt
  logical :: firstisgamma
  real(8), allocatable :: vkloff_kqmtm(:,:), vkloff_mkqmtp(:,:)
  real(8), parameter :: epslat=1.d-6

  write(*,*) "b_xsgeneigveclauncher here at rank", rank
  write(*,*) "use screening parameters = ", tscreen 

  ! Initialize universal variables
  call init0
  ! k-point setup
  ! Also allocated the radial functions (mod_APW_LO)
  call init1
  ! q-point and qmt-point setup
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

  if(NORM2(vqlmt(1:3,1)) > epslat) then 
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

  if(rank .eq. 0) then 
    ! Write Q-points
    call writeqpts
    do iq = 1, nqpt
      if(tscreen) then 
        call genfilname(iqmt=iq, scrtype='', setfilext=.true.)
      else
        call genfilname(iqmt=iq, setfilext=.true.)
      end if
      call writegqpts(iq, filext)
    end do
  end if

  ! Offsets for k+qmt grid are stored in qvkloff

  ! Groundstate offsets for the k-qmt grid
  allocate(vkloff_kqmtm(3,nqpt))

  ! Groundstate offsets for the -(k+qmt) grid
  allocate(vkloff_mkqmtp(3,nqpt))

  ! For each Q-point in the Q-point list generate grids and
  ! save offsets.
  do iq = 1, nqpt

    call xsgrids_init(vqlmt(1:3, iq), gkmax)
    if(mpiglobal%rank == 0) then 
      call xsgrids_write_grids(iq)
    end if

    ! Offset for (k-qmt) grid
    vkloff_kqmtm(1:3,iq) = k_kqmtm%kqmtset%vkloff

    ! Offset for -(k+qmt) grid
    vkloff_mkqmtp(1:3,iq) = mkqmtp%kset%vkloff

    !write(*,*) "iq=", iq
    !write(*,*) "vqlmt(1:3,iq)=", vqlmt(1:3,iq)
    !write(*,*) "off: k+qmt", qvkloff(1:3,iq)
    !write(*,*) "off: k-qmt", vkloff_kqmtm(1:3,iq)
    !write(*,*) "off: -(k+qmt)", vkloff_mkqmtp(1:3,iq)

    call xsgrids_finalize()

  end do

  ! Depending on the BSE hamiltonian to be constructed 
  ! the eigensolutions on differing grids are needed.
  if(input%xs%bse%coupling) then 
    ! Full bse with time inverted anti-resonant basis
    if(input%xs%bse%ti) then 
      tmqmt = .false.
      tminus = .false.
      ! (k+qmt) grid
      call b_xsgeneigvec(1, nqpt, vql(1:3,1:nqpt), qvkloff(1:3,1:nqpt), tscreen, tmqmt, tminus)
      tminus = .true.
      ! -(k+qmt) grid
      call b_xsgeneigvec(1, nqpt, vql(1:3,1:nqpt), vkloff_mkqmtp(1:3,1:nqpt), tscreen, tmqmt, tminus)
    ! Full bse 
    else
      tmqmt = .false.
      tminus = .false.
      ! (k+qmt) grid
      call b_xsgeneigvec(1, nqpt, vql(1:3,1:nqpt), qvkloff(1:3,1:nqpt), tscreen, tmqmt, tminus)
      tmqmt = .true.
      if(nqpt > 1) then 
        ! (k-qmt) grid
        call b_xsgeneigvec(2, nqpt, vql(1:3,1:nqpt), vkloff_kqmtm(1:3,1:nqpt), tscreen, tmqmt, tminus)
      end if
    end if
  ! TDA case
  else
    tmqmt = .false.
    tminus = .false.
    ! (k+qmt) grid
    call b_xsgeneigvec(1, nqpt, vql(1:3,1:nqpt), qvkloff(1:3,1:nqpt), tscreen, tmqmt, tminus)
  end if

  deallocate(vkloff_kqmtm)
  deallocate(vkloff_mkqmtp)

end subroutine b_xsgeneigveclauncher
