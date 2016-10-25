! Copyright (C) 2005-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmatxs
! !INTERFACE:
subroutine writepmatxs
! !USES:
  use modinput, only: input
  use modmpi, only: procs, rank, firstofset, lastofset, barrier
  use mod_misc, only: task, filext
  use mod_kpoint, only: nkpt, vkl
  use mod_Gkvector, only: ngkmax, vgkl, ngk, gkc, tpgkc, sfacgk, igkig, vgkc
  use mod_APW_LO, only: apwordmax, nlotot, nlomax, lolmax
  use mod_muffin_tin, only: lmmaxapw
  use mod_atoms, only: natmtot
  use mod_eigensystem, only: nmatmax 
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use modxas, only: ncg
  use modxs, only: tscreen, fnpmat, fnpmat_t, kpari,&
                  & kparf, hybridhf, ripaa, ripalo,&
                  & riploa, riplolo, apwcmt, locmt,&
                  & unitout
  use m_putpmat
  use m_genfilname
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}, {\tt PMAT\_XS.OUT} or
!   {\tt PMAT\_SCR.OUT} depending on the context of the execution.
!
! !REVISION HISTORY:
!   Created 2006 (S. Sagmeister)
!   Modifications, August 2010 (S. Sagmeister)
!EOP
!BOC

  implicit none

  ! Local variables
  integer :: ik, reclen
  character(32) :: fnam
  complex(8), allocatable :: apwalmt(:, :, :, :)
  complex(8), allocatable :: evecfvt(:, :)
  complex(8), allocatable :: evecsvt(:, :)
  complex(8), allocatable :: pmat(:, :, :)
  character(256) :: string
  logical :: fast

  ! External functions
  logical, external :: tqgamma

  ! Check if fast (default) version of matrix elements is used
  fast=.false.
  if(associated(input%properties)) then
    if(associated(input%properties%momentummatrix)) then
      if(input%properties%momentummatrix%fastpmat) fast=.true.
    end if
  end if

  ! Task 120 is 'writepmat'
  fast=(task.ne.120).or.((task.eq.120).and.fast)

  ! Check if Q-point is gamma point
  if(task .ne. 120) then
    if(.not.tqgamma(1)) return
  end if

  ! Check if this routine is called for the screening
  tscreen=(task .ge. 400) .and. (task .le. 499)

  if((task .eq. 120).or. tscreen) then
    fnam = 'PMAT'
    call genfilname(basename=trim(fnam), appfilext=.true.,&
      & filnam=fnpmat)
    call genfilname(basename=trim(fnam), procs=procs, rank=rank,&
      & appfilext=.true., filnam=fnpmat_t)
  else
    fnam = 'PMAT_XS'
    call genfilname(basename=trim(fnam), filnam=fnpmat)
    call genfilname(basename=trim(fnam), procs=procs, rank=rank,&
      & filnam=fnpmat_t)
  end if

  ! Initialise universal variables
  call init0
  call init1
  if(task .ne. 120) call init2

  ! Generate index ranges for parallel execution
  call genparidxran('k', nkpt)
  ! K-point interval for process
  !kpari = firstofset(rank, nkpt)
  !kparf = lastofset(rank, nkpt)

  allocate(apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
  allocate(evecfvt(nmatmax, nstfv))
  allocate(evecsvt(nstsv, nstsv))

  ! Allocate the momentum matrix elements array
  if((input%xs%bse%xas) .and. (task .le. 400)) then ! Allocation for xas calculation
    allocate(pmat(3, ncg, nstsv))
  else
    allocate(pmat(3, nstsv, nstsv))
  end if  

  if(task .eq. 120) then
    ! Read density and potentials from file
    if(hybridhf) then
      ! In case of hf hybrids use pbe potential
      string=filext
      filext='_PBE.OUT'
      call readstate
      filext=string
    else
      call readstate
    end if
    ! Find the new linearisation energies
    call linengy
    ! Generate the apw radial functions
    call genapwfr
    ! Generate the local-orbital radial functions
    call genlofr
    ! Update potential in case if hf hybrids
    if(hybridhf) call readstate
    ! Find the record length
    inquire(iolength=reclen) pmat
    open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT',&
      & status='REPLACE',recl=reclen)
  end if

  ! Get eigenvectors for q=0
  if((.not. tscreen) .and. (task .ne. 120)) call genfilname(iqmt=0, setfilext=.true.)

  ! Generate band combinations
  if(task .eq. 120) then
    call ematbdcmbs(0)
  else
    call ematbdcmbs(1)
  end if

  if(fast) then
    if(allocated(apwcmt)) deallocate(apwcmt)
    allocate(apwcmt(nstfv, apwordmax, lmmaxapw, natmtot))
    if(allocated(ripaa)) deallocate(ripaa)
    allocate(ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw,natmtot, 3))
    if(nlotot .gt. 0) then
      if(allocated(locmt)) deallocate(locmt)
      allocate(locmt(nstfv, nlomax,-lolmax:lolmax, natmtot))
      if(allocated(ripalo)) deallocate(ripalo)
      allocate(ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, natmtot, 3))
      if(allocated(riploa)) deallocate(riploa)
      allocate(riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, natmtot, 3))
      if(allocated(riplolo)) deallocate(riplolo)
      allocate(riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, natmtot, 3))
    end if

    ! Calculate gradient of radial functions times spherical harmonics
    call pmatrad
  end if

  kloop: do ik = kpari, kparf

    if(task .ne. 120) call chkpt(2, (/ task, ik /), 'ematqk:&
     & task, k - point index; momentum matrix elements')

    ! Get the eigenvectors and values from file
    call getevecfv(vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
    call getevecsv(vkl(1, ik), evecsvt)

    ! Find the matching coefficients
    call match(ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik),&
      & sfacgk(1, 1, 1, ik), apwalmt)

    if(fast) then

      ! Generate apw expansion coefficients for muffin-tin
       call genapwcmt(input%groundstate%lmaxapw, ngk(1, ik), 1,&
         & nstfv, apwalmt, evecfvt, apwcmt)

      ! Generate local orbital expansion coefficients for muffin-tin
      if(nlotot .gt. 0) & 
        & call genlocmt(ngk(1, ik), 1, nstfv, evecfvt, locmt)

      ! Calculate the momentum matrix elements
      if((input%xs%bse%xas) .and. (task .le. 400)) then
        call genpmatcorxs(ik, ngk(1, ik), apwalmt, evecfvt, evecsvt, pmat)
      else
        call genpmatxs(ngk(1, ik), igkig(1, 1, ik),&
         & vgkc(1, 1, 1, ik), evecfvt, evecsvt, pmat)
      end if

    else

      ! Calculate the momentum matrix elements
      call genpmat(ngk(1, ik), igkig(1, 1, ik), vgkc(1, 1, 1, ik),&
       & apwalmt, evecfvt, evecsvt, pmat)

    end if

    if(task .eq. 120) then
      ! Write the matrix elements to direct-access file
      write(50,rec=ik) pmat
    else
      ! Parallel write
      call putpmat(ik, .true., trim(fnpmat), pmat)
    end if

  end do kloop

  call barrier

  inquire(iolength=reclen) vkl(:, ik), nstsv, pmat
  deallocate(apwalmt, evecfvt, evecsvt, pmat)

  if(fast) then
    deallocate(apwcmt)
    deallocate(ripaa)
    if(nlotot .gt. 0) then
      deallocate(locmt)
      deallocate(ripalo, riploa, riplolo)
    end if
  end if

  call barrier

  if(.not. input%sharedfs) call cpfiletonodes(trim(fnpmat))

  if(task .eq. 120) then
    close(50)
    if(rank==0) then
      write(*,*)
      write(*, '("Info(writepmatxs):")')
      write(*, '(" momentum matrix elements written to file PMAT.OUT")')
      write(*,*)
    end if
  else
    write(unitout, '(a)') "Info(writepmatxs): momentum matrix elements finished"
  end if

  ! Reset global file extension to default
  call genfilname(setfilext=.true.)

end subroutine writepmatxs
!EOC
