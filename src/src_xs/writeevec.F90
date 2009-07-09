


! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine writeevec(vq, voff, filxt)
  use modmain
use modinput
  use modmpi
  use modxs
  use m_gndstateq
  use m_filedel
  implicit none
  ! arguments
  real(8), intent(in) :: vq(3), voff(3)
  character(*), intent(in) :: filxt
  ! local variables
  integer :: ik, j
  complex(8), allocatable :: apwalm(:, :, :, :)
  ! read from STATE.OUT exclusively
  isreadstate0=.true.
  ! SCF calculation with one cycle
  call gndstateq(voff, filxt)
  if (allocated(evecfv)) deallocate(evecfv)
  allocate(evecfv(nmatmax, nstfv, nspnfv))
  allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
  allocate(apwcmt(nstfv, apwordmax, lmmaxapw, natmtot))
  allocate(locmt(nstfv, nlomax, -lolmax:lolmax, natmtot))
  ! delete existing coefficients files
  if (rank.eq.0) call filedel('APWCMT'//trim(filxt))
  if (rank.eq.0) call filedel('LOCMT'//trim(filxt))
  call genparidxran('k', nkpt)
  do ik=kpari, kparf
     apwcmt(:, :, :, :)=zzero
     locmt(:, :, :, :)=zzero
     call getevecfv(vkl(1, ik), vgkl(1, 1, 1, ik), evecfv)
     call match(ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), sfacgk(1, 1, 1, ik), &
	  apwalm)
     call genapwcmt(input%groundstate%lmaxapw, ngk(1, ik), 1, nstfv, apwalm, evecfv, apwcmt)
     call genlocmt(ngk(1, ik), 1, nstfv, evecfv, locmt)
     do j=0, procs-1
	if (rank.eq.j) then
	   call putapwcmt('APWCMT'//trim(filxt), ik, vkl(1, ik), vq, apwcmt)
	   call putlocmt('LOCMT'//trim(filxt), ik, vkl(1, ik), vq, locmt)
	end if
	call barrier
     end do
  end do
  call endloopbarrier(nkpt, procs)
  isreadstate0=.false.
  deallocate(evecfv, apwalm, apwcmt, locmt)
end subroutine writeevec
