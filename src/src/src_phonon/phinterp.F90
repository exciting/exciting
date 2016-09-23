
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phinterp(ngridp,vploff,reducep,tfbz,twrev,fname)
  implicit none
  ! arguments
  integer, intent(in) :: ngridp(3)
  real(8), intent(in) :: vploff(3)
  logical, intent(in) :: reducep, tfbz, twrev
  character(*), intent(in) :: fname
  ! local variables
  real(8) :: boxl(3,4)
  Integer :: nqpti
  Integer, Allocatable :: ivqi (:, :)
  Integer, Allocatable :: iqmapi (:, :, :)
  Real (8), Allocatable :: vqli (:, :)
  Real (8), Allocatable :: vqci (:, :)
  Real (8), Allocatable :: wqpti (:)
  If (allocated(vqli)) deallocate (vqli)
  Allocate (vqli(3, ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (ivqi(3, ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (vqci(3, ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (wqpti(ngridp(1)*ngridp(2)*ngridp(3)))
  Allocate (iqmapi(0:ngridp(1)-1, 0:ngridp(2)-1, 0:ngridp(3)-1))
  boxl (:, 1) = vploff(:) / dble(ngridp(:))
  boxl (:, 2) = boxl (:, 1)
  boxl (:, 3) = boxl (:, 1)
  boxl (:, 4) = boxl (:, 1)
  boxl (1, 2) = boxl (1, 2) + 1.d0
  boxl (2, 3) = boxl (2, 3) + 1.d0
  boxl (3, 4) = boxl (3, 4) + 1.d0
  ! generate q-point set for interpolation
  Call genppts (reducep, tfbz, ngridp, boxl, nqpti, iqmapi, ivqi, vqli, vqci, wqpti)
  deallocate(ivqi,vqci,wqpti,iqmapi)
  ! interpolate
  call writephnlist(nqpti,vqli,twrev,fname)
end subroutine
