! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getevecfv0
! !INTERFACE:
subroutine getevecfv0(vpl, vgpl, evecfvt)
! !USES:
  use mod_kpoint, only: vkl, nkpt, nkptnr
  use mod_Gkvector, only: ngkmax, vgkl, ngk
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_spin, only: nspnfv
  use mod_misc, only: filext, task
  use modxs, only: filext0, usefilext0, nmatmax0, ngkmax0, ngk0, vkl0, vgkl0
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   real(8) :: vpl(3)          ! k-point vector in lattice coordinates
!   real(8) :: vgpl(3, ngkmax) ! G+k-vectors in lattice coordinates
! OUT:
!   complex(8) :: evecfvt(nmatmax, nstfv, nspnfv) ! Eigenvectors at that k-point
!
! !DESCRIPTION:
!   This routine is a wrapper for {\tt getevecfv} that changes the $k$ and $G+k$ 
!   quantities in {\tt mod\_kpoint} and {\tt mod\_Gkvector} to the corresponding
!   quantities saved in {\tt modxs} (nmatmax0, vkl0,  ngk0, etc.), changes
!   the file extension in {\tt mod\_misc} accordingly, reads in the Eigenvector
!   and finally restores the original state.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!
!EOP
!BOC
  implicit none

  ! Arguments
  real(8), intent(in) :: vpl(3)
  real(8), intent(in) :: vgpl(3, ngkmax)
  complex(8), intent(out) :: evecfvt(nmatmax, nstfv, nspnfv)

  ! Local variables
  integer :: nmatmaxt, ngkmaxt
  integer, allocatable :: ngkt(:, :)
  real(8), allocatable :: vklt(:, :), vgklt(:, :, :, :)
  character(256) :: filextt

  ! Backup variables for k'-grid quantities
  allocate(ngkt(nspnfv, nkpt))
  allocate(vklt(3, nkptnr))
  allocate(vgklt(3, ngkmax, nspnfv, nkpt))
  nmatmaxt = nmatmax
  ngkmaxt = ngkmax
  ngkt(:, :) = ngk(:, :)
  vklt(:, :) = vkl(:, :)
  vgklt(:, :, :, :) = vgkl(:, :, :, :)
  ! Backup file extension for k'-grid 
  filextt = filext

  ! Copy varialbes of k-grid to default variables
  nmatmax = nmatmax0
  ngkmax = ngkmax0
  ngk(:, :) = ngk0(:, :)
  vkl(:, :) = vkl0(:, :)
  ! Re-allocate array
  deallocate(vgkl)
  allocate(vgkl(3, ngkmax0, nspnfv, nkpt))
  vgkl(:, :, :, :) = vgkl0(:, :, :, :)
  ! Get file extension for k grid
  if(usefilext0) then 
    filext = filext0
  else
    call genfilextread(task)
  end if

  ! Call to getevecfv with changed(G+)k-point sets / matrix size
  call getevecfv(vpl, vgpl, evecfvt)

  ! Restore original k' variables
  nmatmax = nmatmaxt
  ngkmax = ngkmaxt
  ngk(:, :) = ngkt(:, :)
  vkl(:, :) = vklt(:, :)
  ! Re-allocate array
  deallocate(vgkl)
  allocate(vgkl(3, ngkmax, nspnfv, nkpt))
  vgkl(:, :, :, :) = vgklt(:, :, :, :)
  filext = filextt

  ! Free temporary arrays
  deallocate(ngkt, vklt, vgklt)
end subroutine getevecfv0
!EOC
