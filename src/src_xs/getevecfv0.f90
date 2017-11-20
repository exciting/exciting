! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getevecfv0
! !INTERFACE:
subroutine getevecfv0(vpl, vgpl, evecfvt)
! !USES:
  use mod_kpoint, only: vkl_ptr
  use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
  use mod_eigensystem, only: nmatmax_ptr
  use mod_misc, only: filext
  use modxs, only: filext0
  use mod_spin, only: nspnfv
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_ematptr
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
  real(8), intent(in) :: vgpl(3,ngkmax0_ptr,nspnfv)
  complex(8), intent(out) :: evecfvt(nmatmax0_ptr,nstfv,nspnfv)

  ! Local variables
  integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
  integer, pointer :: ngk_ptr_save(:,:)
  real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
  character(256) :: filext_save

  nullify(nmatmax_ptr_save)
  nullify(ngkmax_ptr_save)
  nullify(ngk_ptr_save)
  nullify(vkl_ptr_save)
  nullify(vgkl_ptr_save)

  ! Backup pointers to default locations
  nmatmax_ptr_save => nmatmax_ptr
  ngkmax_ptr_save => ngkmax_ptr
  ngk_ptr_save => ngk_ptr
  vkl_ptr_save => vkl_ptr
  vgkl_ptr_save => vgkl_ptr

  filext_save = filext

  ! Set default pointers to bra state quantities
  nmatmax_ptr => nmatmax0_ptr
  ngkmax_ptr => ngkmax0_ptr
  ngk_ptr => ngk0_ptr
  vkl_ptr => vkl0_ptr
  vgkl_ptr => vgkl0_ptr

  filext = filext0

  ! Call to getevecfv with changed default (G+)k-set pointers / matrix size
  call getevecfv(vpl, vgpl, evecfvt)

  ! Restore default pointers
  nmatmax_ptr => nmatmax_ptr_save
  ngkmax_ptr => ngkmax_ptr_save
  ngk_ptr => ngk_ptr_save
  vkl_ptr => vkl_ptr_save
  vgkl_ptr => vgkl_ptr_save

  filext = filext_save

end subroutine getevecfv0
!EOC
