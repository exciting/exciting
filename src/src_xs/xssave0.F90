
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: xssave0
! !INTERFACE:
subroutine xssave0
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   This routine is be called after init0, init1 and init2xs in order to save
!   variables realted to the k-point set for q=0.
!
! !REVISION HISTORY:
!   Created March 2005 (Sagmeister)
!EOP
!BOC
  implicit none
  ! allocate the k-point arrays
  if (allocated(vkl0)) deallocate(vkl0)
  allocate(vkl0(3,nkptnr))
  ! allocate the G+k-point arrays
  if (allocated(ngk0)) deallocate(ngk0)
  if (allocated(igkig0)) deallocate(igkig0)
  if (allocated(vgkl0)) deallocate(vgkl0)
  if (allocated(vgkc0)) deallocate(vgkc0)
  if (allocated(gkc0)) deallocate(gkc0)
  if (allocated(tpgkc0)) deallocate(tpgkc0)
  if (allocated(sfacgk0)) deallocate(sfacgk0)
  allocate(ngk0(nkpt,nspnfv))
  allocate(igkig0(ngkmax,nspnfv,nkpt))
  allocate(vgkl0(3,ngkmax,nspnfv,nkpt))
  allocate(vgkc0(3,ngkmax,nspnfv,nkpt))
  allocate(gkc0(ngkmax,nspnfv,nkpt))
  allocate(tpgkc0(2,ngkmax,nspnfv,nkpt))
  allocate(sfacgk0(ngkmax,natmtot,nspnfv,nkpt))
  ! overlap and Hamiltonian matrix sizes
  if (allocated(nmat0)) deallocate(nmat0)
  allocate(nmat0(nkpt,nspnfv))
  ! save variables for k-vectors
  nkpt0=nkpt
  vkl0(:,:) = vkl(:,:)
  ! save variables for G+k-vectors
  ngkmax0 = ngkmax
  ngk0(:,:) = ngk(:,:)
  igkig0(:,:,:) = igkig(:,:,:)
  vgkl0(:,:,:,:) = vgkl(:,:,:,:)
  vgkc0(:,:,:,:) = vgkc(:,:,:,:)
  gkc0(:,:,:) = gkc(:,:,:)
  tpgkc0(:,:,:,:) = tpgkc(:,:,:,:)
  sfacgk0(:,:,:,:) = sfacgk(:,:,:,:)
  ! save variables for overlap and Hamiltonian matrix sizes
  nmatmax0 = nmatmax
  nmat0(:,:) = nmat(:,:)
end subroutine xssave0
!EOC
