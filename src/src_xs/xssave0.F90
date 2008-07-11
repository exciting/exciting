
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xssave0
  !
  ! This routine is be called after init0, init1 and init2xs in order to save
  ! variables realted to the k-point set for q=0.
  !
  use modmain
  use modxs
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
  allocate(igkig0(ngkmax,nkpt,nspnfv))
  allocate(vgkl0(3,ngkmax,nkpt,nspnfv))
  allocate(vgkc0(3,ngkmax,nkpt,nspnfv))
  allocate(gkc0(ngkmax,nkpt,nspnfv))
  allocate(tpgkc0(2,ngkmax,nkpt,nspnfv))
  allocate(sfacgk0(ngkmax,natmtot,nkpt,nspnfv))
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
