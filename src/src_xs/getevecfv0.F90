
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevecfv0(vpl,vgpl,evecfvt)
  use modmain
  use modxs
  ! arguments
  real(8), intent(in) :: vpl(3)
  real(8), intent(in) :: vgpl(3,ngkmax)
  complex(8), intent(out) :: evecfvt(nmatmax,nstfv,nspnfv)
  ! local variables
  integer :: nmatmaxt,ngkmaxt
  integer, allocatable :: ngkt(:,:)
  real(8), allocatable :: vklt(:,:), vgklt(:,:,:,:)
  character(256) :: filextt
  ! copy varialbes of k+(q=0) to default variables
  allocate(ngkt(nspnfv,nkpt))
  allocate(vklt(3,nkptnr))
  allocate(vgklt(3,ngkmax,nspnfv,nkpt))
  nmatmaxt=nmatmax
  nmatmax=nmatmax0
  ngkmaxt=ngkmax
  ngkmax=ngkmax0
  ngkt(:,:)=ngk(:,:)
  ngk(:,:)=ngk0(:,:)
  vklt(:,:)=vkl(:,:)
  vkl(:,:)=vkl0(:,:)
  vgklt(:,:,:,:)=vgkl(:,:,:,:)
  vgkl(:,:,:,:)=vgkl0(:,:,:,:)
  filextt=filext
  ! call to getevecfv with changed (G+)k-point sets / matrix size
  call genfilextread(task)
  call getevecfv(vpl,vgpl,evecfvt)
  ! restore original variables
  nmatmax=nmatmaxt
  ngkmax=ngkmaxt
  ngk(:,:)=ngkt(:,:)
  vkl(:,:)=vklt(:,:)
  vgkl(:,:,:,:)=vgklt(:,:,:,:)
  filext=filextt
  deallocate(ngkt,vklt,vgklt)
end subroutine getevecfv0
