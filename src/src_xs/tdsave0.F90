
subroutine tdsave0
  !
  ! to be called after init0, init1 and init2xs to save q=0 variables
  !
  use modmain
  use modtddft
  implicit none

  ! allocate the k-point arrays
  if (allocated(vkl0)) deallocate(vkl0)
  allocate(vkl0(3,nkpt))
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

end subroutine tdsave0
