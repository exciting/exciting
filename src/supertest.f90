subroutine supertest
  use m_ematqk_hack 
  use mod_kpointset
  use mod_lattice
  use mod_Gvector
  use mod_Gkvector, only: gkmax 
  use modinput, only: input
  use m_getunit

  implicit none

  integer(4) :: ngridk(3), un, i
  real(8) :: vkloff(3), gmaxvr
  type(k_set) :: kgrid
  type(G_set) :: Ggrid
  type(kq_set) :: kqgrid
  type(Gk_set) :: Gkgrid

  write(*,*) "Hello, this is a test!"
  
  call init0
  call init1 ! gkmax = input:rgkmax / min(rMT)

  call getunit(un)

  ngridk = input%groundstate%ngridk
  vkloff = input%groundstate%vkloff
  gmaxvr = input%groundstate%gmaxvr
  
  call generate_k_vectors(kgrid, bvec, ngridk, vkloff, .false.)
  open(un, file='k_TEST.out', action='write', status='replace')
  call print_k_vectors(kgrid, un)
  close(un)

  intgv(1,:) = [0, 0]
  intgv(2,:) = [0, 0]
  intgv(3,:) = [0, 0]

  call generate_G_vectors(Ggrid, bvec, intgv, gmaxvr)
  open(un, file='G_TEST.out', action='write', status='replace')
  call print_G_vectors(Ggrid, un)
  close(un)

  call generate_Gk_vectors(Gkgrid,kgrid,Ggrid,gkmax)
  open(un, file='Gk_TEST.out', action='write', status='replace')
  do i=1,8
    call print_Gk_vectors(Gkgrid, i, un)
  end do
  close(un)
  
  call generate_kq_vectors(kqgrid,bvec,ngridk, vkloff,.false.)
  open(un, file='kq_TEST.out', action='write', status='replace')
  call print_kq_vectors(kqgrid, un)
  close(un)

  call ematqk_setup(kqgrid, Gkgrid)
  call ematdummy(1)

end subroutine supertest
