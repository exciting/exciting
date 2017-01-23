subroutine setgqvars(qset, gset, gqset)
  use mod_kpointset
  use mod_kpoint, only: nkpt
  use mod_atoms, only: natmtot
  use mod_qpoint, only: nqpt, vql
  use modxs, only: ngqmax, ngq, vgql, vgqc,&
                 & gqc, tpgqc, sfacgq, igqig,&
                 & ivgigq, intgqv, ikmapikq

  implicit none
  
  type(q_set), intent(in) :: qset
  type(G_set), intent(in) :: gset
  type(Gk_set), intent(in) :: gqset

  integer(4) :: ispin, iq, igq, ig, ivg(3)

  ispin = 1

  ! Set new number of reduced q-points 
  nqpt = qset%qset%nkpt

  ! Set new reduced q-vectors 
  if(allocated(vql)) deallocate(vql)
  allocate(vql(3,nqpt))
  vql(1:3,1:nqpt) = qset%qset%vkl(1:3,1:nqpt)

  ! Set new number of G+q vectors
  ngqmax = gqset%ngkmax
  if(allocated(ngq)) deallocate(ngq)
  allocate(ngq(nqpt))
  ngq(1:nqpt) = gqset%ngk(ispin,1:nqpt)

  ! Set new G+q vectors
  ! G+q in lattice coordinates
  if(allocated(vgql)) deallocate(vgql)
  allocate(vgql(3, ngqmax, nqpt))
  vgql = 0.0d0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      vgql(1:3, igq, iq) = gqset%vgkl(1:3, igq, ispin, iq)
    end do
  end do
  ! G+q in Cartesian coordinates
  if(allocated(vgqc)) deallocate(vgqc)
  allocate(vgqc(3, ngqmax, nqpt))
  vgqc = 0.0d0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      vgqc(1:3, igq, iq) = gqset%vgkc(1:3, igq, ispin, iq)
    end do
  end do
  ! |G+q|
  if(allocated(gqc)) deallocate(gqc)
  allocate(gqc(1:3, gqset%ngkmax))
  gqc = 0.0d0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      gqc(igq, iq) = gqset%gkc(igq, ispin, iq)
    end do
  end do
  ! Theta(G+q), Phi(G+q)
  if(allocated(tpgqc)) deallocate(tpgqc)
  allocate(tpgqc(2, ngqmax, nqpt))
  tpgqc = 0.0d0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      tpgqc(1:2, igq, iq) = gqset%tpgkc(1:2, igq, ispin, iq)
    end do
  end do
  ! Structure factor
  if(allocated(sfacgq)) deallocate(sfacgq)
  allocate(sfacgq(ngqmax, natmtot, nqpt))
  sfacgq = 0.0d0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      sfacgq(igq, :, iq) = gqset%sfacgk(igq, :, ispin, iq)
    end do
  end do
  ! Set new map from G+q index to G index
  if(allocated(igqig)) deallocate(igqig)
  allocate(igqig(ngqmax, nqpt))
  igqig = 0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      igqig(igq, iq) = gqset%igkig(igq,ispin,iq)
    end do
  end do
  write(*,*) intgqv
  ! Set 3d G index + 1d iq --> 1d G+q index map
  if(allocated(ivgigq)) deallocate(ivgigq)
  allocate (ivgigq(intgqv(1,1):intgqv(1,2), &
  &                intgqv(2,1):intgqv(2,2), &
  &                intgqv(3,1):intgqv(3,2), nqpt))
  ivgigq = 0
  do iq = 1, nqpt
    do igq = 1, ngq(iq)
      ig = igqig(igq, iq)
      ivg = gset%ivg(1:3, ig)
      ivgigq(ivg(1), ivg(2), ivg(3), iq) = igq
    end do
  end do

  ! Set mapping between k+q and k'
  if(allocated(ikmapikq)) deallocate(ikmapikq)
  allocate(ikmapikq(nkpt, nqpt))
  do iq = 1, nqpt
    ikmapikq(1:nkpt, iq) = qset%ikiq2ikp_nr(1:nkpt, qset%qset%ikp2ik(iq))
  end do
end subroutine setgqvars
