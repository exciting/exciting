! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2offs(voffk, reduceq)
  use modinput
  use modmpi
  use mod_muffin_tin
  use mod_atoms
  use mod_kpoint
  use mod_qpoint
  use mod_Gkvector
  use modxs

  implicit none

  real(8), intent(in) :: voffk(3) ! q-grid offset in k-grid coordinates
  logical, intent(in) :: reduceq  ! whether or not to reduce the q-grid

  integer(4) :: iq, iv(3)
  real(8) :: boxl(3, 4)



  ! Map offset to first k-parallelepiped
  call r3frac(input%structure%epslat, voffk, iv)

  ! q-point grid with offset
  boxl(:,1) = voffk/dble(ngridq)
  boxl(:, 2) = boxl(:, 1)
  boxl(:, 3) = boxl(:, 1)
  boxl(:, 4) = boxl(:, 1)
  boxl(1, 2) = boxl(1, 2) + 1.d0
  boxl(2, 3) = boxl(2, 3) + 1.d0
  boxl(3, 4) = boxl(3, 4) + 1.d0

  ! Set mod_qpoint variables to a reduced shifted q grid

  if(allocated(ivq)) deallocate(ivq)
  allocate(ivq(3, ngridq(1)*ngridq(2)*ngridq(3)))
  if(allocated(vql)) deallocate(vql)
  allocate(vql(3, ngridq(1)*ngridq(2)*ngridq(3)))
  if(allocated(vqc)) deallocate(vqc)
  allocate(vqc(3, ngridq(1)*ngridq(2)*ngridq(3)))
  if(allocated(wqpt)) deallocate(wqpt)
  allocate(wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
  if(allocated(iqmap)) deallocate(iqmap)
  allocate(iqmap(0:ngridq(1)-1, 0:ngridq(2)-1, 0:ngridq(3)-1))
  ! Generate reduced/non-reduced q-point set
  call genppts(reduceq, input%xs%bse%fbzq, ngridq,&
    & boxl, nqpt, iqmap, ivq, vql, vqc, wqpt)

  ! Derive offset of k+q grid from k-grid offset vkloff and q points
  ! Also generate k,q --> k' mapping
  if(allocated(qvkloff)) deallocate(qvkloff)
  allocate(qvkloff(3, 0:nqpt))
  if(allocated(ikmapikq)) deallocate(ikmapikq)
  allocate(ikmapikq(nkpt, nqpt))
  qvkloff(:, 0) = input%groundstate%vkloff(:)
  do iq = 1, nqpt
    ! offset for k+q-point set derived from q-point
    call genqvkloff(vql(1,iq), qvkloff(1,iq))
    ! map from k-point index to k+q point index for same k
    call findkmapkq(vql(1,iq), qvkloff(1,iq), ikmapikq(1,iq))
  end do

  !---------------------!
  !     G+q-point set   !
  !---------------------!
  ! warning for small gqmax
  if(input%xs%gqmax .ge. gkmax) then
    write(*,'(a, 2g18.10)') 'Warning(init2offs/xs): input%xs%gqmax >= gkmax: ', &
    &  input%xs%gqmax, gkmax
  end if

  ! check consistency with fft G-vector array
  if(input%groundstate%gmaxvr .lt. 2*gkmax + input%xs%gqmax) then
    write(*,*)
    write(*,'("Error(init2offs): gmaxvr < 2 gkmax + gqmax : ",2g18.10)') &
      & input%groundstate%gmaxvr, 2*gkmax+input%xs%gqmax
    write(*,'(" gmaxvr : ",g18.10)') input%groundstate%gmaxvr
    write(*,'(" gkmax  : ",g18.10)') gkmax
    write(*,'(" gqmax  : ",g18.10)') input%xs%gqmax
    write(*,'(" maximum value for gqmax    : ",g18.10)') &
      & input%groundstate%gmaxvr - 2*gkmax
    write(*,'(" increase gmaxvr and re-do scf calculation or decrease gqmax or rgkmax.")')
    write(*,*)
    call terminate
  end if

  !! Write modxs G+q quantities

  ! maximum number of g+q vectors for all q
  call getngqmax
  ! allocate the g+q-vector arrays
  if(allocated(ngq)) deallocate(ngq)
  allocate(ngq(nqpt))
  if(allocated(igqig)) deallocate(igqig)
  allocate(igqig(ngqmax, nqpt))
  if(allocated(vgql)) deallocate(vgql)
  allocate(vgql(3, ngqmax, nqpt))
  if(allocated(vgqc)) deallocate(vgqc)
  allocate(vgqc(3, ngqmax, nqpt))
  if(allocated(gqc)) deallocate(gqc)
  allocate(gqc(ngqmax, nqpt))
  if(allocated(tpgqc)) deallocate(tpgqc)
  allocate(tpgqc(2, ngqmax, nqpt))
  if(allocated(sfacgq)) deallocate(sfacgq)
  allocate(sfacgq(ngqmax, natmtot, nqpt))
  if(allocated(ylmgq)) deallocate(ylmgq)
  allocate(ylmgq(lmmaxapw, ngqmax, nqpt))
  if(allocated(ivgigq)) deallocate(ivgigq)
  allocate(ivgigq(intgqv(1,1):intgqv(1,2), &
      &                intgqv(2,1):intgqv(2,2), &
      &                intgqv(3,1):intgqv(3,2), nqpt))

  do iq = 1, nqpt
    ! generate g+q vectors
    call gengqvec(iq, vql(1,iq), vqc(1,iq), ngq(iq), &
    &             igqig(1,iq), vgql(1,1,iq), vgqc(1,1,iq), &
    &             gqc(1,iq), tpgqc(1,1,iq))
    ! generate structure factors for g-vectors
    call gensfacgp(ngq(iq), vgqc(1,1,iq), ngqmax, sfacgq(1,1,iq))
    ! spherical harmonics for g+q-vectors
    call genylmgq(iq, input%groundstate%lmaxvr)
  end do

  !---------------------------!
  !     coulomb potential     !
  !---------------------------!
  if(allocated(sptclg)) deallocate(sptclg)
  allocate(sptclg(ngqmax, nqpt))
  do iq = 1, nqpt
    ! set up coulomb potential square root
    call genptclg('nocutoff', ngqmax, ngq(iq), vgqc(:, :, iq), &
      & gqc(:, iq), sptclg(:, iq))
  end do

end subroutine
