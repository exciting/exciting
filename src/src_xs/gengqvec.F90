


! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gengqvec
! !INTERFACE:


subroutine gengqvec(iq, vpl, vpc, ngp, igpig, vgpl, vgpc, gpc, tpgpc)
! !USES:
use modinput
  use modmain
  use modxs
! !INPUT/OUTPUT PARAMETERS:
!   vpl   : p-point vector in lattice coordinates (in,real(3))
!   vpc   : p-point vector in Cartesian coordinates (in,real(3))
!   ngp   : number of G+p-vectors returned (out,integer)
!   igpig : index from G+p-vectors to G-vectors (out,integer(ngkmax))
!   vgpl  : G+p-vectors in lattice coordinates (out,real(3,ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (out,real(3,ngkmax))
!   gpc   : length of G+p-vectors (out,real(ngkmax))
!   tpgpc : (theta, phi) coordinates of G+p-vectors (out,real(2,ngkmax))
! !DESCRIPTION:
!   Generates a set of ${\bf G+p}$-vectors for the input ${\bf p}$-point with
!   length less than {\tt gkmax}. These are used as the plane waves in the APW
!   functions. Also computes the spherical coordinates of each vector.
!   Based on {\tt gengpvec}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: iq
  real(8), intent(in) :: vpl(3)
  real(8), intent(in) :: vpc(3)
  integer, intent(out) :: ngp
  integer, intent(out) :: igpig(ngqmax)
  real(8), intent(out) :: vgpl(3, ngqmax)
  real(8), intent(out) :: vgpc(3, ngqmax)
  real(8), intent(out) :: gpc(ngqmax)
  real(8), intent(out) :: tpgpc(2, ngqmax)
  ! local variables
  integer::ig, igp
  real(8)::v(3), t1, t2

  integer :: isym, lspl, igpt, ivlt(3)
  real(8) :: vl(3), vc(3), vlt(3), vct(3), vctl(3), s(3, 3), c(3, 3)

  if (input%xs%gqmax.lt.input%structure%epslat) then
     igp=1
     igpig(igp)=igp
     vgpl(:, igp)=vpl(:)
     vgpc(:, igp)=vpc(:)
     call sphcrd(vgpc(1, igp), gpc(igp), tpgpc(1, igp))
     ivgigq(0, 0, 0, iq)=igp
     ngp=1
     return
  end if
  t1=input%xs%gqmax**2
  ivgigq(:, :, :, iq)=0
  igp=0
  do ig=1, ngvec
     v(:)=vgc(:, ig)+vpc(:)
     t2=v(1)**2+v(2)**2+v(3)**2
     if (t2.lt.t1) then
	igp=igp+1
	if (igp.gt.ngqmax) then
	   write(*, *)
	   write(*, '("Error(gengpvec): number of G+p-vectors exceeds ngqmax")')
	   write(*, *)
	   stop
	end if
        ! index to G-vector
	igpig(igp)=ig
        ! G+p-vector in lattice coordinates
	vgpl(:, igp)=dble(ivg(:, ig))+vpl(:)
        ! G+p-vector in Cartesian coordinates
	vgpc(:, igp)=v(:)
        ! G+p-vector length and (theta, phi) coordinates
	call sphcrd(vgpc(1, igp), gpc(igp), tpgpc(1, igp))
        ! map from grid to G+p-vector
	ivgigq(ivg(1, ig), ivg(2, ig), ivg(3, ig), iq)=igp
     end if
  end do
  ngp=igp
  if (input%xs%dbglev.gt.1) then
     write(*, '(a)') 'Debug(gengqvec): igp, isym, lspl, vl, vlt'
     do igp=1, ngp
	vl(:)=dble(ivg(:, igpig(igp)))
	vc=matmul(bvec, vl)
	do isym=1, nsymcrys
	   lspl=lsplsymc(isym)
	   c(:, :)=symlatc(:, :, lspl)
	   s(:, :)=dble(symlat(:, :, lspl))
	   vlt=matmul(vl, s)
	   ivlt=nint(vlt)
	   vct=matmul(vc, c)
	   vctl=matmul(binv, vct)
	   igpt=ivgigq(ivlt(1), ivlt(2), ivlt(3), iq)
!!$           write(*,'(3i6,15f8.4,2x,f8.4)') igp,isym,lspl,vl,vc,vlt,vct,vctl,&
!!$                sum(abs(vlt-vgpl(:,igpt)))
	   write(*, '(3i6, 5x, 3i5, 3x, 3i5)') igp, isym, lspl, nint(vl), nint(vlt)
	end do
     end do
     write(*, *)
  end if
  return
end subroutine gengqvec
!EOC
