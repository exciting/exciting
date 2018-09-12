subroutine gv2xmt(is, ia,grho, vx, v2xsr)
! i!USES:
use modinput
use mod_Gvector
use mod_muffin_tin
use mod_SHT
use mod_atoms
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: ia
real(8), intent(in) :: grho(lmmaxvr, nrmtmax)
real(8) :: gvrho(lmmaxvr, nrmtmax, 3)
real(8), intent(inout) :: vx(lmmaxvr, nrmtmax)
real(8), intent(in) :: v2xsr(lmmaxvr, nrmtmax)
! local variables
integer::ias,nr, i
! allocatable arrays
real(8), allocatable :: v2xsrmt(:, :), gv2x(:, :), gv2xsri(:,:)
real(8), allocatable :: gv2xsrmt(:, :, :)
allocate(v2xsrmt(lmmaxvr, nrmtmax), gv2x(lmmaxvr, nrmtmax),gv2xsri(lmmaxvr,nrmtmax))
allocate(gv2xsrmt(lmmaxvr, nrmtmax, 3))
ias=idxas(ia, is)
nr=nrmt(is)
! convert v2xsr to spherical harmonics
!v2xsr(:,:)=v2xsr(:,:)*sqrt(grho(:,:))
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, v2xsr, lmmaxvr, 0.d0, &
 v2xsrmt, lmmaxvr)
! compute grad v2xsr
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, v2xsrmt, gv2xsrmt)
gv2x(:, 1:nr)=0.d0
do i=1, 3
!  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grho(:, :), &
!   lmmaxvr, 0.d0, gvrho(:, :, i), lmmaxvr)
  gvrho(:,:,i)=0.d0
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, gv2xsrmt(:, :, i), &
   lmmaxvr, 0.d0, gv2xsri, lmmaxvr)
!this should be changed, verifyyyy!!!
  gv2x(:, 1:nr)=gv2x(:, 1:nr)+gv2xsri(:, 1:nr)* &
                (gvrho(:, 1:nr, i)/grho(:,1:nr))
end do
!gv2xsr(:, 1:nr)=sqrt(gv2xsri(:, 1:nr, 1)**2+gv2xsri(:, 1:nr, 2)**2+gv2xsri(:, 1:nr, 3)**2)
vx(:,1:nr)=vx(:,1:nr)-gv2x(:,1:nr)
return
end subroutine
