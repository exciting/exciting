subroutine gv2xmt_spin(is, ia,grhoin, vx, v2xsr)
! i!USES:
use modinput
use mod_Gvector
use mod_muffin_tin
use mod_SHT
use mod_atoms
use mod_potential_and_density
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: ia
real(8), intent(in) :: grhoin(lmmaxvr, nrmtmax)
real(8), allocatable :: grho(:,:)
real(8) :: gvrho(lmmaxvr, nrmtmax, 3), grfmt(lmmaxvr, nrmtmax, 3)
real(8), intent(inout) :: vx(lmmaxvr, nrmtmax)
real(8), intent(in) :: v2xsr(lmmaxvr, nrmtmax)
! local variables
integer::ias,nr, i, j
! allocatable arrays
real(8), allocatable :: v2xmt(:, :,:), gv2x(:, :)
real(8), allocatable :: gv2xsrmt(:, :, :), gv2xsri(:,:,:), v2x(:,:,:) 
allocate(v2xmt(lmmaxvr, nrmtmax,3), gv2x(lmmaxvr, nrmtmax))
allocate(gv2xsrmt(lmmaxvr, nrmtmax, 3),gv2xsri(lmmaxvr,nrmtmax,3),v2x(lmmaxvr,nrmtmax,3))
if (allocated(grho)) deallocate(grho)
allocate(grho(lmmaxvr,nrmtmax))
grho(:,:)=2.d0*grhoin(:,:)
ias=idxas(ia, is)
nr=nrmt(is)
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr,&
& nrmtmax, rhomt(:, :, ias), grfmt)
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, gvrho(:, :, i), lmmaxvr)
  v2x(:,:,i)=v2xsr(:,:)*gvrho(:,:,i)
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, v2x(:,:,i), lmmaxvr, 0.d0, &
       v2xmt(:,:,i), lmmaxvr)
  call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, v2xmt(:,:,i), gv2xsrmt(:,:,:))
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, gv2xsrmt(:, :, i), &
   lmmaxvr, 0.d0, gv2xsri(:,:,i), lmmaxvr)    
end do
gv2x(:, 1:nr)=gv2xsri(:, 1:nr,1) +gv2xsri(:, 1:nr,2) +gv2xsri(:, 1:nr,3) !* &
vx(:,1:nr)=vx(:,1:nr)-gv2x(:,1:nr)

return
end subroutine
