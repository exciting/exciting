subroutine gv2xir(grho, vx, v2xsr)
! !USES:
use modinput
use mod_Gvector
use mod_potential_and_density
implicit none
! arguments
real(8), intent(in) :: grho(ngrtot)
real(8), intent(inout) :: vx(ngrtot)
real(8), intent(in) :: v2xsr(ngrtot)
! local variables
integer::ig, ifg, i,j
real(8) :: g2rho(ngrtot)
real(8) :: grho2(ngrtot)
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfft1_v2x(:), zfft2_v2x(:)
real(8), allocatable :: gvrho(:, :)
complex(8), allocatable :: zfft1_rho(:), zfft2_rho(:), rhog(:)
allocate(rfir(ngrtot))
allocate(zfft1_v2x(ngrtot), zfft2_v2x(ngrtot))
allocate(gvrho(ngrtot, 3))
allocate(zfft1_rho(ngrtot), zfft2_rho(ngrtot),rhog(ngrtot))
zfft1_rho(:)=rhoir(:)
call zfftifc(3, ngrid, -1, zfft1_rho)
rhog(:)=zfft1_rho(:)
! |grad rho|
rfir(:)=0.d0
do i=1,3
  zfft2_rho(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2_rho(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1_rho(ifg)), dble(zfft1_rho(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2_rho)
  gvrho(:, i)=dble(zfft2_rho(:))
  zfft1_v2x(:)=v2xsr(:)*gvrho(:,i)
  call zfftifc(3, ngrid, -1, zfft1_v2x)
  ! (grad dxdg2).(grad rho)
  !!add gvrho
  do j=1, 3
     zfft2_v2x(:)=0.d0
     do ig=1, ngvec
        ifg=igfft(ig)
        zfft2_v2x(ifg)=vgc(j, ig)*cmplx(-aimag(zfft1_v2x(ifg)), dble(zfft1_v2x(ifg)), 8)
     end do
     call zfftifc(3, ngrid, 1, zfft2_v2x)
     if (i==j) then
        rfir(:)=rfir(:)+dble(zfft2_v2x(:))
     endif
  end do
enddo
vx(:)=vx(:)-rfir(:) 
return
end subroutine

