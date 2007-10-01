subroutine gradrf(rfmt,rfir,grfmt,grfir)
use modmain
implicit none
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: rfir(ngrtot)
real(8), intent(out) :: grfmt(lmmaxvr,nrmtmax,natmtot,3)
real(8), intent(out) :: grfir(ngrtot,3)
! local variables
integer is,ia,ias,i,ig,ifg
! allocatable arrays
real(8), allocatable :: grfmt1(:,:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(grfmt1(lmmaxvr,nrmtmax,3))
allocate(zfft1(ngrtot),zfft2(ngrtot))
! muffin-tin gradient
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call gradrfmt(lmaxvr,nrmt(is),spr(1,is),lmmaxvr,nrmtmax,rfmt(1,1,ias), &
     grfmt1)
    do i=1,3
      grfmt(:,1:nrmt(is),ias,i)=grfmt1(:,1:nrmt(is),i)
    end do
  end do
end do
! interstitial gradient
zfft1(:)=rfir(:)
call zfftifc(3,ngrid,-1,zfft1)
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
  end do
  call zfftifc(3,ngrid,1,zfft2)
  grfir(:,i)=dble(zfft2(:))
end do
deallocate(grfmt1,zfft1,zfft2)
return
end subroutine

