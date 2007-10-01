subroutine dbxcplot
use modmain
implicit none
! local variables
integer idm,is,ia,ias,ir
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:,:)
real(8), allocatable :: rvfir(:,:)
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: rfir(:)
real(8), allocatable :: grfmt(:,:,:,:)
real(8), allocatable :: grfir(:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(dbxcplot): spin-unpolarised field is zero")')
  write(*,*)
  stop
end if
! initialise universal variables
call init0
! read magnetisation from file
call readstate
allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,3))
allocate(rvfir(ngrtot,3))
allocate(rfmt(lmmaxvr,nrmtmax,natmtot))
allocate(rfir(ngrtot))
allocate(grfmt(lmmaxvr,nrmtmax,natmtot,3))
allocate(grfir(ngrtot,3))
if (ndmag.eq.3) then
! non-collinear
  rvfmt(:,:,:,:)=bxcmt(:,:,:,:)
  rvfir(:,:)=bxcir(:,:)
else
! collinear
  rvfmt(:,:,:,1:2)=0.d0
  rvfir(:,1:2)=0.d0
  rvfmt(:,:,:,3)=bxcmt(:,:,:,1)
  rvfir(:,3)=bxcir(:,1)
end if
rfmt(:,:,:)=0.d0
rfir(:)=0.d0
do idm=1,3
  call gradrf(rvfmt(1,1,1,idm),rvfir(1,idm),grfmt,grfir)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is)
        rfmt(:,ir,ias)=rfmt(:,ir,ias)+grfmt(:,ir,ias,idm)
      end do
    end do
  end do
  rfir(:)=rfir(:)+grfir(:,idm)
end do
select case(task)
case(91)
  open(50,file='DBXC1D.OUT',action='WRITE',form='FORMATTED')
  open(51,file='DBXCLINES.OUT',action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,lmaxvr,lmmaxvr,rfmt,rfir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(dbxcplot):")')
  write(*,'(" 1D divergence of exchange-correlation field written to &
   &DBXC1D.OUT")')
  write(*,'(" vertex location lines written to DBXCLINES.OUT")')
case(92)
  open(50,file='DBXC2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,lmaxvr,lmmaxvr,rfmt,rfir)
  close(50)
  write(*,'("Info(dbxcplot):")')
  write(*,'(" 2D divergence of exchange-correlation field written to &
   &DBXC2D.OUT")')
case(93)
  open(50,file='DBXC3D.OUT',action='WRITE',form='FORMATTED')
  call plot3d(50,1,lmaxvr,lmmaxvr,rfmt,rfir)
  close(50)
  write(*,'("Info(dbxcplot):")')
  write(*,'(" 3D divergence of exchange-correlation field written to &
   &DBXC3D.OUT")')
end select
write(*,*)
deallocate(rvfmt,rvfir,rfmt,rfir,grfmt,grfir)
return
end subroutine

