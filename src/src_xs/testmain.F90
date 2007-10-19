
subroutine testmain
  use modmain
!  use m_ftfun
  implicit none




!!$  real(8), allocatable :: fmt(:,:,:)
!!$  real(8), allocatable :: fir(:)
!!$  complex(8),allocatable :: gft(:)
!!$  integer :: ir
!!$
!!$  call init0
!!$
!!$  allocate(fmt(lmmaxvr,nrmtmax,natmtot))
!!$  allocate(fir(ngrtot),gft(ngrtot))
!!$
!!$  fmt(:,:,:)=0.d0
!!$  fmt(1,:,:)=1.d0/y00
!!$  fir(:)=1.d0
!!$
!!$  call ftfun(.true.,.true.,fir,fmt,gft)
!!$
!!$  do ir=1,ngrtot
!!$     write(2000,'(3g18.10)') gft(ir), abs(gft(ir))
!!$  end do
!!$
!!$  write(*,*) 'normalization:',sum(abs(gft)**2)
!!$
!!$  deallocate(fir,fmt,gft)


end subroutine

