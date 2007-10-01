subroutine dynsym(vpl,dynp)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(inout) :: dynp(3*natmtot,3*natmtot)
! local variables
logical tapp
integer isym,i,j,n
real(8) v(3),s(3,3),t1
! allocatable arrays
complex(8), allocatable :: dyns(:,:)
! external functions
real(8) r3taxi
external r3taxi
allocate(dyns(3*natmtot,3*natmtot))
n=0
dyns(:,:)=0.d0
! use the symmetries which leave vpl invariant
do isym=1,nsymlat
  s(:,:)=dble(symlat(:,:,isym))
  call r3mtv(s,vpl,v)
  if (r3taxi(vpl,v).lt.epslat) then
    call dynsymapp(symlat(1,1,isym),vpl,dynp,tapp,dyns)
    if (tapp) n=n+1
  end if
end do
if (n.eq.0) then
  write(*,*)
  write(*,'("Error(dynsym): no symmetries leave vpl invariant")')
  write(*,*)
  stop
end if
t1=1.d0/dble(n)
dynp(:,:)=t1*dyns(:,:)
! make the matrix Hermitian
do i=1,3*natmtot
  do j=i,3*natmtot
    dynp(i,j)=0.5d0*(dynp(i,j)+conjg(dynp(j,i)))
    dynp(j,i)=conjg(dynp(i,j))
  end do
end do
deallocate(dyns)
return
end subroutine

