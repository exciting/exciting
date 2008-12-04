subroutine setuphsvect(n,m,system,evecfv,ldv,h,s)
  use diisinterfaces
  use  modfvsystem
  use modmain, only : nstfv,zone,zzero
  implicit none
  integer ,intent(in):: n,m,ldv
  complex(8), intent(inout):: evecfv(ldv,nstfv)
  type(evsystem)::system
  complex(8), intent(out)::h(n,nstfv),s(n,nstfv)
integer ::i
complex(8)::z
real(8)::t

  call zhemm('L','U',n,m,zone,system%overlap%za(1,1),n,evecfv(1,1),ldv,&
       zzero,s(1,1),n)
  Do i=1,m
  	t=(dble (zdotc(n,evecfv(1,i),1,s(1,i),1)))
    z=1.d0/sqrt(t)
   ! write(*,*)"z",z
  	call zscal(n,z,evecfv(1,i),1)
  	call zscal(n,z,s(1,i),1)
  end do
  h=h
   !call zhemm('L','U',n,m,zone,system%hamilton%za(1,1),n,evecfv(1,1),ldv,&
   !  zzero,h(1,1),n)

end subroutine setuphsvect
