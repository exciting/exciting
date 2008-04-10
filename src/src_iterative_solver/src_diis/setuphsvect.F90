subroutine setuphsvect(n,m,system,evecfv,ldv,h,s)
  use diisinterfaces
  use  modfvsystem
  use modmain, only : nstfv,zone,zzero
  implicit none
  integer ,intent(in):: n,m,ldv
  complex(8), intent(in):: evecfv(ldv,nstfv)
  type(evsystem)::system
  complex(8), intent(out)::h(n,nstfv),s(n,nstfv)

  call zhemm('L','U',n,m,zone,system%hamilton%za(1,1),n,evecfv(1,1),ldv,&
     zzero,h(1,1),n)
  call zhemm('L','U',n,m,zone,system%overlap%za(1,1),n,evecfv(1,1),ldv,&
      zzero,s(1,1),n)
    
end subroutine setuphsvect
