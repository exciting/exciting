subroutine setuphsvect(n,m,hamilton,overlap,evecfv,h,s)
  use diisinterfaces
  use modmain, only : nstfv,zone,zzero
  implicit none
  integer ,intent(in):: n,m
  complex(8), intent(in):: hamilton(n,n),overlap(n,n),evecfv(n,m)
  complex(8), intent(out)::h(n,m),s(n,m)

  call zhemm('L','U',n,m,zone,hamilton(1,1),n,evecfv(1,1),n,&
     zzero,h(1,1),n)
  call zhemm('L','U',n,m,zone,overlap(1,1),n,evecfv(1,1),n,&
      zzero,s(1,1),n)
       


end subroutine setuphsvect
