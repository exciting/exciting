subroutine setuphsvect(n,m,hamilton,overlap,evecfv,h,s)
use modmain, only : nmatmax,nstfv
implicit none
integer ,intent(in):: n,m
complex(8), intent(in):: hamilton(n,n),overlap(n,n),evecfv(nmatmax,m)
complex(8), intent(out)::h(n,m),s(n,m)
call zhemm('L','U',n,m,complex(1,0),hamilton,evecfv,nmatmax,h,n)
call zhemm('L','U',n,m,complex(1,0),overlap,evecfv,nmatmax,s,n)
end subroutine
