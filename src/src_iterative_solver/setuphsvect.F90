subroutine setuphsvect(n,hamilton,overlap,evecfv,h,s)
use modmain, only : nmatmax,nstfv
implicit none
integer ,intent(in):: n
complex, intent(in):: hamilton(n,n),overlap(n,n),evecfv(nmatmax,nstfv)
complex, intent(out)::h(n,nstfv),s(n,nstfv)
integer m 
m=nstfv
call zhemm('L','U',n,m,complex(1,0),hamilton,evecfv,nmatmax,h,n)
call zhemm('L','U',n,m,complex(1,0),overlap,evecfv,nmatmax,s,n)
end subroutine
