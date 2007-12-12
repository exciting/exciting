subroutine calcupdatevectors(n,iunconverged,P,w,r,evalfv,phi) 
use modmain, only:nstfv,nmatmax
integer ,intent (in)::n , iunconverged
complex(8),intent(in)::P(nmatmax,nmatmax),r(n,nstfv)
complex(8),intent(out)::phi(n,nstfv)
complex(8):: v(n,nstfv)

integer m,i
m=nstfv
call zgemm('C','N',n,n,n,m,complex(1,0),P,nmatmax,r,n,complex(0,0),v,m)

do i=1,m
if(abs(w(i)-evalfv(i)).lt.1e-6)then
call zscal(n,complex(1.0/(w(i)-evalfv(i)),0),v,1)
else
v(:,i)=0
endif
end do

call zgemm('N','N',n,n,n,m,complex(1,0),P,nmatmax,r,n,complex(0,0),phi,m)

end subroutine
