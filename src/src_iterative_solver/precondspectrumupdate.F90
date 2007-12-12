subroutine precondspectrumupdate(n,m,hamilton,overlap,P,w)
use modmain ,only:nstfv,nmatmax
implicit none
integer, intent(in) :: n,m
complex(8), intent(in)::hamilton(n,n),overlap(n,n)
complex(8),intent(in)::P(nmatmax,nmatmax)
complex(8),intent(inout)::w(nmatmax) 
complex(8)::h(n,m),s(n,m)
call setuphsvect(n,m,hamilton,overlap,P,h(:,:),s(:,:))
call rayleighqotient(n,m,P, h(:,:),s(:,:),w)
end subroutine