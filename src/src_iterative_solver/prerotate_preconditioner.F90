subroutine prerotate_preconditioner(n,h,evecfv,X)

use modmain ,only:nstfv,nmatmax
implicit none

integer, intent(in) :: n
complex(8), intent(in)::h(n,n),evecfv(nstfv)
complex(8),intent(inout)::X(n,n)

complex(8):: hs(2*nstfv,2*nstfv) 
complex(8):: tmp(nmatmax,2*nstfv)
integer m,mfound,i
real(8):: v,c(2*nstfv,2*nstfv)
!work arrays
complex(8):: work(4*nstfv)
real(8):: rwork(14*nstfv),abstol
integer:: iwork(10*nstfv),ifail(2*nstfv),info
real(8) dlamch,eval(n)
external dlamch
abstol=2.d0*dlamch('S')
m=2*nstfv
call zhemm('L','U',n,m,complex(1,0),h,X,nmatmax,complex(0,0),tmp,nmatmax)
call zgemm('C','N',n,m,n,m,complex(1,0),X,nmatmax,tmp,nmatmax,&
     complex(0,0),hs,m)
call ZHEEVX( 'V', 'A', 'U',m , hs, m, v,v, i, i,&
                        ABSTOL, mfound, eval, c, m, WORK, 2*m, RWORK,&
                        IWORK, IFAIL, INFO )
tmp=0.0
call zgemm('N','N',n,m,m,m,complex(1,0),X,nmatmax,c,m,complex(0,0),tmp,nmatmax)
X(:,1:m)=tmp(:,1:m)
end subroutine
