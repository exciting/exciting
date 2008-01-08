subroutine prerotate_preconditioner(n,m,h,evecfv,P)
  use modmain ,only:nstfv,nmatmax,zone,zzero
  implicit none
  integer, intent(in) :: n,m
  complex(8), intent(in)::h(n,n),evecfv(nstfv)
  complex(8),intent(inout)::P(nmatmax,nmatmax)
  complex(8):: hs(m,m) 
  complex(8):: tmp(nmatmax,m),c(m,m)
  integer:: mfound,i
  real(8):: v
  !work arrays
  complex(8):: work(2*m)
  real(8):: rwork(7*m),abstol
  integer:: iwork(5*m),ifail(m),info
  real(8) dlamch,eval(n)
  external dlamch
  abstol=2.d0*dlamch('S')
#ifdef DEBUG
write(*,*) "prerotate zgemm"
#endif
  call zhemm('L','U',n,m,zone,h,P,nmatmax,zzero,tmp,nmatmax)
  call zgemm('C','N',n,m,n,m,zone,P,nmatmax,tmp,nmatmax,&
       cmplx(0,0),hs,m)
  call ZHEEVX( 'V', 'A', 'U',m , hs, m, v,v, i, i,&
       ABSTOL, mfound, eval, c, m, WORK, 2*m, RWORK,&
       IWORK, IFAIL, INFO )
 
 ! call zgemm('N','N',n,m,m,m,zone,P,nmatmax,c,m,zzero,tmp,nmatmax)
  
  !P(:,1:m)=tmp(:,1:m)
  !zcopy!
end subroutine prerotate_preconditioner
