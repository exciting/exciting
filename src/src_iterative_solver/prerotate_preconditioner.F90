subroutine prerotate_preconditioner(n,h,o,evecfv,X,w)
  use diisinterfaces
  use modmain ,only:nstfv,nmatmax
  implicit none

  integer, intent(in) :: n
  complex(8), intent(in)::h(n,n),o(n,n),evecfv(nstfv)
  complex(8),intent(inout)::X(nmatmax,nmatmax)
  real(8),intent(inout)::w(n)
  complex(8) zdotc
  external zdotc
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

  call zhemm('L','U',n,m,complex(1,0),h,X,nmatmax,complex(0,0),tmp,nmatmax)
  do i=1,m
     c(i,1)=zdotc(n,X(1,i),1,tmp(1,i),1)
  end do
  call zhemm('L','U',n,m,complex(1,0),o,X,nmatmax,complex(0,0),tmp,nmatmax)
  do i=1,m
     c(i,2)=zdotc(n,X(1,i),1,tmp(1,i),1)
  end do
  w(1:m)=c(1:m,1)/c(1:m,2)

end subroutine prerotate_preconditioner
