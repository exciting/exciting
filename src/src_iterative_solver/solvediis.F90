subroutine solvediis(m,Pmatrix,Qmatrix,c)
  use diisinterfaces
  implicit none

  integer, intent(in)::m

  complex(8), intent(in)::Pmatrix(m,m),Qmatrix(m,m)
  complex(8), intent(out)::c(m)
  complex(8):: work(2*m)
  real(8):: rwork(7*m),abstol,v
  integer:: iwork(5*m),ifail(m),info,mfound,lwork
  real(8) dlamch 
  external dlamch
  abstol=2.d0*dlamch('S')
  lwork =2*m
  call zhegvx(1,'V','I','U',m,Pmatrix,m,Qmatrix,m,&
       v,v,1,1,abstol,mfound,v,c,m,work,lwork, &
       rwork,iwork,ifail,info)
end subroutine solvediis
