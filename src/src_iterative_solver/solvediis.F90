subroutine solvediis(m,Pmatrix,Qmatrix,c)
  use diisinterfaces
implicit none

integer, intent(in)::m

complex(8), intent(in)::Pmatrix(m,m),Qmatrix(m,m)
complex(8), intent(out)::c(m)
complex(8):: work(2*m)
real(8):: rwork(7*m),abstol,v
integer:: iwork(5*m),ifail(m),info,mfound
real(8) dlamch ,eval(m)
external dlamch
abstol=2.d0*dlamch('S')
call zhpgvx(1,'V','I','U',2*m,Pmatrix,Qmatrix,v,v,1,1,abstol,mfound,eval,c,m,work, &
 rwork,iwork,ifail,info)
end subroutine