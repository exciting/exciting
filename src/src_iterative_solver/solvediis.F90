subroutine solvediis(m,Pmatrix,Qmatrix,c)
  use diisinterfaces
  implicit none

  integer, intent(in)::m

  real(8), intent(inout)::Pmatrix(m,m),Qmatrix(m,m)
  real(8), intent(out)::c(m)
  real(8):: work(8*m)
  real(8):: rwork(7*m),abstol,v
  integer:: iwork(5*m),ifail(m),info,mfound,lwork,i
  
  abstol=2.d0*dlamch('S')
  lwork =8*m
  i=1   
  call dsygvx(1,'V','I','U',m,Pmatrix,m,Qmatrix,m,&
       v,v,i,i,abstol,mfound,v,c,m,work,lwork, &
      iwork,ifail,info)
       
    if (info.ne.0) then
     write(*,*)
     write(*,'("Error(solvediis): diagonalisation failed")')
     write(*,'(" ZHEGVX returned INFO = ",I8)') info
     if (info.gt.m) then
        i=info-m
        write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
        write(*,'("  is not positive definite")')
        write(*,'(" Order of overlap matrix : ",I8)') m
        write(*,*) 
#ifdef DEBUG               
        write(775,*)(Pmatrix)
        write(776,*)(Qmatrix)
        stop
#endif
c=0.0
c(m)=1.0
     end if
     
  end if    
       
end subroutine solvediis
