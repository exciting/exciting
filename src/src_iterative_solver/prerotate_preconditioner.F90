subroutine prerotate_preconditioner(n,m,h,P)
  use modmain ,only:nstfv,nmatmax,zone,zzero
  implicit none
  integer, intent(in) :: n,m
  complex(8), intent(in)::h(n,n)
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

write(*,*) "prerotate zgemm n,m,nmatmax" ,n,m,nmatmax

  call zhemm('L','U',n,m,zone,h,nmatmax,P,nmatmax,zzero,tmp,nmatmax)
  call zgemm('C','N',m,m,n,zone,P,nmatmax,tmp,nmatmax,zzero,hs,m)
  call ZHEEVX( 'V', 'A', 'U',m , hs, m, v,v, i, i,&
       ABSTOL, mfound, eval, c, m, WORK, 2*m, RWORK,&
       IWORK, IFAIL, INFO )
           if (info.ne.0) then
     write(*,*)
     write(*,'("Error(prerotate): diagonalisation failed")')
     write(*,'(" ZHEGVX returned INFO = ",I8)') info
     if (info.gt.m) then
        i=info-m
        write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
        write(*,'("  is not positive definite")')
        write(*,'(" Order of overlap matrix : ",I8)') m
        write(*,*) 
#ifdef DEBUG               
        write(775,*)hs
        write(776,*)c
        stop
#endif
     end if
 end if
  call zgemm('N','N',n,m,m,zone,P,nmatmax,c,m,zzero,tmp,nmatmax)
  call zcopy(n*m,tmp(1,1),1,P(1,1),1)
 

end subroutine prerotate_preconditioner
