subroutine seceqfvprecond  (n,system,X,w,evalfv,evecfv)
  use modmain, only: nmatmax,nstfv
  use  modfvsystem
  use diisinterfaces
  implicit none
  integer, intent(in)::n
type(evsystem)::system
  complex(8),intent(out)::evecfv(nmatmax,nstfv)
  real(8),intent(out)::evalfv(nstfv), w(nmatmax)
  complex(8),intent (OUT)::X(nmatmax,nmatmax)	

  !local var
  integer i
  !workarrays for lapaack
  integer::info,lwork
  integer  :: iwork(5*n)
  integer  :: ifail(n),mfound
  real(8)  ::v,abstol
  real(8)  :: rwork(7*n)
  complex(8):: work(2*n)
  lwork=2*n
   abstol=2.d0*dlamch('S')
  call zhegvx(1,'V','A','U',n,system%hamilton%za,n,system%overlap%za,n,v,v,1,nstfv,abstol,mfound,w,X,nmatmax, &
       work,lwork,rwork,iwork,ifail,info)
  if (info.ne.0) then
     write(*,*)
     write(*,'("Error(seceqnfv): diagonalisation failed")')
     write(*,'(" ZHPGVX returned INFO = ",I8)') info
     if (info.gt.n) then
        i=info-n
        write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
        write(*,'("  is not positive definite")')
        write(*,'(" Order of overlap matrix : ",I8)') n
        write(*,*)
     end if
     stop
  end if
  call dcopy(nstfv,w(1),1,evalfv(1),1)
  call zcopy(nstfv*nmatmax,X,1,evecfv(1,1),1)


end subroutine seceqfvprecond
