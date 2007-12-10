subroutine seceqfvprecond  (ik,n,h,o,X,evalfv,evecfv)
  use modmain
  integer, intent(in)::n,ik
  complex(8),intent(in)::h(n*(n+1)/2),o(n*(n+1)/2)
  complex(8),intent(out)::evecfv(nmatmax,nstfv)
  real(8),intent(out)::evalfv(nstfv)
  complex(8),intent (OUT)::X(nmatmax,nmatmax)	
!local var
 
!workarrays for lapaack
integer::info
integer  :: iwork(5*n)
integer  :: ifail(n)
real(8)  :: w(n)
real(8)  :: rwork(7*n)
complex(8):: work(2*n)

  call zhpgvx(1,'V','A','U',n,h,o,vl,vu,1,nstfv,abstol,m,w,X,nmatmax, &
       work,rwork,iwork,ifail,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnfv): diagonalisation failed")')
  write(*,'(" ZHPGVX returned INFO = ",I8)') info
  if (info.gt.nmatp) then
    i=info-nmatp
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') nmatp
    write(*,*)
  end if
  stop
end if

evalfv(1:nstfv)=w(1:nstfv)
evecfv(:,1:nstfv)=X(:,1:nstfv)
  call writeprecond(ik,n,X,w)
end subroutine seceqfvprecond
