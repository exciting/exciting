subroutine setuphsvect(n,m,hamilton,overlap,evecfv,h,s)
  use modmain, only : nmatmax,nstfv
  implicit none
  integer ,intent(in):: n,m
  complex(8), intent(in):: hamilton(n,n),overlap(n,n),evecfv(nmatmax,m)
  complex(8), intent(out)::h(n,m),s(n,m)

  call zhemm('L','U',n,m,cmplx(1,0),hamilton,n,evecfv,nmatmax,&
       cmplx(0,0),h,n)
  call zhemm('L','U',n,m,cmplx(1,0),overlap,n,evecfv,nmatmax,&
       cmplx(0,0),s,n)

end subroutine setuphsvect
