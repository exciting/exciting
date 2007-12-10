!BOP
! !ROUTINE: seceqn
subroutine residualvector(n,h,o,evecfv,evalfv,r,rnorm)
use modmain, only: nmatmax
	
! !INPUT/OUTPUT PARAMETERS:

! !DESCRIPTION:
!  the residual is
!  \newcommand{\bra}[1]{\langle #1|}
!\newcommand{\ket}[1]{|#1\rangle}
!\newcommand{\braket}[2]{\langle #1|#2\rangle}
!$$
!\ket{\mathbf{R}\left(\ket{\mathbf{A}^{ap}},E^{ap}\right)}=(\mathbf{H}-E^{ap}\mathbf{S})\ket{ \mathbf{A}^{ap}}
!$$
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
!use modmain
implicit none
integer , intent (in)::n
!packed ut
complex(8),intent(in)::evecfv(nmatmax),h(n*(n+1)/2),o(n*(n+1)/2) !vector
complex(8),intent(out)::r(n)
real(8),intent(in)::evalfv
real(8),intent(out)::rnorm
complex(8) zdotc
external zdotc
integer:: i ,np
complex(8)::HeS(n*(n+1)/2) 
np=n*(n+1)/2
  ! blas call means : HminuseS(:)=h(:)-evalfv(ievec,ispn)*o(:)
           call zcopy(np,h,1,hes,1)
           call zaxpy(np,dcmplx(-evalfv, 0),o,1,hes,1)
r(:)=0.0
call zhpmv("U",n,dcmplx(1.0,0.0),HeS,evecfv, 1, dcmplx(0,0), r(1), 1)
#ifdef DEBUG
write(441,*)"Hes in residualvector",HeS
write(442,*)"evecfv in residualvector",evecfv
write(443,*)"n,np,r in residualvector",n,np,r

#endif
rnorm=0
do i=1,n
rnorm=rnorm+conjg(r(i))*r(i)
end do
rnorm=sqrt(rnorm)
end subroutine
