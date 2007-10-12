!BOP
! !ROUTINE: seceqn
subroutine residualvector(n,np,HeS,evecfv,r,rnorm)

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
use modmain
integer , intent (in)::n,np
complex(8),intent(in)::HeS(np) !packed ut
complex(8),intent(in)::evecfv(nmatmax) !vector
complex(8),intent(out)::r(n)
real(8),intent(out)::rnorm
real(8) zdotc
external zdotc

r(:)=0.0
call zhpmv("U",n,(1,1),HeS(:),evecfv(:), 1, (0,0), r(:), 1)
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
