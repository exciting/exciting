!BOP
! !ROUTINE: seceqn
subroutine residualvectors(n,iunconverged,h,s,evalfv,r,rnorms)
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
use modmain, only:nstfv
implicit none
integer , intent (in)::n,iunconverged
!packed ut
complex(8),intent(in)::h(n,nstfv),s(n,nstfv) 
complex(8),intent(out)::r(n,nstfv)
real(8),intent(in)::evalfv(nstfv)
real(8),intent(out)::rnorms(nstfv)
integer i
complex(8) zdotc
external zdotc

do i=1,nstfv
call zcopy(n,h(1,i),1,r(1,i),1)
call zaxpy(n,complex(-evalfv(i),0),s(1,i),1,r(1,i),1)
rnorms(i)= sqrt( zdotc(n,r(1,i),1,r(1,i),1))
end do

end subroutine
