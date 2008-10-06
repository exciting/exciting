subroutine    mixmsec(iscl,v,dv,n)
use modmain,only:
integer, intent(in)::iscl,n
real(8), intent(inout)::v(n)
!local variables
real(8),allocatable::S(:,:),Y(:,:),YY(:,:),STEP(:)
integer noldsteps
noldsteps=min(iscl-1,8)

if(iscl .gt. 2) call readbroydsteps(Vprev,Residprev,n,noldsteps)
call appent_current_to_broyden_file(Vprev,Residprev,v,dv,n,noldsteps)
if (iscl.eq.2) then
call MixPratt()
else

!setup S,Y,YY,F

!rescale things

call  MSEC1(Y,S,YY,F,STEP,MAXMIX,MEMORY,DMIX,IFAIL,DELTA,MUSE)
!          Y,S:            Conventional Y and S arrays
!          YY:             Matrix of Y*Y values
!          F:              -Grad(MAXMIX) at the current point (residue)
!		   MAXMIX:         Length of the variable vector
!          MEMORY:         Total number of memory values to use
!                          Most recent is last
!          DMIX:           Scaler for initial matrix
!
!          Output
!          STEP            Multi-Secant Step

v=v+STEP
endif

end subroutine
