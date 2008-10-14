subroutine    mixmsec(iscl,potential,residual,n)
use modmain,only:beta0,betainc,betadec
implicit none
integer, intent(in)::iscl,n
real(8), intent(inout)::potential(n)
real(8), intent(in)::residual(n)
!local variables
real(8),allocatable::S(:,:),Y(:,:),YY(:,:),STEP(:),work(:)
integer noldsteps,nwork
noldsteps=min(iscl-1,8)

if(iscl .gt. 2) call readbroydsteps(work,residual,n,noldsteps)
call appent_current_to_broyden_file(work,residual,potential,residual,n,noldsteps)
if (iscl.lt.2) then
allocate(work(3*n))
nwork=3*n
write(*,*) "n:",n
!call mixadapt(iscl,beta0,betainc,betadec,n,potential,work,work(n+1),work(2*n+1),residual)
deallocate(work)
else

!setup S,Y,YY,F
write(*,*) "n2:",n
allocate (S(n,n),Y(n,n),YY(n,n),STEP(n))
!rescale things

!call  MSEC1(Y,S,YY,F,STEP,MAXMIX,MEMORY,DMIX,IFAIL,DELTA,MUSE)
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

potential=potential+STEP
deallocate (S,Y,YY,STEP)
endif

end subroutine
