subroutine    mixmsec(iscl,potential,residualnorm,n)
use modmain,only:beta0,betainc,betadec
use modmixermsec,only:residual,last_outputp,last_inputp,initmixermsec
implicit none
integer, intent(in)::iscl,n
real(8), intent(inout)::potential(n) ! input/output potential
real(8), intent(in)::residualnorm    ! residual norm
!local variables
real(8),allocatable::S(:,:),Y(:,:),YY(:,:),STEP(:)
integer noldsteps,nwork
noldsteps=min(iscl-1,8)
if(iscl .eq. 1)then

call initmixermsec(n)
call mixadapt(iscl,beta0,betainc,betadec,n,potential,\
		last_inputp,last_outputp,residual,residualnorm)

!setup S,Y,YY,F
write(*,*) "n2:",n
else
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
