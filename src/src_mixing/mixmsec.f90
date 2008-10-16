subroutine    mixmsec(iscl,potential,residualnorm,n)
	use modmain,only:beta0,betainc,betadec
	use modmixermsec,only:residual,last_outputp,last_inputp,initmixermsec,\
	freearraysmixermsec
	implicit none
	integer, intent(in)::iscl,n
	real(8), intent(inout)::potential(n) ! input/output potential
	real(8), intent(out)::residualnorm    ! residual norm
	!local variables
	real(8),allocatable::S(:,:),Y(:,:),YY(:,:),STEP(:)
	integer noldsteps,nwork
	integer, parameter::noldstepsmax=8
	real(8),parameter::DELTA=1e-3,DMIX=.5
	integer:: ifail
	noldsteps=min(iscl-1,noldstepsmax)
	if(iscl .eq. 1)then
		call mixadapt(iscl,beta0,betainc,betadec,n,potential,\
			last_inputp,last_outputp,residual,residualnorm)
	else
		allocate (S(n,noldstepsmax),Y(n,noldstepsmax),YY(noldstepsmax,noldstepsmax),STEP(n))
		residual=potential-last_outputp
		call appent_current_to_broyden_file(n,potential,residual)
		call readbroydsteps_and_init_SY(noldsteps,noldstepsmax,n,S,Y)
		call setup_YY(n,noldstepsmax,S,Y,YY)
		call MSEC1(Y,S,YY,residual,STEP,n,noldstepsmax,DMIX,IFAIL,DELTA,noldsteps)
		!          Y,S:            Conventional Y and S arrays
		!          YY:             Matrix of Y*Y values
		!          residual:       -Grad(MAXMIX) at the current point (residue)
		!		   n:              Length of the variable vector
		!          noldsteps :     Total number of memory values to use
		!                          Most recent is last
		!          DMIX:           Scaler for initial matrix
		!
		!          Output
		!          STEP            Multi-Secant Step
		potential=potential+STEP
		last_outputp=potential
		deallocate (S,Y,YY,STEP)
	endif
end subroutine
