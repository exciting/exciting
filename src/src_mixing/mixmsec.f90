
! Copyright (C) 2002-2008 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixmsec
! !INTERFACE:
subroutine    mixmsec(iscl,potential,residualnorm,n)
!INPUT/OUTPUT PARAMETERS:
! iscl    	: self-consistent loop number (in,integer)
! potential	: potential coefitients packed in one array
! residualnorm: measure for convergence
! n			: length of mixing vector
!config params:
	use modmain,only:beta0,betainc,betadec
! persistent arrays and create/desdruct functions
	use modmixermsec,only:residual,last_outputp,last_inputp,initmixermsec,\
	freearraysmixermsec,noldstepsmax
	implicit none
	integer, intent(in)::iscl,n
	real(8), intent(inout)::potential(n) ! input/output potential
	real(8), intent(out)::residualnorm    ! residual norm
	!local variables
	real(8),allocatable::S(:,:),Y(:,:),YY(:,:),STEP(:)

	integer noldsteps,nwork
	real(8),parameter::DELTA=1e-3,DMIX=.5
	integer:: ifail
	noldsteps=min(iscl-2,noldstepsmax)
	if(iscl .le. 2)then
		if(iscl .eq. 2)then
			residual=potential-last_outputp
			call write_current_to_broyden_file(n,iscl,potential,residual)
		endif
		call mixadapt(iscl,beta0,betainc,betadec,n,potential,\
			last_outputp,last_inputp,residual,residualnorm)


	else
		allocate (S(n,noldstepsmax),Y(n,noldstepsmax),YY(noldstepsmax,noldstepsmax),STEP(n))
		residual=potential-last_outputp

		call readbroydsteps_and_init_SY(noldsteps,n,S,Y,potential,residual)
		call write_current_to_broyden_file(n,iscl,potential,residual)
		call rescaleYS(noldsteps,n,S,Y,potential,residual)
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
