
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
	use modmixermsec,only:residual,last_outputp,last_inputp,initmixermsec,&
freearraysmixermsec,noldstepsmax,noldstepsin_file,&
noldsteps,qmx

	implicit none
	integer, intent(in)::iscl,n
	real(8), intent(inout)::potential(n) ! input/output potential
	real(8), intent(out)::residualnorm    ! residual norm
	!local variables
	real(8),allocatable::S(:,:),Y(:,:),YY(:,:),broydenstep(:)
	real(8),parameter::DELTA=1e-3,DMIX=.5
	integer:: ifail
	real sreduction
	noldsteps=noldstepsin_file
	if(iscl .le. 2)then
		if(iscl .eq. 2)then
			residual=potential-last_outputp
			call write_current_to_broyden_file(n,iscl,potential,residual)
		endif
		call mixadapt(iscl,beta0,betainc,betadec,n,potential,\
			last_outputp,last_inputp,residual,residualnorm)


	else
		allocate (S(n,noldstepsmax),Y(n,noldstepsmax))
		allocate(YY(noldstepsmax,noldstepsmax))
		allocate(broydenstep(n))
		residual=potential-last_outputp
		call check_msecparameters()
		call readbroydsteps_and_init_SY(noldsteps,n,S,Y,potential,residual)
		call write_current_to_broyden_file(n,iscl,potential,residual)
	    call stepbound(sreduction)
        write(21,4141)sreduction,qmx
		call rescaleYS(noldsteps,n,S,Y,potential,residual)
		call setup_YY(iscl,n,S,Y,YY)
		call MSEC1(Y,S,YY,residual,broydenstep,n,noldstepsmax,DMIX,IFAIL,DELTA,noldsteps)
		!          Y,S:            Conventional Y and S arrays
		!          YY:             Matrix of Y*Y values
		!          residual:       -Grad(MAXMIX) at the current point (residue)
		!		   n:              Length of the variable vector
		!          noldsteps :     Total number of memory values to use
		!                          Most recent is last
		!          DMIX:           Scaler for initial matrix
		!
		!          Output
		!          broydenstep            Multi-Secant Step
	   ! call stepbound(sreduction)

		potential(1:n)=potential(1:n)+broydenstep(1:n)
		last_outputp=potential

if(.not. allocated(broydenstep))write(*,*)"errore malloc broydenstep"
		deallocate (S,Y,YY,broydenstep)
4141    format(':REDuction and DMIX in Broyd:',3f10.4,E14.5)

	endif

end subroutine

