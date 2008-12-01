
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
	use mod_potential_and_density,only:beta0,betainc,betadec
	use mod_charge_and_moment,only:chgir,chgmttot,chgtot
! persistent arrays and create/desdruct functions
	use modmixermsec,only:residual,last_outputp,work2,work3,initmixermsec,&
freearraysmixermsec,noldstepsmax,noldstepsin_file,&
noldsteps,qmx,dmix,dmixout,TCharge,SCharge,splane,tplane,  qmx_input,qtot

	implicit none
	integer, intent(in)::iscl,n
	real(8), intent(inout)::potential(n) ! input/output potential
	real(8), intent(out)::residualnorm    ! residual norm
	!local variables
	real(8),allocatable::S(:,:),Y(:,:),YY(:,:),broydenstep(:)
	real(8),parameter::DELTA=1e-3
	integer:: ifail
	real(8) sreduction,dmixm
	real(8),external::dnrm2


	noldsteps=noldstepsin_file
	sreduction=1.2
	if(iscl .le. 2)then

		if(iscl .ge. 2)then
			residual=potential-last_outputp
			call write_current_to_broyden_file(n,iscl,potential,residual)
		endif

		call mixadapt(iscl,beta0,betainc,betadec,n,potential,&
			last_outputp,work3,work2,residualnorm)
		last_outputp=potential
		if(iscl.eq.2 .and. allocated(work2) .and.allocated(work3))&
			 deallocate(work2,work3)
	else
		allocate (S(n,noldstepsmax),Y(n,noldstepsmax))
		allocate(YY(noldstepsmax,noldstepsmax))
		allocate(broydenstep(n))
		residual=potential-last_outputp
		SCHARGE=chgir
		TCharge= chgtot
		call check_msecparameters()
		call readbroydsteps_and_init_SY(noldsteps,n,S,Y,potential,residual)
		call write_current_to_broyden_file(n,iscl,potential,residual)
		!write(*,210)':PLANE:  INTERSTITIAL TOTAL ',Tplane, ' DISTAN ',Splane
     	! write(*,210)':CHARG:  CLM CHARGE   TOTAL ',TCharge,' DISTAN ',SCharge
210 	FORMAT(A,F12.5,A,F11.7)
	    call stepbound(sreduction)
        write(60,4141)sreduction,qmx
		call rescaleYS(noldsteps,n,S,Y,potential,residual)
		call setup_YY(iscl,n,S,Y,YY)

DMIXM=0.1




		call MSEC1(Y,S,YY,residual,broydenstep,n,noldstepsmax,DMIXM,IFAIL,DELTA,noldsteps)
		!          Y,S:            Conventional Y and S arrays
		!          YY:             Matrix of Y*Y values
		!          residual:       -Grad(MAXMIX) at the current point (residue)
		!		   n:              Length of the variable vector
		!          noldsteps :     Total number of memory values to use
		!                          Most recent is last
		!          DMIXM:           Scaler for initial matrix
		!
		!          Output
		!          broydenstep            Multi-Secant Step
		!call MSEC2(Y,S,YY,residual,broydenstep,n,noldstepsmax,DMIXM,IFAIL,DELTA)

	   if(IFAIL .ne. 0)then
                write(21,*)':WARNING: Inversion of Multi-Secant Matrix Failed'

              stop
       endif
       !call stepbound(sreduction)

		potential=potential+broydenstep
		last_outputp=potential


		deallocate (S,Y,YY,broydenstep)
4141    format('REDuction and DMIX in Broyd:',3f10.4,E14.5)
residualnorm=dnrm2(n,residual,1)
 !after /sqrt(n) its not the residual norm anny more
 !nore is it residual mean square but thats how it is in mixadapt
 !
 residualnorm=residualnorm/sqrt(dble(n))
	endif

 qtot= residualnorm
end subroutine

