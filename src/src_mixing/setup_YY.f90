



subroutine setup_YY(iscl, n, S, Y, YY)
use modmixermsec, only:noldstepsmax, PWHIST, FHIST, CLMHIST, &
YHIST, qmx, RedPred, RedOld, PM1, qmx_input, noldsteps, &
MUSE, MSECINFO, IDSCALE, residual, dmix_last, DIAG, dmixout, rtrap
use mod_Gvector, only:ngrtot
implicit none
integer, intent(in)::iscl, n
real(8), intent(in)::S(n, noldstepsmax), Y(n, noldstepsmax)
real(8), intent(out)::YY(noldstepsmax, noldstepsmax)
real(8)::SS(noldstepsmax, noldstepsmax)
real(8), parameter:: dbase=0.005D0
real(8)::ascl1, DMIXM,	 DMIX
integer::k, j
!--------------------------------------------------------------------
!       Generate the noldstepsmax x noldstepsmax Matrices
!       Also generate scaling information
!--------------------------------------------------------------------
!
!       Note: should use dgemms here -- for later
!
	DO J=1, noldstepsmax
	  DO K=1, J
		SS(J, K) =dot_product(S(1:n, J), S(1:n, K))
		YY(J, K) =dot_product(Y(1:n, J), Y(1:n, K))
	  enddo

	enddo
!       Do transpose part as well
	DO J=1, noldstepsmax
		DO K=1, J-1
			SS(K, J)=SS(J, K)
			YY(K, J)=YY(J, K)
		enddo
	enddo
	 DIAG=max(DIAG, 1D-12)
	DIAG=min(DIAG, 0.1D0)

	RedPred=1
        !experiment:

	MUSE=noldsteps
	 call LimitDMIX(Y, S, YY, residual, FHIST, YHIST, PWHIST, CLMHIST, n, noldsteps, noldstepsmax, &
			     ascl1, qmx_input, dmix_last, qmx, dmixout, &
			     IDSCALE, ngrtot, PM1, rtrap, Dbase, &
			     ISCL, RedPred, RedOld, DIAG, MUSE, MSECINFO)


       DMIXM=dmixout(1)
      DMIX=DMIXM
       qmx=dmixm
	if(noldstepsmax .lt. 10)then
#ifdef DEBUG
	write(*, 8001)MUSE, noldstepsmax, MSECINFO(1:4)
8001	format(':DIRM :  MEMORY ', i1, '/', i1, ' RESCALE ', F6.3, ' RED ', F6.3, ' PRED ', F6.3, ' NEXT ', F6.3)
	else
	write(*, 8002)MUSE, noldstepsmax, MSECINFO(1:4)
8002	format(':DIRM :  MEMORY ', i2, '/', i2, ' RESCALE ', F6.3, ' RED ', F6.3, ' PRED ', F6.3, ' NEXT ', F6.3)
#endif
	endif

end subroutine
