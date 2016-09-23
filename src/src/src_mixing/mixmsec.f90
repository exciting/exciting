
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixmsec
! !INTERFACE:
!
!
Subroutine mixmsec (iscl, potential, residualnorm, n)
!INPUT/OUTPUT PARAMETERS:
      Use modinput
! iscl    	: self-consistent loop number (in,integer)
! potential	: potential coefitients packed in one array
! residualnorm: measure for convergence
! n			: length of mixing vector
!config params:
      Use mod_potential_and_density, Only:
      Use mod_charge_and_moment, Only: chgir, chgmttot, chgtot
! persistent arrays and create/desdruct functions
      Use modmixermsec, Only: residual, last_outputp, work2, work3, &
     & initmixermsec, freearraysmixermsec, noldstepsmax, &
     & noldstepsin_file, noldsteps, qmx, dmix, dmixout, TCharge, &
     & SCharge, splane, tplane, qmx_input, qtot
!EOP
!BOC
      Implicit None
      Integer, Intent (In) :: iscl, n
      Real (8), Intent (Inout) :: potential (n)! input/output potential
      Real (8), Intent (Out) :: residualnorm ! residual norm
	!local variables
      Real (8), Allocatable :: S (:, :), Y (:, :), YY (:, :), &
     & broydenstep (:)
      Real (8), Parameter :: DELTA = 1e-3
      Integer :: ifail
      Real (8) :: sreduction, dmixm
      Real (8), External :: dnrm2
!
!
      noldsteps = noldstepsin_file
      sreduction = 1.2
!      write(*,*) iscl
      If (iscl .Le. input%groundstate%PrelimLinSteps) Then
!
         If (iscl .Ge. 0) Then
!          if (iscl.eq.0) last_outputp=0d0
            residual(1:n) = potential - last_outputp(1:n)
            Call write_current_to_broyden_file (n, iscl, potential, &
           & residual)
         End If
         Call mixadapt (iscl, input%groundstate%beta0, &
        & input%groundstate%betainc, input%groundstate%betadec, n, &
        & potential, last_outputp, work3, work2, residualnorm)
         last_outputp(1:n) = potential
         If (iscl .Eq. input%groundstate%PrelimLinSteps .And. allocated(work2) .And. allocated(work3)) &
        & deallocate (work2, work3)
      Else
         Allocate (S(n, noldstepsmax), Y(n, noldstepsmax))
         Allocate (YY(noldstepsmax, noldstepsmax))
         Allocate (broydenstep(n))
         residual(1:n) = potential - last_outputp(1:n)
         SCharge = chgir
         TCharge = chgtot
         Call check_msecparameters ()
         Call readbroydsteps_and_init_SY (noldsteps, n, S, Y, &
        & potential, residual)
         Call write_current_to_broyden_file (n, iscl, potential, &
        & residual)
         !write(*,210)':PLANE:  INTERSTITIAL TOTAL ',Tplane, ' DISTAN ',Splane
         ! write(*,210)':CHARG:  CLM CHARGE   TOTAL ',TCharge,' DISTAN ',SCharge
         Call stepbound (sreduction)
         !Write (60, 4141) sreduction, qmx
         Call rescaleYS (noldsteps, n, S, Y, potential, residual)
         Call setup_YY (iscl, n, S, Y, YY)
!
         dmixm = 0.1
!
!
!
!
         Call MSEC1 (Y, S, YY, residual, broydenstep, n, noldstepsmax, &
        & dmixm, ifail, DELTA, noldsteps)
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
!
         If (ifail .Ne. 0) Then
            Write (*,*) ':WARNING: Inversion of Multi-Secant Matrix Fa&
           &iled'
!
            Stop
         End If
       !call stepbound(sreduction)
!
         potential = potential + broydenstep
         last_outputp(1:n) = potential
!
!
         Deallocate (S, Y, YY, broydenstep)
4141     Format ('REDuction and DMIX in Broyd:', 3 f10.4, E14.5)
         residualnorm = dnrm2 (n, residual, 1)
 !after /sqrt(n) its not the residual norm anny more
 !nore is it residual mean square but thats how it is in mixadapt
 !
         residualnorm = residualnorm / Sqrt (dble(n))
      End If
!
      qtot = residualnorm
End Subroutine
!EOC
