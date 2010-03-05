
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine setup_YY (iscl, n, S, Y, YY)
      Use modmixermsec, Only: noldstepsmax, PWHIST, FHIST, CLMHIST, &
     & YHIST, qmx, RedPred, RedOld, PM1, qmx_input, noldsteps, MUSE, &
     & MSECINFO, IDSCALE, residual, dmix_last, DIAG, dmixout, rtrap
      Use mod_Gvector, Only: ngrtot
      Implicit None
      Integer, Intent (In) :: iscl, n
      Real (8), Intent (In) :: S (n, noldstepsmax), Y (n, noldstepsmax)
      Real (8), Intent (Out) :: YY (noldstepsmax, noldstepsmax)
      Real (8) :: SS (noldstepsmax, noldstepsmax)
      Real (8), Parameter :: dbase = 0.005D0
      Real (8) :: ascl1, DMIXM, DMIX
      Integer :: k, j
!--------------------------------------------------------------------
!       Generate the noldstepsmax x noldstepsmax Matrices
!       Also generate scaling information
!--------------------------------------------------------------------
!
!       Note: should use dgemms here -- for later
!
      Do j = 1, noldstepsmax
         Do k = 1, j
            SS (j, k) = dot_product (S(1:n, j), S(1:n, k))
            YY (j, k) = dot_product (Y(1:n, j), Y(1:n, k))
         End Do
!
      End Do
!       Do transpose part as well
      Do j = 1, noldstepsmax
         Do k = 1, j - 1
            SS (k, j) = SS (j, k)
            YY (k, j) = YY (j, k)
         End Do
      End Do
      DIAG = Max (DIAG, 1D-12)
      DIAG = Min (DIAG, 0.1D0)
!
      RedPred = 1
        !experiment:
!
      MUSE = noldsteps
      Call LimitDMIX (Y, S, YY, residual, FHIST, YHIST, PWHIST, &
     & CLMHIST, n, noldsteps, noldstepsmax, ascl1, qmx_input, &
     & dmix_last, qmx, dmixout, IDSCALE, ngrtot, PM1, rtrap, dbase, &
     & iscl, RedPred, RedOld, DIAG, MUSE, MSECINFO)
!
!
      DMIXM = dmixout (1)
      DMIX = DMIXM
      qmx = DMIXM
      If (noldstepsmax .Lt. 10) Then
#ifdef DEBUG
         Write (*, 8001) MUSE, noldstepsmax, MSECINFO (1:4)
8001     Format (':DIRM :  MEMORY ', i1, '/', i1, ' RESCALE ', F6.3, ' &
        &RED ', F6.3, ' PRED ', F6.3, ' NEXT ', F6.3)
      Else
         Write (*, 8002) MUSE, noldstepsmax, MSECINFO (1:4)
8002     Format (':DIRM :  MEMORY ', i2, '/', i2, ' RESCALE ', F6.3, ' &
        &RED ', F6.3, ' PRED ', F6.3, ' NEXT ', F6.3)
#endif
      End If
!
End Subroutine
