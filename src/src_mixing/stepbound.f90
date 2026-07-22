
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine stepbound (reduction)
      Use modmixermsec, Only: SCHARGE, splane, dbase, qmx_input, qmx, &
     & qtot
      Use precision, Only: i32, dp
      Implicit None
      Real (dp), Intent (Out) :: reduction
      Real (dp) :: limit, DSlope, PFACT
!
!       Simpler form
!       Set the limiting term based upon the maximum of
!               qtot:           The charge difference
!               splane:         The PW difference
!               Scharge:        The CLM difference
!               dbase:          Lower Bound
!
      Parameter (DSlope=2.0_dp)! How much to reduce exponentially
      Parameter (PFACT=3.5_dp)! Controls reduction in terms of limit
      limit = DSlope * Max (qtot, splane/PFACT)
      reduction = 0.1_dp + Exp (-limit)
      qmx = qmx_input * reduction
      If (qmx .Lt. dbase) qmx = dbase
      qmx = Min (qmx, qmx_input, 1.0_dp)
!
      Return
End
