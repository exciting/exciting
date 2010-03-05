
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine stepbound (reduction)
      Use modmixermsec, Only: SCHARGE, splane, dbase, qmx_input, qmx, &
     & qtot
      Implicit None
      Real (8), Intent (Out) :: reduction
      Real (8) :: limit, DSlope, PFACT
!
!       Simpler form
!       Set the limiting term based upon the maximum of
!               qtot:           The charge difference
!               splane:         The PW difference
!               Scharge:        The CLM difference
!               dbase:          Lower Bound
!
      Parameter (DSlope=2.0D0)! How much to reduce exponentially
      Parameter (PFACT=3.5D0)! Controls reduction in terms of limit
      limit = DSlope * Max (qtot, splane/PFACT)
      reduction = 0.1 + Exp (-limit)
      qmx = qmx_input * reduction
      If (qmx .Lt. dbase) qmx = dbase
      qmx = Min (qmx, qmx_input, 1.0D0)
!
      Return
End
