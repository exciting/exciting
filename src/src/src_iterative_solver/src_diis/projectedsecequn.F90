
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine projectedsecequn (m, hp, op, evecp, evalp)
!
      Use modmain
      Implicit None
!update evecfv,eval  acording to sv
      Integer, Intent (In) :: m
      Complex (8), Intent (In) :: hp (2*m*(2*m+1)/2), op &
     & (2*m*(2*m+1)/2)
      Complex (8), Intent (Out) :: evecp (2*m, m)
      Real (8), Intent (Out) :: evalp (m)
      Complex (8) :: work (2*2*m)
      Real (8) :: abstol, rwork (7*2*m), vl, vu
      Integer :: nfound, iwork (5*2*m), ifail (2*m), info
      Real (8) :: dlamch
      External dlamch
!
      abstol = 20.d0 * dlamch ('S')
!
      Call zhpgvx (1, 'V', 'I', 'U', 2*m, hp, op, vl, vu, 1, m, abstol, &
     & nfound, evalp, evecp, 2*m, work, rwork, iwork, ifail, info)
#ifndef DEBUG2
      If (info .Gt. 0) Then
#endif	
         Write (*,*) "ifail: ", ifail, "\ninfo: ", info, "\nnfound: ", &
        & nfound
!
#ifndef DEBUG2	
         Stop
      End If
#endif
End Subroutine
