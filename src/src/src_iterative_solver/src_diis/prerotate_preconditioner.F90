
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine prerotate_preconditioner (n, m, h, P)
      Use modmain, Only: nstfv, nmatmax, zone, zzero
      Implicit None
      Integer, Intent (In) :: n, m
      Complex (8), Intent (In) :: h (n, n)
      Complex (8), Intent (Inout) :: P (nmatmax, nmatmax)
      Complex (8) :: hs (m, m)
      Complex (8) :: tmp (nmatmax, m), c (m, m)
      Integer :: mfound, i
      Real (8) :: v
  !work arrays
      Complex (8) :: work (2*m)
      Real (8) :: rwork (7*m), abstol
      Integer :: iwork (5*m), ifail (m), info
      Real (8) :: dlamch, eval (n)
      External dlamch
      abstol = 2.d0 * dlamch ('S')
!
!
!
      Call zhemm ('L', 'U', n, m, zone, h, nmatmax, P, nmatmax, zzero, &
     & tmp, nmatmax)
      Call zgemm ('C', 'N', m, m, n, zone, P, nmatmax, tmp, nmatmax, &
     & zzero, hs, m)
      Call ZHEEVX ('V', 'A', 'U', m, hs, m, v, v, i, i, abstol, mfound, &
     & eval, c, m, work, 2*m, rwork, iwork, ifail, info)
      Write (*,*) "prerotate zgemm n, m, nmatmax, mfound", n, m, &
     & nmatmax, mfound
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(prerotate): diagonalisation failed")')
         Write (*, '(" ZHEGVX returned INFO = ", I8)') info
         If (info .Gt. m) Then
            i = info - m
            Write (*, '(" The leading minor of the overlap matrix of or&
           &der ", I8)') i
            Write (*, '("  is not positive definite")')
            Write (*, '(" Order of overlap matrix : ", I8)') m
            Write (*,*)
#ifdef DEBUG		
            Write (775,*) hs
            Write (776,*) c
            Stop
#endif
         End If
      End If
      Call zgemm ('N', 'N', n, m, m, zone, P, nmatmax, c, m, zzero, &
     & tmp, nmatmax)
      Call zcopy (n*m, tmp(1, 1), 1, P(1, 1), 1)
!
!
End Subroutine prerotate_preconditioner
