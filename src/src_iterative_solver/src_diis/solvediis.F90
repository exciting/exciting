
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine solvediis (m, Pmatrix, Qmatrix, c)
      Use diisinterfaces
      Implicit None
      Integer, Intent (In) :: m
!
      Real (8), Intent (Inout) :: Pmatrix (m+1, m+1), Qmatrix (m+1, &
     & m+1)
      Real (8), Intent (Out) :: c (m+1)
      Real (8) :: work (8*m)
      Real (8) :: rwork (7*m), abstol, v
      Integer :: iwork (5*m), ifail (m), info, mfound, lwork, i
!
      abstol = 2.d0 * dlamch ('S')
      lwork = 8 * m
      i = 1
      Call dsygvx (1, 'V', 'I', 'U', m, Pmatrix, m+1, Qmatrix, m+1, v, &
     & v, i, i, abstol, mfound, v, c, m+1, work, lwork, iwork, ifail, &
     & info)
!
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(solvediis): diagonalisation failed")')
         Write (*, '(" ZHEGVX returned INFO = ", I8)') info
         If (info .Gt. m) Then
            i = info - m
            Write (*, '(" The leading minor of the overlap matrix of or&
           &der ", I8)') i
            Write (*, '("  is not positive definite")')
            Write (*, '(" Order of overlap matrix : ", I8)') m
            Write (*,*)
#ifdef DEBUG		
            Write (775,*) (Pmatrix)
            Write (776,*) (Qmatrix)
            Stop
#endif
            c = 0.0
            c (m) = 1.0
         End If
!
      End If
!
End Subroutine solvediis
