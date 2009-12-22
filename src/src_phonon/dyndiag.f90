!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dyndiag (dynp, w, ev)
      Use modmain
      Implicit None
! arguments
      Complex (8), Intent (In) :: dynp (3*natmtot, 3*natmtot)
      Real (8), Intent (Out) :: w (3*natmtot)
      Complex (8), Intent (Out) :: ev (3*natmtot, 3*natmtot)
! local variables
      Integer :: is, ia, ip, js, ja, jp
      Integer :: i, j, n, lwork, info
      Real (8) :: t1
! allocatable arrays
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: work (:)
! number of phonon branches
      n = 3 * natmtot
      ev (:, :) = 0.d0
      i = 0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Do ip = 1, 3
               i = i + 1
               j = 0
               Do js = 1, nspecies
! mass factor
                  If ((spmass(is) .Le. 0.d0) .Or. (spmass(js) .Le. &
                 & 0.d0)) Then
! infinite mass
                     t1 = 0.d0
                  Else
                     t1 = 1.d0 / Sqrt (spmass(is)*spmass(js))
                  End If
                  Do ja = 1, natoms (js)
                     Do jp = 1, 3
                        j = j + 1
                        If (i .Le. j) Then
! use Hermitian average of dynamical matrix
                           ev (i, j) = 0.5d0 * t1 * (dynp(i, &
                          & j)+conjg(dynp(j, i)))
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      Allocate (rwork(3*n))
      lwork = 2 * n
      Allocate (work(lwork))
      Call zheev ('V', 'U', n, ev, n, w, work, lwork, rwork, info)
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(dyndiag): diagonalisation failed")')
         Write (*, '(" ZHEEV returned INFO = ", I8)') info
         Write (*,*)
         Stop
      End If
      Do i = 1, n
         If (w(i) .Ge. 0.d0) Then
            w (i) = Sqrt (w(i))
         Else
            w (i) = - Sqrt (Abs(w(i)))
         End If
      End Do
      Deallocate (rwork, work)
      Return
End Subroutine
