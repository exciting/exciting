!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine genveelu (l, u, j, lmax, vee)
      Implicit None
! arguments
      Integer, Intent (In) :: l
      Real (8), Intent (In) :: u
      Real (8), Intent (In) :: j
      Integer, Intent (In) :: lmax
      Real (8), Intent (Out) :: vee &
     & (-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax)
! local variables
      Integer :: m1, m2, m3, m4, k, q
      Real (8), Parameter :: fourpi = 12.566370614359172954d0
      Real (8) :: r1, r2, f (0:6)
      Real (8) :: sum1, sum2, t1
! external functions
      Real (8) :: gaunt
      External gaunt
! Slater integrals F(k) for d and f electrons in Ry, to be converted in Htr
      f (:) = 0.d0
      f (0) = u
      Select Case (l)
      Case (0)
! s electrons only f(0)=u
      Case (1)
! p electrons
         f (2) = 5.d0 * j
      Case (2)
! d electrons: ratio r1 = F(4)/F(2), see PRB 52, R5467 (1995)
         r1 = 0.625d0
         f (2) = (14.d0*j) / (1.d0+r1)
         f (4) = f (2) * r1
      Case (3)
! f electrons: r2 = F(6)/F(2), r1 = F(4)/F(2), see PRB 50, 16861 (1994)
         r1 = 451.d0 / 675.d0
         r2 = 1001.d0 / 2025.d0
         f (2) = 6435.d0 * j / (286.d0+195.d0*r1+250.d0*r2)
         f (4) = f (2) * r1
         f (6) = f (2) * r2
      Case Default
         Write (*,*)
         Write (*, '("Error(genveelu): invalid l : ", I8)') l
         Write (*,*)
         Stop
      End Select
      Do m1 = - l, l
         Do m2 = - l, l
            Do m3 = - l, l
               Do m4 = - l, l
                  sum1 = 0.d0
                  Do k = 0, 2 * l, 2
                     sum2 = 0.d0
                     Do q = - k, k
                        t1 = gaunt (l, k, l, m1, q, m2) * gaunt (l, k, &
                       & l, m3,-q, m4)
                        If (Mod(q, 2) .Eq. 0) Then
                           sum2 = sum2 + t1
                        Else
                           sum2 = sum2 - t1
                        End If
                     End Do
                     sum1 = sum1 + f (k) * sum2 / dble (2*k+1)
                  End Do
                  vee (m1, m3, m2, m4) = fourpi * sum1
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
