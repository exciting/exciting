!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: polynom
! !INTERFACE:
Real (8) Function polynom (m, np, xa, ya, c, x)
! !INPUT/OUTPUT PARAMETERS:
!   m  : order of derivative (in,integer)
!   np : number of points to fit (in,integer)
!   xa : abscissa array (in,real(np))
!   ya : ordinate array (in,real(np))
!   c  : work array (out,real(np))
!   x  : evaluation abscissa (in,real)
! !DESCRIPTION:
!   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
!   function returns the $m$th derviative of the polynomial at $x$, while for
!   $m<0$ the integral of the polynomial from the first point in the array to
!   $x$ is returned.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! argmuments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: np
      Real (8), Intent (In) :: xa (np)
      Real (8), Intent (In) :: ya (np)
      Real (8), Intent (Out) :: c (np)
      Real (8), Intent (In) :: x
! local variables
      Integer :: i, j, k
      Real (8) :: x0, x1, x2, x3, y1, y2, y3
      Real (8) :: t1, t2, t3, t4, t5, t6, t7, sum
! fast evaluations for small np
      Select Case (np)
      Case (1)
         Select Case (m)
         Case (:-1)
            polynom = ya (1) * (x-xa(1))
         Case (0)
            polynom = ya (1)
         Case Default
            polynom = 0.d0
         End Select
         Return
      Case (2)
         c (2) = (ya(2)-ya(1)) / (xa(2)-xa(1))
         t1 = x - xa (1)
         Select Case (m)
         Case (:-1)
            polynom = t1 * (ya(1)+0.5d0*c(2)*t1)
         Case (0)
            polynom = c (2) * t1 + ya (1)
         Case (1)
            polynom = c (2)
         Case Default
            polynom = 0.d0
         End Select
         Return
      Case (3)
         x1 = xa (2) - xa (1)
         x2 = xa (3) - xa (1)
         y1 = ya (2) - ya (1)
         y2 = ya (3) - ya (1)
         t1 = 1.d0 / (x1*x2*(x2-x1))
         t2 = x1 * y2
         t3 = x2 * y1
         c (2) = t1 * (x2*t3-x1*t2)
         c (3) = t1 * (t2-t3)
         t1 = x - xa (1)
         Select Case (m)
         Case (:-1)
            polynom = t1 * &
           & (ya(1)+t1*(0.5d0*c(2)+0.3333333333333333333d0*c(3)*t1))
         Case (0)
            polynom = ya (1) + t1 * (c(2)+c(3)*t1)
         Case (1)
            polynom = 2.d0 * c (3) * t1 + c (2)
         Case (2)
            polynom = 2.d0 * c (3)
         Case Default
            polynom = 0.d0
         End Select
         Return
      Case (4)
         x1 = xa (2) - xa (1)
         x2 = xa (3) - xa (1)
         x3 = xa (4) - xa (1)
         y1 = ya (2) - ya (1)
         y2 = ya (3) - ya (1)
         y3 = ya (4) - ya (1)
         t1 = 1.d0 / (x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
         t2 = x1 * x2 * y3
         t3 = x2 * x3 * y1
         t4 = x3 * x1 * y2
         t5 = x1 ** 2
         t6 = x2 ** 2
         t7 = x3 ** 2
         c (2) = t1 * &
        & (t3*(x3*t6-x2*t7)+t4*(x1*t7-x3*t5)+t2*(x2*t5-x1*t6))
         c (3) = t1 * (t3*(t7-t6)+t4*(t5-t7)+t2*(t6-t5))
         c (4) = t1 * (t3*(x2-x3)+t4*(x3-x1)+t2*(x1-x2))
         t1 = x - xa (1)
         Select Case (m)
         Case (:-1)
            polynom = t1 * (ya(1)+t1*(0.5d0*c(2)+&
           & t1*(0.3333333333333333333d0*c(3)+0.25d0*c(4)*t1)))
         Case (0)
            polynom = ya (1) + t1 * (c(2)+t1*(c(3)+c(4)*t1))
         Case (1)
            polynom = c (2) + t1 * (2.d0*c(3)+3.d0*c(4)*t1)
         Case (2)
            polynom = 6.d0 * c (4) * t1 + 2.d0 * c (3)
         Case (3)
            polynom = 6.d0 * c (4)
         Case Default
            polynom = 0.d0
         End Select
         Return
      End Select
      If (np .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(polynom): np <= 0 : ", I8)') np
         Write (*,*)
         Stop
      End If
      If (m .Ge. np) Then
         polynom = 0.d0
         Return
      End If
! find the polynomial coefficients in divided differences form
      c (:) = ya (:)
      Do i = 2, np
         Do j = np, i, - 1
            c (j) = (c(j)-c(j-1)) / (xa(j)-xa(j+1-i))
         End Do
      End Do
! special case m=0
      If (m .Eq. 0) Then
         sum = c (1)
         t1 = 1.d0
         Do i = 2, np
            t1 = t1 * (x-xa(i-1))
            sum = sum + c (i) * t1
         End Do
         polynom = sum
         Return
      End If
      x0 = xa (1)
! convert to standard form
      Do j = 1, np - 1
         Do i = 1, np - j
            k = np - i
            c (k) = c (k) + (x0-xa(k-j+1)) * c (k+1)
         End Do
      End Do
      If (m .Gt. 0) Then
! take the m'th derivative
         Do j = 1, m
            Do i = m + 1, np
               c (i) = c (i) * dble (i-j)
            End Do
         End Do
         t1 = c (np)
         t2 = x - x0
         Do i = np - 1, m + 1, - 1
            t1 = t1 * t2 + c (i)
         End Do
         polynom = t1
      Else
! find the integral
         t1 = c (np) / dble (np)
         t2 = x - x0
         Do i = np - 1, 1, - 1
            t1 = t1 * t2 + c (i) / dble (i)
         End Do
         polynom = t1 * t2
      End If
      Return
End Function
!EOC
