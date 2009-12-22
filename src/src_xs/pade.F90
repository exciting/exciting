!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_pade
      Implicit None
Contains
!
!
      Subroutine pade (m, z, n, iw, ih, h)
    !
    ! Implementation of a Pade approximant using Thiele's method.
    ! Expressions taken from K. Lee, Phys. Rev. B 54 R8286 (1996)
    !
    ! Created Dec. 2006 (Sagmeister)
    !
         Use m_ctdfrac
         Implicit None
    ! arguments
         Integer, Intent (In) :: m
         Complex (8), Intent (In) :: z (m)
         Integer, Intent (In) :: n
         Complex (8), Intent (In) :: iw (n), ih (n)
         Complex (8), Intent (Out) :: h (m)
    ! local variables
         Character (*), Parameter :: thisnam = 'pade'
         Complex (8), Allocatable :: acoef (:), bcoef (:), a (:), aa &
        & (:), bb (:), c (:, :)
         Complex (8) :: zz
         Integer :: j, l, k
!
    ! require order higher than two
         If (n < 2) Then
            Write (*,*) 'Error(' // thisnam // '): approximant order to&
           &o small (< 2)'
            Call terminate
         End If
!
    ! allocate
         Allocate (acoef(n), bcoef(0:n), a(n), aa(0:n), bb(0:n), c(n, &
        & n))
!
    ! coefficients for numerator and denominator
         a (1) = ih (1)
         c (1, :) = ih (:)
         Do j = 2, n
            Do l = 2, n
               c (j, l) = (a(j-1)-c(j-1, l)) / ((iw(l)-iw(j-1))*c(j-1, &
              & l))
            End Do
            a (j) = c (j, j)
         End Do
!
    ! calculate values at given frequencies
         Do k = 1, m
            zz = z (k)
            bcoef (:) = (1.d0, 0.d0)
            bcoef (0) = (0.d0, 0.d0)
            acoef (1) = a (1)
            Do j = 2, n
               acoef (j) = a (j) * (zz-iw(j-1))
            End Do
       ! continued fraction evaluation of the Pade approximant
            Call ctdfrac (n, acoef, bcoef, h(k))
         End Do
!
         Deallocate (acoef, bcoef, a, aa, bb, c)
!
      End Subroutine
!
End Module
