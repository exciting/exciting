
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Module m_pade
      Implicit None
Contains

!BOP
! !ROUTINE: pade
! !INTERFACE:
      Subroutine pade (m, z, n, iw, ih, h)
! !INPUT/OUTPUT PARAMETERS:
!   m     : number of values in array below (in,integer)
!   z     : values for which analytic continuation is to be performed
!           (in,complex(m))
!   n     : number of values in arrays below (in,integer)
!   iw    : values for which function is evaluated
!           (in,complex(n))
!   ih    : function evaluated at values "iw" (in,complex(n))
!   h     : function evaluated at values "z" (out,complex(n))
!
! !USES:
         Use m_ctdfrac
! !DESCRIPTION:
!   Implementation of a Pade approximant using Thiele's reciprocal-difference method.
!   This routine takes a complex function evaluated at an initial set of arguments,
!   approximates the function with the help of Pade approximants, and evaluates (extrapolates/rotates)
!   this approximation at a given set of arguments.
!   Expressions taken from K. Lee and K. Chang, Phys. Rev. B 54 R8285 (1996).
!
! !REVISION HISTORY:
!   Created December 2006 (S. Sagmeister)
!   Documentation added, July 2011 (S. Sagmeister)
!
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: m
         Complex (8), Intent (In) :: z (m)
         Integer, Intent (In) :: n
         Complex (8), Intent (In) :: iw (n), ih (n)
         Complex (8), Intent (Out) :: h (m)
    ! local variables
         Character (*), Parameter :: thisnam = 'pade'
         Complex (8), Allocatable :: acoef (:), bcoef (:), a (:), c (:, :)
         Complex (8) :: zz
         Integer :: j, l, k

    ! require order higher than two
         If (n < 2) Then
            Write (*,*) 'Error(' // thisnam // '): approximant order to&
           &o small (< 2)'
            Call terminate
         End If

    ! allocate
         Allocate (acoef(n), bcoef(0:n), a(n), c(n,n))

    ! coefficients for numerator and denominator (see reference cited above)
         a (1) = ih (1)
         c (1, :) = ih (:)
         Do j = 2, n
            Do l = 2, n
               c (j, l) = (a(j-1)-c(j-1, l)) / ((iw(l)-iw(j-1))*c(j-1, &
              & l))
            End Do
            a (j) = c (j, j)
         End Do

    ! calculate values at given frequencies
    ! loop over frequencies for which the analytic continuation is to be calculated
         Do k = 1, m
            zz = z (k)
            bcoef (:) = (1.d0, 0.d0)
            bcoef (0) = (0.d0, 0.d0)
    ! the coefficients for the continued fraction expression of the Pade approximant
    ! are set up below
    ! loop over frequencies for which the function has already been evaluated
            acoef (1) = a (1)
            Do j = 2, n
               acoef (j) = a (j) * (zz-iw(j-1))
            End Do
    ! continued fraction evaluation of the Pade approximant
            Call ctdfrac (n, acoef, bcoef, h(k))
         End Do

         Deallocate (acoef, bcoef, a, c)

      End Subroutine
!EOC

End Module
