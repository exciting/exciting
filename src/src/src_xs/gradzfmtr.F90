!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!BOP
! !ROUTINE: gradzfmtr
! !INTERFACE:
!
!
Subroutine gradzfmtr (lmax, nr, r, l1, m1, ld1, ld2, fmt, gfmt)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   ld1   : leading dimension 1 (in,integer)
!   ld2   : leading dimension 2 (in,integer)
!   fmt  : real muffin-tin function (in,real(nr))
!   gfmt : gradient of zfmt (out,real(ld1,ld2,3))
! !DESCRIPTION:
!   Calculates the gradient of a muffin-tin function with real spherical
!   harmonics expansion coefficients, $f(r)$, corresponding to a specific
!   $lm$-combination. The gradient is given in a spherical harmonics
!   representation.
!   The $y$-component is divided by $i$ to be expressed as a real number.
!   See routine {\tt gradzfmt}.
!
! !REVISION HISTORY:
!   Created April 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Integer, Intent (In) :: l1, m1
      Integer, Intent (In) :: ld1
      Integer, Intent (In) :: ld2
      Real (8), Intent (In) :: fmt (nr)
      Real (8), Intent (Out) :: gfmt (ld1, ld2, 3)
  ! local variables
      Integer :: ir, l2, m2, lm2
  ! square root of two
      Real (8), Parameter :: sqtwo = 1.4142135623730950488d0
      Real (8) :: t1, t2, t3, t4, t5
      Real (8) :: tt1, tt2
  ! automatic arrays
      Real (8) :: f (nr), g1 (nr), cf (3, nr)
  ! external functions
      Real (8) :: clebgor
      External clebgor
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(gradzfmtr): lmax < 0 : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      Do ir = 1, nr
         gfmt (:, ir, :) = 0.d0
      End Do
  ! compute the radial derivatives
      f (1:nr) = fmt (1:nr)
      Call fderiv (1, nr, r, f, g1, cf)
      t1 = Sqrt (dble(l1+1)/dble(2*l1+1))
      t2 = Sqrt (dble(l1)/dble(2*l1+1))
      l2 = l1 + 1
      If (l2 .Le. lmax) Then
         lm2 = l2 ** 2
         Do m2 = - l2, l2
            lm2 = lm2 + 1
            t3 = clebgor (l2, 1, l1, m2,-1, m1)
            t4 = clebgor (l2, 1, l1, m2, 0, m1)
            t5 = clebgor (l2, 1, l1, m2, 1, m1)
            Do ir = 1, nr
               tt1 = g1 (ir)
               tt2 = t1 * (fmt(ir)*dble(l1)/r(ir)-tt1)
               gfmt (lm2, ir, 1) = gfmt (lm2, ir, 1) + ((t3-t5)/sqtwo) &
              & * tt2
               gfmt (lm2, ir, 2) = gfmt (lm2, ir, 2) + ((-t3-t5)/sqtwo) &
              & * tt2
               gfmt (lm2, ir, 3) = gfmt (lm2, ir, 3) + t4 * tt2
            End Do
         End Do
      End If
      l2 = l1 - 1
      If (l2 .Ge. 0) Then
         lm2 = l2 ** 2
         Do m2 = - l2, l2
            lm2 = lm2 + 1
            t3 = clebgor (l2, 1, l1, m2,-1, m1)
            t4 = clebgor (l2, 1, l1, m2, 0, m1)
            t5 = clebgor (l2, 1, l1, m2, 1, m1)
            Do ir = 1, nr
               tt1 = g1 (ir)
               tt2 = t2 * (fmt(ir)*dble(l1+1)/r(ir)+tt1)
               gfmt (lm2, ir, 1) = gfmt (lm2, ir, 1) + ((t3-t5)/sqtwo) &
              & * tt2
               gfmt (lm2, ir, 2) = gfmt (lm2, ir, 2) + ((-t3-t5)/sqtwo) &
              & * tt2
               gfmt (lm2, ir, 3) = gfmt (lm2, ir, 3) + t4 * tt2
            End Do
         End Do
      End If
End Subroutine gradzfmtr
!EOC
