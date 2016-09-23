!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: symrfmt
! !INTERFACE:
!
!
Subroutine symrfmt (lrstp, is, rot, rfmt, srfmt)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   is    : species number (in,integer)
!   rot   : rotation matrix (in,real(3,3))
!   rfmt  : input muffin-tin function (in,real(lmmaxvr,nrmtmax))
!   srfmt : output muffin-tin function (out,real(lmmaxvr,nrmtmax))
! !DESCRIPTION:
!   Applies a symmetry operation (in the form of a rotation matrix) to a real
!   muffin-tin function. The input function can also be the output function. See
!   the routines {\tt rtozflm} and {\tt rotzflm}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Integer, Intent (In) :: is
      Real (8), Intent (In) :: rot (3, 3)
      Real (8), Intent (In) :: rfmt (lmmaxvr, nrmtmax)
      Real (8), Intent (Out) :: srfmt (lmmaxvr, nrmtmax)
! local variables
      Integer :: ir, irc, nri, nro, iro
! allocatable arrays
      Complex (8), Allocatable :: zfmt1 (:, :), zfmt2 (:, :)
      Allocate (zfmt1(lmmaxvr, nrmtmax), zfmt2(lmmaxvr, nrmtmax))
! convert real function to complex spherical harmonic expansion
      nri = 0
      irc = 0
      Do ir = 1, nrmt (is), lrstp
         irc = irc + 1
         If (ir .Le. nrmtinr(is)) Then
            Call rtozflm (input%groundstate%lmaxinr, rfmt(:, ir), &
           & zfmt1(:, irc))
            srfmt (lmmaxinr+1:, ir) = 0.d0
            nri = irc
         Else
            Call rtozflm (input%groundstate%lmaxvr, rfmt(:, ir), &
           & zfmt1(:, irc))
         End If
      End Do
! first point in the outer point of the muffin-tin
      iro = nri + 1
! number of points in the outer part
      nro = irc - nri
! rotate the complex function
      Call rotzflm (rot, input%groundstate%lmaxinr, nri, lmmaxvr, &
     & zfmt1, zfmt2)
      Call rotzflm (rot, input%groundstate%lmaxvr, nro, lmmaxvr, &
     & zfmt1(:, iro), zfmt2(:, iro))
! convert complex function to real spherical harmonic expansion
      irc = 0
      Do ir = 1, nrmt (is), lrstp
         irc = irc + 1
         If (ir .Le. nrmtinr(is)) Then
            Call ztorflm (input%groundstate%lmaxinr, zfmt2(:, irc), &
           & srfmt(:, ir))
         Else
            Call ztorflm (input%groundstate%lmaxvr, zfmt2(:, irc), &
           & srfmt(:, ir))
         End If
      End Do
      Deallocate (zfmt1, zfmt2)
      Return
End Subroutine
!EOC
