!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhonorm
! !INTERFACE:
!
!
Subroutine rhonorm
! !USES:
      Use modmain
! !DESCRIPTION:
!   Loss of precision of the calculated total charge can result because the
!   muffin-tin density is computed on a set of $(\theta,\phi)$ points and then
!   transformed to a spherical harmonic representation. This routine adds a
!   constant to the density so that the total charge is correct. If the error in
!   total charge exceeds a certain tolerance then a warning is issued.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Changed from rescaling to adding, September 2006 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir
      Real (8) :: t1, t2
! error in average density
      t1 = (chgtot-chgcalc) / omega
! add the constant difference to the density
      t2 = t1 / y00
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               rhomt (1, ir, ias) = rhomt (1, ir, ias) + t2
            End Do
         End Do
      End Do
      rhoir (:) = rhoir (:) + t1
! add the difference to the charges
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t2 = t1 * (4.d0*pi/3.d0) * rmt (is) ** 3
            chgmt (ias) = chgmt (ias) + t2
            chgmttot = chgmttot + t2
         End Do
      End Do
      chgir = chgtot - chgmttot
      chgcalc = chgtot
      Return
End Subroutine
!EOC
