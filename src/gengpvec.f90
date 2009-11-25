!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gengpvec
! !INTERFACE:
!
!
Subroutine gengpvec (vpl, vpc, ngp, igpig, vgpl, vgpc, gpc, tpgpc)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   vpl   : p-point vector in lattice coordinates (in,real(3))
!   vpc   : p-point vector in Cartesian coordinates (in,real(3))
!   ngp   : number of G+p-vectors returned (out,integer)
!   igpig : index from G+p-vectors to G-vectors (out,integer(ngkmax))
!   vgpl  : G+p-vectors in lattice coordinates (out,real(3,ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (out,real(3,ngkmax))
!   gpc   : length of G+p-vectors (out,real(ngkmax))
!   tpgpc : (theta, phi) coordinates of G+p-vectors (out,real(2,ngkmax))
! !DESCRIPTION:
!   Generates a set of ${\bf G+p}$-vectors for the input $p$-point with length
!   less than {\tt gkmax}. These are used as the plane waves in the APW
!   functions. Also computes the spherical coordinates of each vector.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (In) :: vpc (3)
      Integer, Intent (Out) :: ngp
      Integer, Intent (Out) :: igpig (ngkmax)
      Real (8), Intent (Out) :: vgpl (3, ngkmax)
      Real (8), Intent (Out) :: vgpc (3, ngkmax)
      Real (8), Intent (Out) :: gpc (ngkmax)
      Real (8), Intent (Out) :: tpgpc (2, ngkmax)
! local variables
      Integer :: ig, igp
      Real (8) :: v (3), t1, t2
      t1 = gkmax ** 2
      igp = 0
      Do ig = 1, ngvec
         v (:) = vgc (:, ig) + vpc (:)
         t2 = v (1) ** 2 + v (2) ** 2 + v (3) ** 2
         If (t2 .Lt. t1) Then
            igp = igp + 1
            If (igp .Gt. ngkmax) Then
               Write (*,*)
               Write (*, '("Error(gengpvec): number of G+p-vectors exce&
              &eds ngkmax")')
               Write (*,*)
               Stop
            End If
! index to G-vector
            igpig (igp) = ig
! G+p-vector in lattice coordinates
            vgpl (:, igp) = dble (ivg(:, ig)) + vpl (:)
! G+p-vector in Cartesian coordinates
            vgpc (:, igp) = v (:)
! G+p-vector length and (theta, phi) coordinates
            Call sphcrd (vgpc(:, igp), gpc(igp), tpgpc(:, igp))
         End If
      End Do
      ngp = igp
      Return
End Subroutine
!EOC
