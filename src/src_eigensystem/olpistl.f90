!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: olpistl
! !INTERFACE:
!
!
Subroutine olpistl (tapp, ngp, igpig, v, o)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   v     : input vector to which O is applied if tapp is .true., otherwise
!           not referenced (in,complex(nmatmax))
!   o     : O applied to v if tapp is .true., otherwise it is the overlap
!           matrix in packed form (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the overlap matrix for the APW
!   basis functions. The overlap is given by
!   $$ O^{\rm I}({\bf G+k,G'+k})=\tilde{\Theta}({\bf G-G'}), $$
!   where $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tapp
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Complex (8), Intent (In) :: v (nmatmax)
      Complex (8), Intent (Inout) :: o (*)
! local variables
      Integer :: i, j, k, iv (3), ig
      Complex (8) zt1
      If (tapp) Then
! apply the overlap operator to v
         Do i = 1, ngp
            Do j = i, ngp
               iv (:) = ivg (:, igpig(i)) - ivg (:, igpig(j))
               ig = ivgig (iv(1), iv(2), iv(3))
               zt1 = cfunig (ig)
               o (i) = o (i) + zt1 * v (j)
               If (i .Ne. j) o (j) = o (j) + conjg (zt1) * v (i)
            End Do
         End Do
      Else
! calculate the matrix elements
         k = 0
         Do j = 1, ngp
            Do i = 1, j
               k = k + 1
               iv (:) = ivg (:, igpig(i)) - ivg (:, igpig(j))
               ig = ivgig (iv(1), iv(2), iv(3))
               o (k) = o (k) + cfunig (ig)
            End Do
         End Do
      End If
      Return
End Subroutine
!EOC
