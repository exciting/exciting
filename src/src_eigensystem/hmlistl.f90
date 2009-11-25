!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlistl
! !INTERFACE:
!
!
Subroutine hmlistl (tapp, ngp, igpig, vgpc, v, h)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tapp  : .true. if the Hamiltonian is to be applied to the input vector,
!           .false. if the full matrix is to be calculated (in,logical)
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   v     : input vector to which H is applied if tapp is .true., otherwise
!           not referenced (in,complex(nmatmax))
!   h     : H applied to v if tapp is .true., otherwise it is the Hamiltonian
!           matrix in packed form (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the Hamiltonian matrix for the APW
!   basis functions. The Hamiltonian is given by
!   $$ H^{\rm I}({\bf G+k,G'+k})=\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}({\bf G-G'})+V^{\sigma}({\bf G-G'}), $$
!   where $V^{\sigma}$ is the effective interstitial potential for spin $\sigma$
!   and $\tilde{\Theta}$ is the characteristic function. See routine
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
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: v (nmatmax)
      Complex (8), Intent (Inout) :: h (*)
! local variables
      Integer :: i, j, k, ig, iv (3)
      Real (8) :: t1
      Complex (8) zt1
      If (tapp) Then
! apply the Hamiltonian operator to v
         Do i = 1, ngp
            Do j = i, ngp
               iv (:) = ivg (:, igpig(i)) - ivg (:, igpig(j))
               ig = ivgig (iv(1), iv(2), iv(3))
               t1 = 0.5d0 * dot_product (vgpc(:, i), vgpc(:, j))
               zt1 = veffig (ig) + t1 * cfunig (ig)
               h (i) = h (i) + zt1 * v (j)
               If (i .Ne. j) h (j) = h (j) + conjg (zt1) * v (i)
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
               t1 = 0.5d0 * dot_product (vgpc(:, i), vgpc(:, j))
               h (k) = h (k) + veffig (ig) + t1 * cfunig (ig)
            End Do
         End Do
      End If
      Return
End Subroutine
!EOC
