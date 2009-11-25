!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: updatpos
! !INTERFACE:
!
!
Subroutine updatpos
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Updates the current atomic positions according to the force on each atom. If
!   ${\bf r}_{ij}^m$ is the position and ${\bf F}_{ij}^m$ is the force acting on
!   it for atom $j$ of species $i$ and after time step $m$, then the new
!   position is calculated by
!   $$ {\bf r}_{ij}^{m+1}={\bf r}_{ij}^m+\tau_{ij}^m\left({\bf F}_{ij}^m
!    +{\bf F}_{ij}^{m-1}\right), $$
!   where $\tau_{ij}^m$ is a parameter governing the size of the displacement.
!   If ${\bf F}_{ij}^m\cdot{\bf F}_{ij}^{m-1}>0$ then $\tau_{ij}^m$ is
!   increased, otherwise it is decreased.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, ispn, is, ia, ias
      Real (8) :: t1
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! compute the dot-product between the current and previous total force
            t1 = dot_product (forcetot(:, ias), forcetp(:, ias))
! if the force is in the same direction then increase step size parameter
            If (t1 .Gt. 0.d0) Then
               tauatm (ias) = tauatm (ias) + &
              & input%structureoptimization%tau0atm
            Else
               tauatm (ias) = input%structureoptimization%tau0atm
            End If
! check for negative mass
            If (spmass(is) .Gt. 0.d0) Then
               atposc (:, ia, is) = atposc (:, ia, is) + tauatm (ias) * &
              & (forcetot(:, ias)+forcetp(:, ias))
            End If
         End Do
      End Do
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! compute the lattice coordinates of the atomic positions
            Call r3mv (ainv, atposc(:, ia, is), input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
! set the previous to the current total force
            forcetp (:, ias) = forcetot (:, ias)
         End Do
      End Do
! write lattice vectors and optimised atomic positions to file
      Call writegeom (.True.)
! write the optimised interatomic distances to file
      Call writeiad (.True.)
! check for overlapping muffin-tins
      Call checkmt
! generate structure factors for G-vectors
      Call gensfacgp (ngvec, vgc, ngvec, sfacg)
! generate the characteristic function
      Call gencfun
! generate structure factors for G+k-vectors
      Do ik = 1, nkpt
         Do ispn = 1, nspnfv
            Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), &
           & ngkmax, sfacgk(:, :, ispn, ik))
         End Do
      End Do
! determine the new nuclear-nuclear energy
      Call energynn
      Return
End Subroutine
!EOC
