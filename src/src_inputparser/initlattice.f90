
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine initlattice
      Use modinput
      Use mod_atoms
      Use mod_lattice
      Implicit None
      Integer :: is, ia, js, ja, i
      Real (8) :: t1, v (3)
! scale the lattice vectors (scaling not referenced again in code)
      input%structure%crystal%basevect(:, 1) = &
     & input%structure%crystal%stretch(1) * &
     & input%structure%crystal%basevect(:, 1)
      input%structure%crystal%basevect(:, 2) = &
     & input%structure%crystal%stretch(2) * &
     & input%structure%crystal%basevect(:, 2)
      input%structure%crystal%basevect(:, 3) = &
     & input%structure%crystal%stretch(3) * &
     & input%structure%crystal%basevect(:, 3)
      input%structure%crystal%basevect(:, :) = &
     & input%structure%crystal%scale * &
     & input%structure%crystal%basevect(:, :)
!         Do is = 1, nspecies
!            Do ia = 1, natoms (is)
!               Do js = 1, nspecies
!                  Do ja = 1, natoms (is)
!                     Do i = 1, 3
!                        t1 = Abs (input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(i)-&
!                       & input%structure%speciesarray(js)%species%atomarray(ja)%atom%coord(i))
!                        If (t1 .Gt. input%structure%crystal%basevect(i, &
!                       & i)) input%structure%crystal%basevect(i, i) = &
!                       & t1
!                     End Do
!                  End Do
!               End Do
!            End Do
!         End Do
         Do i = 1, 3
            input%structure%crystal%basevect(i, i) = &
           & input%structure%crystal%basevect(i, i)
         End Do
! convert atomic positions from Cartesian to lattice coordinates

!GB 3.10.2012  check if system has a cartesian input
      If (input%structure%cartesian) Then
         Call r3minv (input%structure%crystal%basevect, ainv)
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               Call r3mv (ainv, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), v)
               input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:) = v (:)
            End Do
         End Do
      End If
End Subroutine
