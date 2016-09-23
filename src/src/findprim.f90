!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: findprim
! !INTERFACE:
!
!
Subroutine findprim
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   This routine finds the smallest primitive cell which produces the same
!   crystal structure as the conventional cell. This is done by searching
!   through all the vectors which connect atomic positions and finding those
!   which leave the crystal structure invariant. Of these, the three shortest
!   which produce a non-zero unit cell volume are chosen.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, js, ia, ja, ka
      Integer :: i1, i2, i3, iv (3)
      Integer :: i, j, n, na
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Real (8) :: t1, t2
      Real (8) :: apl (3, maxatoms)
! allocatable arrays
      Real (8), Allocatable :: dp (:)
      Real (8), Allocatable :: vp (:, :)
! external functions
      Real (8) :: r3taxi
      External r3taxi
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
! make sure all atomic coordinates are in [0,1)
            Call r3frac (input%structure%epslat, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), iv)
! determine atomic Cartesian coordinates
            Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
           & atposc(:, ia, is))
         End Do
      End Do
! find the smallest set of atoms
      is = 1
      Do js = 1, nspecies
! if a species has only one atom the cell must be primitive
         If (natoms(js) .Eq. 1) Return
         If (natoms(js) .Lt. natoms(is)) is = js
      End Do
      n = 27 * natoms (is)
      Allocate (dp(n), vp(3, n))
! generate set of possible lattice vectors
      n = 0
      Do ia = 1, natoms (is)
         v1 (:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:) - &
        & input%structure%speciesarray(is)%species%atomarray(1)%atom%coord(:)
         Do i1 = - 1, 1
            v2 (1) = v1 (1) + dble (i1)
            Do i2 = - 1, 1
               v2 (2) = v1 (2) + dble (i2)
               Do i3 = - 1, 1
                  v2 (3) = v1 (3) + dble (i3)
                  t1 = Sqrt (v2(1)**2+v2(2)**2+v2(3)**2)
                  If (t1 .Lt. input%structure%epslat) Go To 20
! check if vector v2 leaves conventional cell invariant
                  Do js = 1, nspecies
                     Do ja = 1, natoms (js)
                        v3 (:) = input%structure%speciesarray(js)%species%atomarray(ja)%atom%coord(:) + v2 (:)
                        Call r3frac (input%structure%epslat, v3, iv)
                        Do ka = 1, natoms (js)
! check both positions and magnetic fields
                           t1 = r3taxi (input%structure%speciesarray(js)%species%atomarray(ka)%atom%coord(:), v3)
                           t2 = r3taxi (input%structure%speciesarray(js)%species%atomarray(ja)%atom%bfcmt(:), &
                          & input%structure%speciesarray(js)%species%atomarray(ka)%atom%bfcmt(:))
                           If ((t1 .Lt. input%structure%epslat) .And. &
                          & (t2 .Lt. input%structure%epslat)) Go To 10
                        End Do
! atom ja has no equivalent under translation by v2
                        Go To 20
10                      Continue
                     End Do
                  End Do
! cell invariant under translation by v2, so add to list
                  n = n + 1
                  Call r3mv (input%structure%crystal%basevect, v2, &
                 & vp(:, n))
                  dp (n) = Sqrt (vp(1, n)**2+vp(2, n)**2+vp(3, n)**2)
20                Continue
               End Do
            End Do
         End Do
      End Do
! find the shortest lattice vector
      j = 1
      t1 = 1.d8
      Do i = 1, n
         If (dp(i) .Lt. t1+input%structure%epslat) Then
            j = i
            t1 = dp (i)
         End If
      End Do
      input%structure%crystal%basevect(:, 1) = vp (:, j)
! find the next shortest lattice vector not parallel to the first
      j = 1
      t1 = 1.d8
      Do i = 1, n
         Call r3cross (input%structure%crystal%basevect(:, 1), vp(:, &
        & i), v1)
         t2 = Sqrt (v1(1)**2+v1(2)**2+v1(3)**2)
         If (t2 .Gt. input%structure%epslat) Then
            If (dp(i) .Lt. t1+input%structure%epslat) Then
               j = i
               t1 = dp (i)
            End If
         End If
      End Do
      input%structure%crystal%basevect(:, 2) = vp (:, j)
! find the next shortest lattice vector which gives non-zero unit cell volume
      Call r3cross (input%structure%crystal%basevect(:, 1), &
     & input%structure%crystal%basevect(:, 2), v1)
      j = 1
      t1 = 1.d8
      Do i = 1, n
         t2 = dot_product (vp(:, i), v1(:))
         If (Abs(t2) .Gt. input%structure%epslat) Then
            If (dp(i) .Lt. t1+input%structure%epslat) Then
               j = i
               t1 = dp (i)
            End If
         End If
      End Do
      input%structure%crystal%basevect(:, 3) = vp (:, j)
      Call r3minv (input%structure%crystal%basevect, ainv)
! remove redundant atoms
      Do is = 1, nspecies
         na = 0
         Do ia = 1, natoms (is)
            Call r3mv (ainv, atposc(:, ia, is), v1)
            Call r3frac (input%structure%epslat, v1, iv)
            Do ja = 1, na
               t1 = r3taxi (apl(:, ja), v1)
               If (t1 .Lt. input%structure%epslat) Go To 30
            End Do
            na = na + 1
            apl (:, na) = v1 (:)
30          Continue
         End Do
         Do ia = 1, na
            input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:) = apl (:, ia)
            Call r3mv (input%structure%crystal%basevect, apl(:, ia), &
           & atposc(:, ia, is))
         End Do
         natoms (is) = na
      End Do
      Deallocate (dp, vp)
      Return
End Subroutine
!EOC
