!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: findsymlat
!
!
Subroutine findsymlat
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Finds the point group symmetries which leave the Bravais lattice invariant.
!   Let $A$ be the matrix consisting of the lattice vectors in columns, then
!   $$ g=A^{\rm T}A $$
!   is the metric tensor. Any $3\times 3$ matrix $S$ with elements $-1$, 0 or 1
!   is a point group symmetry of the lattice if $\det(S)$ is $-1$ or 1, and
!   $$ S^{\rm T}gS=g. $$
!   The first matrix in the set returned is the identity.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Removed arguments and simplified, April 2007 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: i, j, md, sym (3, 3)
      Integer :: i11, i12, i13, i21, i22, i23, i31, i32, i33
      Real (8) :: s (3, 3), g (3, 3), sgs (3, 3)
      Real (8) :: c (3, 3), v (3), t1
! external functions
      Integer :: i3mdet
      Real (8) :: r3taxi
      External i3mdet, r3taxi
! determine metric tensor
      Call r3mtm (input%structure%crystal%basevect, &
     & input%structure%crystal%basevect, g)
! loop over all possible symmetry matrices
      nsymlat = 0
      Do i11 = - 1, 1
         Do i12 = - 1, 1
            Do i13 = - 1, 1
               Do i21 = - 1, 1
                  Do i22 = - 1, 1
                     Do i23 = - 1, 1
                        Do i31 = - 1, 1
                           Do i32 = - 1, 1
                              Do i33 = - 1, 1
                                 sym (1, 1) = i11
                                 sym (1, 2) = i12
                                 sym (1, 3) = i13
                                 sym (2, 1) = i21
                                 sym (2, 2) = i22
                                 sym (2, 3) = i23
                                 sym (3, 1) = i31
                                 sym (3, 2) = i32
                                 sym (3, 3) = i33
! determinant of matrix
                                 md = i3mdet (sym)
! matrix should be orthogonal
                                 If (Abs(md) .Ne. 1) Go To 10
! check invariance of metric tensor
                                 s (:, :) = dble (sym(:, :))
                                 Call r3mtm (s, g, c)
                                 Call r3mm (c, s, sgs)
                                 Do i = 1, 3
                                    Do j = 1, 3
                                       If (Abs(sgs(i, j)-g(i, j)) .Gt. &
                                      & input%structure%epslat) Go To &
                                      & 10
                                    End Do
                                 End Do
! check invariance of spin-spiral q-vector if required
                                 If (isspinspiral()) Then
                                    Call r3mtv (s, &
                                   & input%groundstate%spin%vqlss, v)
                                    t1 = r3taxi &
                                   & (input%groundstate%spin%vqlss, v)
                                    If (t1 .Gt. input%structure%epslat) &
                                   & Go To 10
                                 End If
                                 nsymlat = nsymlat + 1
                                 If (nsymlat .Gt. 48) Then
                                    Write (*,*)
                                    Write (*, '("Error(findsymlat): mor&
                                   &e than 48 symmetries found")')
                                    Write (*, '(" (lattice vectors may &
                                   &be linearly dependent)")')
                                    Write (*,*)
                                    Stop
                                 End If
! store the symmetry and determinant in global arrays
                                 symlat (:, :, nsymlat) = sym (:, :)
                                 symlatd (nsymlat) = md
10                               Continue
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      If (nsymlat .Eq. 0) Then
         Write (*,*)
         Write (*, '("Error(findsymlat): no symmetries found")')
         Write (*,*)
         Stop
      End If
! make the first symmetry the identity
      Do i = 1, nsymlat
         If ((symlat(1, 1, i) .Eq. 1) .And. (symlat(1, 2, i) .Eq. 0) &
        & .And. (symlat(1, 3, i) .Eq. 0) .And. (symlat(2, 1, i) .Eq. 0) &
        & .And. (symlat(2, 2, i) .Eq. 1) .And. (symlat(2, 3, i) .Eq. 0) &
        & .And. (symlat(3, 1, i) .Eq. 0) .And. (symlat(3, 2, i) .Eq. 0) &
        & .And. (symlat(3, 3, i) .Eq. 1)) Then
            sym (:, :) = symlat (:, :, 1)
            symlat (:, :, 1) = symlat (:, :, i)
            symlat (:, :, i) = sym (:, :)
            md = symlatd (1)
            symlatd (1) = symlatd (i)
            symlatd (i) = md
            Go To 20
         End If
      End Do
20    Continue
! index to the inverse of each operation
      Do i = 1, nsymlat
         Call i3minv (symlat(:, :, i), sym)
         Do j = 1, nsymlat
            If ((symlat(1, 1, j) .Eq. sym(1, 1)) .And. (symlat(1, 2, j) &
           & .Eq. sym(1, 2)) .And. (symlat(1, 3, j) .Eq. sym(1, 3)) &
           & .And. (symlat(2, 1, j) .Eq. sym(2, 1)) .And. (symlat(2, 2, &
           & j) .Eq. sym(2, 2)) .And. (symlat(2, 3, j) .Eq. sym(2, 3)) &
           & .And. (symlat(3, 1, j) .Eq. sym(3, 1)) .And. (symlat(3, 2, &
           & j) .Eq. sym(3, 2)) .And. (symlat(3, 3, j) .Eq. sym(3, 3))) &
           & Then
               isymlat (i) = j
               Go To 30
            End If
         End Do
         Write (*,*)
         Write (*, '("Error(findsymlat): inverse operation not found")')
         Write (*, '(" for lattice symmetry ", I2)') i
         Write (*,*)
         Stop
30       Continue
      End Do
! determine the lattice symmetries in Cartesian coordinates
      Do i = 1, nsymlat
         s (:, :) = dble (symlat(:, :, i))
         Call r3mm (s, ainv, c)
         Call r3mm (input%structure%crystal%basevect, c, symlatc(:, :, &
        & i))
! warn for almost vanishing numers in Cartesian symmetry matrices
         if (any((abs(symlatc(:,:,i)).lt.input%structure%epslat).and.(abs(symlatc(:,:,i)).gt.0.d0))) then
            write(*,*)
            write(*,'("Warning(findsymlat): almost vanishing numbers in Cartesian symmetry matrix")')
            write(*,'(" for lattice symmetry: ",i6," ; symmetry matrix below")') i
            write(*,'("  ",3g18.10)') symlatc(1,:,i)
            write(*,'("  ",3g18.10)') symlatc(2,:,i)
            write(*,'("  ",3g18.10)') symlatc(3,:,i)
            write(*,*)
         end if
      End Do
      Return
End Subroutine
!EOC
