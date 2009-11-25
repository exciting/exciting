!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine symdmat (lmax, ld, dmat)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: ld
      Complex (8), Intent (Inout) :: dmat (ld, ld, nspinor, nspinor, &
     & natmtot)
! local variables
      Integer :: isym, lspl, lspn
      Integer :: ispn, jspn
      Integer :: lmmax, lm1, lm2
      Integer :: is, ia, ja, ias, jas
      Real (8) :: det, v (3), th, t1
! automatic arrays
      Logical :: done (natmmax)
! allocatable arrays
      Complex (8), Allocatable :: zflm (:, :)
      Complex (8), Allocatable :: ulm (:, :, :)
      Complex (8), Allocatable :: su2 (:, :, :)
      Complex (8), Allocatable :: dm1 (:, :, :, :, :)
      Complex (8), Allocatable :: dm2 (:, :)
      Complex (8), Allocatable :: dm3 (:, :, :, :)
      Complex (8), Allocatable :: dm4 (:, :)
      Complex (8), Allocatable :: dm5 (:, :)
      lmmax = (lmax+1) ** 2
! allocate local arrays
      Allocate (zflm(lmmax, lmmax))
      Allocate (ulm(lmmax, lmmax, nsymlat))
      Allocate (su2(nspinor, nspinor, nsymlat))
      Allocate (dm1(lmmax, lmmax, nspinor, nspinor, natmmax))
      Allocate (dm2(lmmax, lmmax))
      Allocate (dm3(lmmax, lmmax, nspinor, nspinor))
      Allocate (dm4(nspinor, nspinor), dm5(nspinor, nspinor))
! setup a complex unit matrix for (l,m) components
      zflm (:, :) = 0.d0
      Do lm1 = 1, lmmax
         zflm (lm1, lm1) = 1.d0
      End Do
      Do isym = 1, nsymlat
! construct (l,m) rotation matrix for each lattice symmetry
         Call rotzflm (symlatc(:, :, isym), lmax, lmmax, lmmax, zflm, &
        & ulm(:, :, isym))
! construct SU(2) matrix for proper rotation of spinor components
! (note that rotsu2 uses only the proper part of the rotation matrix)
         If (associated(input%groundstate%spin)) Then
            Call rotaxang (input%structure%epslat, symlatc(:, :, isym), &
           & det, v, th)
            Call axangsu2 (v, th, su2(:, :, isym))
         End If
      End Do
      t1 = 1.d0 / dble (nsymcrys)
      Do is = 1, nspecies
! make copy of the density matrices
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            dm1 (1:lmmax, 1:lmmax, :, :, ia) = dmat (1:lmmax, 1:lmmax, &
           & :, :, ias)
         End Do
         done (:) = .False.
         Do ia = 1, natoms (is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               dmat (:, :, :, :, ias) = 0.d0
               Do isym = 1, nsymcrys
                  lspl = lsplsymc (isym)
                  lspn = lspnsymc (isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
                  ja = ieqatom (ia, is, isym)
                  jas = idxas (ja, is)
! apply (l,m) symmetry matrix as U*D*conjg(U')
                  Do ispn = 1, nspinor
                     Do jspn = 1, nspinor
                        Call zgemm ('N', 'N', lmmax, lmmax, lmmax, &
                       & zone, ulm(:, :, lspl), lmmax, dm1(:, :, ispn, &
                       & jspn, ja), lmmax, zzero, dm2, lmmax)
                        Call zgemm ('N', 'C', lmmax, lmmax, lmmax, &
                       & zone, dm2, lmmax, ulm(:, :, lspl), lmmax, &
                       & zzero, dm3(:, :, ispn, jspn), lmmax)
                     End Do
                  End Do
! apply SU(2) symmetry matrix as U*D*conjg(U') and add
                  If (associated(input%groundstate%spin)) Then
                     Do lm1 = 1, lmmax
                        Do lm2 = 1, lmmax
                           dm4 (:, :) = dm3 (lm1, lm2, :, :)
                           Call z2mm (su2(:, :, lspn), dm4, dm5)
                           Call z2mmct (dm5, su2(:, :, lspn), dm4)
                           dmat (lm1, lm2, :, :, ias) = dmat (lm1, lm2, &
                          & :, :, ias) + dm4 (:, :)
                        End Do
                     End Do
                  Else
                     dmat (1:lmmax, 1:lmmax, 1, 1, ias) = dmat &
                    & (1:lmmax, 1:lmmax, 1, 1, ias) + dm3 (1:lmmax, &
                    & 1:lmmax, 1, 1)
                  End If
! end loop over crystal symmetries
               End Do
! normalise
               dmat (:, :, :, :, ias) = t1 * dmat (:, :, :, :, ias)
               done (ia) = .True.
! rotate into equivalent atoms
               Do isym = 1, nsymcrys
                  ja = ieqatom (ia, is, isym)
                  If ( .Not. done(ja)) Then
                     jas = idxas (ja, is)
                     lspl = lsplsymc (isym)
                     lspn = lspnsymc (isym)
! apply (l,m) symmetry matrix as conjg(U')*D*U (rotates atom ia into atom ja)
                     Do ispn = 1, nspinor
                        Do jspn = 1, nspinor
                           Call zgemm ('C', 'N', lmmax, lmmax, lmmax, &
                          & zone, ulm(:, :, lspl), lmmax, dmat(:, :, &
                          & ispn, jspn, ias), ld, zzero, dm2, lmmax)
                           Call zgemm ('N', 'N', lmmax, lmmax, lmmax, &
                          & zone, dm2, lmmax, ulm(:, :, lspl), lmmax, &
                          & zzero, dm3(:, :, ispn, jspn), lmmax)
                        End Do
                     End Do
! apply SU(2) symmetry matrix as conjg(U')*D*U
                     If (associated(input%groundstate%spin)) Then
                        Do lm1 = 1, lmmax
                           Do lm2 = 1, lmmax
                              dm4 (:, :) = dm3 (lm1, lm2, :, :)
                              Call z2mctm (su2(:, :, lspn), dm4, dm5)
                              Call z2mm (dm5, su2(:, :, lspn), dm4)
                              dmat (lm1, lm2, :, :, jas) = dm4 (:, :)
                           End Do
                        End Do
                     Else
                        dmat (1:lmmax, 1:lmmax, 1, 1, jas) = dm3 &
                       & (1:lmmax, 1:lmmax, 1, 1)
                     End If
                     done (ja) = .True.
                  End If
               End Do
            End If
         End Do
      End Do
      Deallocate (zflm, ulm, su2)
      Deallocate (dm1, dm2, dm3, dm4, dm5)
      Return
End Subroutine
