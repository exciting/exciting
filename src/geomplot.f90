!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine geomplot
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia
! Bohr to Angstroms (CODATA 2002)
      Real (8), Parameter :: au_to_ang = 0.5291772108d0
      Real (8) :: v1 (3), v2 (3), v3 (3), v4 (3), t1
      Real (8) :: dxx, dyx, dyy, dzx, dzy, dzz
! initialise universal variables
      Call init0
!------------------------------------------------!
!     write the XCrysden file to crystal.xsf     !
!------------------------------------------------!
      Open (50, File='CRYSTAL.xsf', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("CRYSTAL")')
      Write (50,*)
      Write (50, '("PRIMVEC")')
      Write (50, '(3G18.10)') input%structure%crystal%basevect(:, 1) * &
     & au_to_ang
      Write (50, '(3G18.10)') input%structure%crystal%basevect(:, 2) * &
     & au_to_ang
      Write (50, '(3G18.10)') input%structure%crystal%basevect(:, 3) * &
     & au_to_ang
      Write (50,*)
      Write (50, '("PRIMCOORD")')
      Write (50, '(2I8)') natmtot, 1
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), v1)
            Write (50, '(A, 3G18.10)') trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol), &
           & v1 (:) * au_to_ang
         End Do
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writexsf):")')
      Write (*, '(" XCrysDen file written to crystal.xsf")')
!-----------------------------------------------!
!     write the V_Sim file to crystal.ascii     !
!-----------------------------------------------!
! determine coordinate system vectors
      t1 = Sqrt (input%structure%crystal%basevect(1, &
     & 1)**2+input%structure%crystal%basevect(2, &
     & 1)**2+input%structure%crystal%basevect(3, 1)**2)
      v1 (:) = input%structure%crystal%basevect(:, 1) / t1
      t1 = Sqrt (input%structure%crystal%basevect(1, &
     & 2)**2+input%structure%crystal%basevect(2, &
     & 2)**2+input%structure%crystal%basevect(3, 2)**2)
      v2 (:) = input%structure%crystal%basevect(:, 2) / t1
      Call r3cross (v1, v2, v3)
      t1 = Sqrt (v3(1)**2+v3(2)**2+v3(3)**2)
      v3 (:) = v3 (:) / t1
      Call r3cross (v3, v1, v2)
      t1 = Sqrt (v2(1)**2+v2(2)**2+v2(3)**2)
      v2 (:) = v2 (:) / t1
      dxx = dot_product (input%structure%crystal%basevect(:, 1), v1(:))
      dyx = dot_product (input%structure%crystal%basevect(:, 2), v1(:))
      dyy = dot_product (input%structure%crystal%basevect(:, 2), v2(:))
      dzx = dot_product (input%structure%crystal%basevect(:, 3), v1(:))
      dzy = dot_product (input%structure%crystal%basevect(:, 3), v2(:))
      dzz = dot_product (input%structure%crystal%basevect(:, 3), v3(:))
      Open (50, File='crystal.ascii', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '(3G18.10)') dxx, dyx, dyy
      Write (50, '(3G18.10)') dzx, dzy, dzz
      Write (50,*)
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            v4 (1) = dot_product (atposc(:, ia, is), v1(:))
            v4 (2) = dot_product (atposc(:, ia, is), v2(:))
            v4 (3) = dot_product (atposc(:, ia, is), v3(:))
            Write (50, '(3G18.10, " ", A)') v4, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol)
         End Do
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writevsim):")')
      Write (*, '(" V_Sim file written to crystal.ascii")')
      Return
End Subroutine
