!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writesymi
! !INTERFACE:
!
!
Subroutine writesymi
! !USES:
      Use modinput
      Use modmain
      Use modxs
! !DESCRIPTION:
!   Outputs the crystal and symmetries including their inverse
!   elements to file {\tt SYMINV.OUT}
!
! !REVISION HISTORY:
!   Created December 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Integer :: i, isym, isymi, lspl, lspli, lspn, lspni
      Real (8) :: vtlc (3), vtlic (3)
  ! output the crystal symmetries and inverses
      Open (50, File='SYMINV'//trim(filext), Action='WRITE', Form='FORMATTED')
      Write (50, '(I4, " : nsymcrys")') nsymcrys
      Do isym = 1, nsymcrys
     ! inverse crystal symmetry
         isymi = scimap (isym)
     ! spatial rotation
         lspl = lsplsymc (isym)
         lspn = lspnsymc (isym)
     ! global spin rotation
         lspli = lsplsymc (isymi)
         lspni = lspnsymc (isymi)
     ! spatial translation
         Call r3mv (input%structure%crystal%basevect, vtlsymc(1, isym), vtlc)
         Call r3mv (input%structure%crystal%basevect, vtlsymc(1, isymi), vtlic)
         Write (50,*)
         Write (50, '("Crystal symmetry, Bravais symmetry : ", 2I4)') isym, lspl
         Write (50, '(" inverse operations		     : ", 2I4)') isymi, lspli
         Write (50, '(" spatial rotation and inverse (lattice coordinates):")')
         Do i = 1, 3
            Write (50, '(3I4, 5x, 3I4)') symlat (i, :, lspl), symlat (i, :, lspli)
         End Do
         Write (50, '(" spatial rotation and inverse (Cartesian coordinates):")')
         Do i = 1, 3
            Write (50, '(3G18.10, 5x, 3G18.10)') symlatc (i, :, lspl), symlatc (i, :, lspli)
         End Do
         Write (50, '(" spatial translation and inverse (lattice coordinates) :")')
         Write (50, '(3G18.10, 5x, 3G18.10)') vtlsymc (:, isym), vtlsymc (:, isymi)
         Write (50, '(" spatial translation and inverse (Cartesian coordinates) :")')
         Write (50, '(3G18.10, 5x, 3G18.10)') vtlc, vtlic
         Write (50, '(" global spin rotation and inverse (lattice coordinates) :")')
         Do i = 1, 3
            Write (50, '(3I4, 5x, 3I4)') symlat (i, :, lspn), symlat (i, :, lspni)
         End Do
         Write (50, '(" global spin rotation and inverse (Cartesian coordinates) :")')
         Do i = 1, 3
            Write (50, '(3G18.10, 5x, 3G18.10)') symlatc (i, :, lspn), symlatc (i, :, lspni)
         End Do
      End Do
      Close (50)
End Subroutine writesymi
!EOC
