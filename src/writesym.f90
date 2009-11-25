!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writesym
! !INTERFACE:
!
!
Subroutine writesym
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Outputs the Bravais, crystal and site symmetry matrices to files
!   {\tt SYMLAT.OUT}, {\tt SYMCRYS.OUT} and {\tt SYMSITE.OUT}, respectively.
!   Also writes out equivalent atoms and related crystal symmetries to
!   {\tt EQATOMS.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ja, ias, i
      Integer :: isym, lspl, lspn
! output the Bravais lattice symmetries
      Open (50, File='SYMLAT'//trim(filext), Action='WRITE', Form='FORM&
     &ATTED')
      Write (50, '(I4, " : nsymlat")') nsymlat
      Do isym = 1, nsymlat
         Write (50,*)
         Write (50, '(I4)') isym
         Do i = 1, 3
            Write (50, '(3I4)') symlat (i, :, isym)
         End Do
      End Do
      Close (50)
! output the crystal symmetries
      Open (50, File='SYMCRYS'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Write (50,*)
      Write (50, '("(translation vectors and rotation matrices are in l&
     &attice coordinates)")')
      Write (50,*)
      Write (50, '(I4, " : nsymcrys")') nsymcrys
      Do isym = 1, nsymcrys
         Write (50,*)
         Write (50, '("Crystal symmetry : ", I4)') isym
         Write (50, '(" spatial translation :")')
         Write (50, '(3G18.10)') vtlsymc (:, isym)
         Write (50, '(" spatial rotation :")')
         lspl = lsplsymc (isym)
         Do i = 1, 3
            Write (50, '(3I4)') symlat (i, :, lspl)
         End Do
         Write (50, '(" global spin rotation :")')
         lspn = lspnsymc (isym)
         Do i = 1, 3
            Write (50, '(3I4)') symlat (i, :, lspn)
         End Do
      End Do
      Close (50)
! output the site symmetries
      Open (50, File='SYMSITE'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Write (50,*)
      Write (50, '("(rotation matrices are in lattice coordinates)")')
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Write (50,*)
            Write (50,*)
            Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') &
           & is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol), &
           & ia
            Write (50, '(I4, " : nsymsite")') nsymsite (ias)
            Do isym = 1, nsymsite (ias)
               Write (50,*)
               Write (50, '(" Site symmetry : ", I4)') isym
               Write (50, '("  spatial rotation :")')
               lspl = lsplsyms (isym, ias)
               Do i = 1, 3
                  Write (50, '(3I4)') symlat (i, :, lspl)
               End Do
               Write (50, '("  global spin rotation :")')
               lspn = lspnsyms (isym, ias)
               Do i = 1, 3
                  Write (50, '(3I4)') symlat (i, :, lspn)
               End Do
            End Do
         End Do
      End Do
      Close (50)
! output the equivalent atoms and related symmetries
      Open (50, File='EQATOMS'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Do is = 1, nspecies
         Write (50,*)
         Write (50, '("Species : ", I4, " (", A, ")")') is, trim &
        & (input%structure%speciesarray(is)%species%chemicalSymbol)
         Do ia = 1, natoms (is)
            Write (50, '(" atom ", I4, " is equivalent to atom(s)")') &
           & ia
            i = 0
            Do ja = 1, natoms (is)
               If (eqatoms(ia, ja, is)) Then
                  If ((i .Gt. 0) .And. (Mod(i, 20) .Eq. 0)) write &
                 & (50,*)
                  Write (50, '(I4)', Advance='NO') ja
                  i = i + 1
               End If
            End Do
            Write (50,*)
         End Do
      End Do
      Close (50)
      Return
End Subroutine
!EOC
