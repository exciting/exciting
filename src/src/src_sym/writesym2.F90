!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writesym2
      Use modmain
      Use modsym
      Implicit None
  ! local varaibles
      Integer :: isym, jsym, igenr
      Character (32) :: str1
  ! write out multiplication table
      Write (str1,*) maxsymcrys
      Open (50, File='SYMMULT'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Write (50,*)
      If (abelsg) write (50, '("The symmetry group is Abelian (commutat&
     &ive)")')
      If (spainvsym) write (50, '("The symmetry group contains spatial inversion symmetry")')
      Write (50, '(" (first and second group element and product below)&
     &")')
      Do isym = 1, nsymcrys
         Do jsym = 1, nsymcrys
            Write (50, '(2i6, 5x, i6)') isym, jsym, sgmut (isym, jsym)
         End Do
      End Do
      Close (50)
#ifdef SYMMULTTABLE
      Open (50, File='SYMMULT_TABLE'//trim(filext), Action='WRITE', &
     & Form='FORMATTED')
      Write (50,*)
      Write (50, '(" (symmetry group multiplication table below)")')
      Do isym = 1, nsymcrys
         Write (50, '('//trim(str1)//'i3.2)') sgmut (isym, :)
      End Do
      Close (50)
#endif
  ! write out generators
      Open (50, File='SYMGENR'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Write (50,*)
      Write (50, '("Number of elements in symmetry group    : ", i6)') &
     & nsymcrys
      Write (50, '("Number of generators for symmetry group : ", i6)') &
     & ngenr
      Write (50,*)
      Do igenr = 1, ngenr
         Write (50, '("generating element:", i4, " , number of elemnts &
        &in orbit:", i4)') genr (igenr), negenr (igenr)
         Write (50, '(" (orbit of generator below)")')
         Write (50, '('//trim(str1)//'i4)') orbgenr (igenr, &
        & :negenr(igenr))
         Write (50,*)
      End Do
      Close (50)
End Subroutine writesym2
