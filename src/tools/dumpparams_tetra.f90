!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: dumpparams_tetra
! !INTERFACE:
!
!
Subroutine dumpparams_tetra (string, comment)
! !USES:
      Use modinput
      Use modmain
      Use modtetra
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to the tetrahedron method parameters of Exciting.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Character (*), Intent (In) :: string, comment
      Open (Unit=77, File=trim(string), Action='write', Position='appen&
     &d')
      Write (77,*)
      Write (77, '("! EXCITING version ", I1.1, ".", I1.1, ".", I3.3)') &
     & version
      Write (77, '(a)') trim (comment)
      Write (77,*)
      Write (77, '("tetraocc")')
      Write (77,*) input%xs%tetra%tetraocc
      Write (77,*)
      Write (77, '("tetraopt")')
      Write (77,*) tetraopt
      Write (77,*)
      Write (77, '("tetrakordexc")')
      Write (77,*) input%xs%tetra%kordexc
      Close (77)
End Subroutine dumpparams_tetra
!EOC
