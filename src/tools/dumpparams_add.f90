!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: dumpparams_add
! !INTERFACE:
!
!
Subroutine dumpparams_add (string, comment)
! !USES:
      Use modmain
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to additional parameters of Exciting.
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
      Write (77, '("optswidth")')
      Write (77,*) optswidth
      Close (77)
End Subroutine dumpparams_add
!EOC
