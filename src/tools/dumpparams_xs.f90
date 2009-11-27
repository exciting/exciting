!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: dumpparams_xs
! !INTERFACE:
!
!
Subroutine dumpparams_xs (string, comment)
! !USES:
      Use modinput
      Use modmain
      Use modxs
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to the parameters related to the excited states
!   implementation.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Character (*), Intent (In) :: string, comment
  ! local variables
      Integer :: i
      Open (Unit=77, File=trim(string), Action='write', Position='appen&
     &d')
      Write (77,*)
      Write (77, '("! EXCITING version ", I1.1, ".", I1.1, ".", I3.3)') &
     & version
      Write (77, '(a)') trim (comment)
      Write (77,*)
      Write (77, '("vgqlmt")')
      Write (77,*) size (input%xs%qpointset%qpoint, 2)
      Do i = 1, size (input%xs%qpointset%qpoint, 2)
         Write (77,*) input%xs%qpointset%qpoint(:, i)
      End Do
      Write (77,*)
      Write (77, '("mdfqtype")')
      Write (77,*) input%xs%tddft%mdfqtype
      Write (77,*)
      Write (77, '("gqmax")')
      Write (77,*) input%xs%gqmax
      Write (77,*)
      Write (77, '("lmaxapwwf")')
      Write (77,*) input%xs%lmaxapwwf
      Write (77,*)
      Write (77, '("fastpmat")')
      Write (77,*) input%xs%fastpmat
      Write (77,*)
      Write (77, '("fastemat")')
      Write (77,*) input%xs%fastemat
      Write (77,*)
      Write (77, '("emattype")')
      Write (77,*) input%xs%emattype
      Write (77,*)
      Write (77, '("lmaxemat")')
      Write (77,*) input%xs%lmaxemat
      Write (77,*)
      Write (77, '("torddf")')
      Write (77,*) input%xs%tddft%torddf
      Write (77,*)
      Write (77, '("tordfxc")')
      Write (77,*) input%xs%tddft%tordfxc
      Write (77,*)
      Write (77, '("acont")')
      Write (77,*) input%xs%tddft%acont
      Write (77,*)
      Write (77, '("nwacont")')
      Write (77,*) input%xs%tddft%nwacont
      Write (77,*)
      Write (77, '("broad")')
      Write (77,*) input%xs%broad
      Write (77,*)
      Write (77, '("aresdf")')
      Write (77,*) input%xs%tddft%aresdf
      Write (77,*)
      Write (77, '("epsdfde")')
      Write (77,*) input%xs%tddft%epsdfde
      Write (77,*)
      Write (77, '("emaxdf")')
      Write (77,*) input%xs%emaxdf
      Write (77,*)
      Write (77, '("dfoffdiag")')
      Write (77,*) input%xs%dfoffdiag
      Write (77,*)
      Write (77, '("tetradf")')
      Write (77,*) input%xs%tetra%tetradf
      Write (77,*)
      Write (77, '("kerndiag")')
      Write (77,*) input%xs%tddft%kerndiag
      Write (77,*)
      Write (77, '("fxctype")')
      Write (77,*) input%xs%tddft%fxctypenumber
      Write (77,*)
      Write (77, '("nexcitmax")')
      Write (77,*) input%xs%BSE%nexcitmax
      Write (77,*)
      Write (77, '("alphalrc")')
      Write (77,*) input%xs%tddft%alphalrc
      Write (77,*)
      Write (77, '("alphalrcdyn")')
      Write (77,*) input%xs%tddft%alphalrcdyn
      Write (77,*)
      Write (77, '("betalrcdyn")')
      Write (77,*) input%xs%tddft%betalrcdyn
      Write (77,*)
      Write (77, '("dftrans")')
      Write (77,*) ndftrans
      Do i = 1, ndftrans
         Write (77,*) dftrans (:, i)
      End Do
      Write (77,*)
      Write (77, '("gather")')
      Write (77,*) input%xs%gather
      Write (77,*)
      Write (77, '("symmorph")')
      Write (77,*) input%xs%symmorph
      Write (77,*)
      Write (77, '("tevout")')
      Write (77,*) input%xs%tevout
      Write (77,*)
      Write (77, '("appinfo")')
      Write (77,*) input%xs%tappinfo
      Write (77,*)
      Write (77, '("dbglev")')
      Write (77,*) input%xs%dbglev
      Write (77,*)
      Write (77, '("screentype")')
      Write (77,*) "'" // trim (input%xs%screening%screentype) // "'"
      Write (77,*)
      Write (77, '("nosymscr")')
      Write (77,*) input%xs%screening%nosym
      Write (77,*)
      Write (77, '("reducekscr")')
      Write (77,*) input%xs%screening%reducek
      Write (77,*)
      Write (77, '("ngridkscr")')
      Write (77,*) input%xs%screening%ngridk
      Write (77,*)
      Write (77, '("vkloffscr")')
      Write (77,*) input%xs%screening%vkloff
      Write (77,*)
      Write (77, '("rgkmaxscr")')
      Write (77,*) input%xs%screening%rgkmax
      Write (77,*)
      Write (77, '("nemptyscr")')
      Write (77,*) input%xs%screening%nempty
      Write (77,*)
      Write (77, '("scrherm")')
      Write (77,*) input%xs%BSE%scrherm
      Write (77,*)
      Write (77, '("bsetype")')
      Write (77,*) "'" // trim (input%xs%BSE%bsetype) // "'"
      Write (77,*)
!
      Write (77,*)
!
      Write (77,*)
!
!
      Write (77,*)
      Write (77, '("nstlce")')
      Write (77,*) nbfce, nafce
      Write (77,*)
      Write (77, '("nstlbse")')
      Write (77,*) nbfbse, nafbse
      Close (77)
End Subroutine dumpparams_xs
!EOC
