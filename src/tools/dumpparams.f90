!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: dumpparams
! !INTERFACE:
!
!
Subroutine dumpparams (string, comment, sppath_, sc_, sc1_, sc2_, sc3_, &
& vacuum_)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to parameters of the main version of Exciting.
!
! !REVISION HISTORY:
!   Created 2007 (Sagmeister)
!   Added parameters, July 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Character (*), Intent (In) :: string, comment
      Character (*), Intent (In) :: sppath_
      Real (8), Intent (In) :: sc_, sc1_, sc2_, sc3_
      Real (8), Intent (In) :: vacuum_
  ! local variables
      Integer :: j, ia, is
      Open (Unit=77, File=trim(string), Action='write', Status='replace&
     &')
      Write (77,*)
      Write (77, '("! EXCITING version ", I1.1, ".", I1.1, ".", I3.3)') &
     & version
      Write (77, '(a)') trim (comment)
      Write (77,*)
      Write (77, '("tasks")')
      Do j = 1, ntasks
         Write (77,*) tasks (j)
      End Do
      Write (77,*)
      Write (77, '("avec")')
      Write (77,*) input%structure%crystal%basevect(1, 1), &
     & input%structure%crystal%basevect(2, 1), &
     & input%structure%crystal%basevect(3, 1)
      Write (77,*) input%structure%crystal%basevect(1, 2), &
     & input%structure%crystal%basevect(2, 2), &
     & input%structure%crystal%basevect(3, 2)
      Write (77,*) input%structure%crystal%basevect(1, 3), &
     & input%structure%crystal%basevect(2, 3), &
     & input%structure%crystal%basevect(3, 3)
      Write (77,*)
      Write (77, '("scale")')
      Write (77,*) sc_
      Write (77,*)
      Write (77, '("scale1")')
      Write (77,*) sc1_
      Write (77,*)
      Write (77, '("scale2")')
      Write (77,*) sc2_
      Write (77,*)
      Write (77, '("scale3")')
      Write (77,*) sc3_
      Write (77,*)
      Write (77, '("epslat")')
      Write (77,*) input%structure%epslat
      Write (77,*)
      Write (77, '("primcell")')
      Write (77,*) input%structure%primcell
      Write (77,*)
      Write (77, '("tshift")')
      Write (77,*) input%structure%tshift
      Write (77,*)
      Write (77, '("autokpt")')
      Write (77,*) input%groundstate%autokpt
      Write (77,*)
      Write (77, '("ngridk")')
      Write (77,*) input%groundstate%ngridk(1), &
     & input%groundstate%ngridk(2), input%groundstate%ngridk(3)
      Write (77,*)
      Write (77, '("vkloff")')
      Write (77,*) input%groundstate%vkloff(1), &
     & input%groundstate%vkloff(2), input%groundstate%vkloff(3)
      Write (77,*)
      Write (77, '("reducek")')
      Write (77,*) input%groundstate%reducek
      Write (77,*)
      Write (77, '("ngridq")')
      Write (77,*) ngridq (1), ngridq (2), ngridq (3)
      Write (77,*)
      Write (77, '("reduceq")')
      Write (77,*) input%phonons%reduceq
      Write (77,*)
      Write (77, '("rgkmax")')
      Write (77,*) input%groundstate%rgkmax
      Write (77,*)
      Write (77, '("gmaxvr")')
      Write (77,*) input%groundstate%gmaxvr
      Write (77,*)
      Write (77, '("lmaxapw")')
      Write (77,*) input%groundstate%lmaxapw
      Write (77,*)
      Write (77, '("lmaxvr")')
      Write (77,*) input%groundstate%lmaxvr
      Write (77,*)
      Write (77, '("lmaxmat")')
      Write (77,*) input%groundstate%lmaxmat
      Write (77,*)
      Write (77, '("lmaxinr")')
      Write (77,*) input%groundstate%lmaxinr
      Write (77,*)
      Write (77, '("fracinr")')
      Write (77,*) input%groundstate%fracinr
      Write (77,*)
      Write (77, '("npsden")')
      Write (77,*) input%groundstate%npsden
      Write (77,*)
      Write (77, '("spinpol")')
      Write (77,*) associated (input%groundstate%spin)
      Write (77,*)
      Write (77, '("spinorb")')
      Write (77,*) input%groundstate%spin%spinorb
      Write (77,*)
      Write (77, '("xctype")')
      Write (77,*) input%groundstate%xctypenumber
      Write (77,*)
      Write (77, '("stype")')
      Write (77,*) input%groundstate%stypenumber
      Write (77,*)
      Write (77, '("swidth")')
      Write (77,*) input%groundstate%swidth
      Write (77,*)
      Write (77, '("epsocc")')
      Write (77,*) input%groundstate%epsocc
      Write (77,*)
      Write (77, '("epschg")')
      Write (77,*) input%groundstate%epschg
      Write (77,*)
      Write (77, '("nempty")')
      Write (77,*) input%groundstate%nempty
      Write (77,*)
      Write (77, '("beta0")')
      Write (77,*) input%groundstate%beta0
      Write (77,*)
      Write (77, '("maxscl")')
      Write (77,*) input%groundstate%maxscl
      Write (77,*)
      Write (77, '("epspot")')
      Write (77,*) input%groundstate%epspot
      Write (77,*)
      Write (77, '("epsengy")')
      Write (77,*) input%groundstate%HartreeFock%epsengy
      Write (77,*)
      Write (77, '("epsforce")')
      Write (77,*) input%structureoptimization%epsforce
      Write (77,*)
      Write (77, '("cfdamp")')
      Write (77,*) input%groundstate%cfdamp
      Write (77,*)
      Write (77, '("sppath")')
      Write (77,*) "'" // trim (sppath_) // "'"
      Write (77,*)
      Write (77, '("scrpath")')
      Write (77,*) "'" // trim (scrpath) // "'"
      Write (77,*)
      Write (77, '("molecule")')
      Write (77,*) input%structure%molecule
      Write (77,*)
      Write (77, '("vacuum")')
      Write (77,*) vacuum_
      Write (77,*)
      Write (77, '("atoms")')
      Write (77,*) nspecies
      Do is = 1, nspecies
         Write (77,*) "'" // trim &
        & (input%structure%speciesarray(is)%species%speciesfile) // "'"
         Write (77,*) natoms (is)
         Do ia = 1, natoms (is)
            Write (77,*) input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(1), &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(2), &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(3), &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(1), &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(2), &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)
         End Do
      End Do
      Write (77,*)
      Write (77, '("plot1d")')
      Write (77,*) nvp1d, npp1d
      Do j = 1, nvp1d
         Write (77,*) vvlp1d (1, j), vvlp1d (2, j), vvlp1d (3, j)
      End Do
      Write (77,*)
      Write (77, '("plot2d")')
      Write (77,*) vclp2d (1, 1), vclp2d (2, 1), vclp2d (3, 1)
      Write (77,*) vclp2d (1, 2), vclp2d (2, 2), vclp2d (3, 2)
      Write (77,*) vclp2d (1, 3), vclp2d (2, 3), vclp2d (3, 3)
      Write (77,*) np2d (1), np2d (2)
      Write (77,*)
!
      Write (77, '("dos")')
      Write (77,*) input%properties%dos%nwdos, &
     & input%properties%dos%ngrdos, input%properties%dos%nsmdos
      Write (77,*) wdos (1), wdos (2)
      Write (77,*)
      Write (77, '("tau0atm")')
      Write (77,*) input%structureoptimization%tau0atm
      Write (77,*)
      Write (77, '("nstfsp")')
      Write (77,*) input%properties%fermisurfaceplot%nstfsp
      Write (77,*)
      Write (77, '("lradstp")')
      Write (77,*) input%groundstate%lradstep
      Write (77,*)
      Write (77, '("chgexs")')
      Write (77,*) input%groundstate%chgexs
      Write (77,*)
      Write (77, '("nprad")')
      Write (77,*) input%groundstate%nprad
      Write (77,*)
      Write (77, '("scissor")')
      Write (77,*) input%properties%bandstructure%scissor
      Write (77,*)
      Write (77, '("optcomp")')
      Do j = 1, noptcomp
         Write (77,*) input%properties%linresponsetensor%optcomp(:, j)
      End Do
      Write (77,*)
      Write (77, '("usegdft")')
      Write (77,*) input%xs%usegdft
      Write (77,*)
      Write (77, '("intraband")')
      Write (77,*) input%xs%tddft%intraband
      Write (77,*)
      Write (77, '("evalmin")')
      Write (77,*) input%groundstate%evalmin
      Write (77,*)
      Write (77, '("deband")')
      Write (77,*) input%groundstate%deband
      Write (77,*)
      Write (77, '("bfieldc")')
      Write (77,*) input%groundstate%spin%bfieldc
      Write (77,*)
      Write (77, '("fixspin")')
      Write (77,*) input%groundstate%spin%fixspinnumber
      Write (77,*)
      Write (77, '("momfix")')
      Write (77,*) input%groundstate%spin%momfix
      Write (77,*)
      Write (77, '("mommtfix")')
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Write (77,*) is, ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%mommtfix(:)
         End Do
      End Do
      Write (77,*)
      Write (77, '("taufsm")')
      Write (77,*) input%groundstate%spin%taufsm
      Write (77,*)
      Write (77, '("autormt")')
      Write (77,*) input%structure%autormt
      Write (77,*)
      Write (77, '("rmtapm")')
      Write (77,*) input%groundstate%rmtapm
      Write (77,*)
      Write (77, '("nosym")')
      Write (77,*) input%groundstate%nosym
      Write (77,*)
      Write (77, '("deltaph")')
      Write (77,*) input%phonons%deltaph
      Write (77,*)
      Write (77, '("phwrite")')
      Write (77,*) nphwrt
      Do j = 1, nphwrt
         Write (77,*) vqlwrt (:, j)
      End Do
      Write (77,*)
      Write (77, '("notes")')
      Do j = 1, notelns
         Write (77, '(A80)') notes (j)
      End Do
      Write (77,*)
      Write (77, '("tforce")')
      Write (77,*) input%groundstate%tforce
      Write (77,*)
      Write (77, '("tfibs")')
      Write (77,*) input%groundstate%tfibs
      Write (77,*)
      Write (77, '("maxitoep")')
      Write (77,*) input%groundstate%OEP%maxitoep
      Write (77,*)
      Write (77, '("tauoep")')
      Write (77,*) input%groundstate%OEP%tauoep
      Write (77,*)
      Write (77, '("kstlist")')
      Do j = 1, nkstlist
         Write (77,*) kstlist (:, j)
      End Do
      Write (77,*)
      Write (77, '("vklem")')
      Write (77,*) input%properties%masstensor%vklem
      Write (77,*)
      Write (77, '("deltaem")')
      Write (77,*) input%properties%masstensor%deltaem
      Write (77,*)
      Write (77, '("ndspem")')
      Write (77,*) input%properties%masstensor%ndspem
      Write (77,*)
      Write (77, '("nosource")')
      Write (77,*) input%groundstate%nosource
      Write (77,*)
      Write (77, '("spinsprl")')
      Write (77,*) input%groundstate%spin%spinsprl
      Write (77,*)
      Write (77, '("vqlss")')
      Write (77,*) input%groundstate%spin%vqlss
      Write (77,*)
      Write (77, '("nwrite")')
      Write (77,*) input%groundstate%nwrite
      Write (77,*)
      Write (77, '("tevecsv")')
      Write (77,*) input%groundstate%tevecsv
      Write (77,*)
      Write (77, '("lda+u")')
      Write (77,*) ldapu
      Do is = 1, nspecies
         Write (77,*) is, llu (is), ujlu (1, is), ujlu (2, is)
      End Do
      Write (77,*)
      Write (77, '("rdmxctype")')
      Write (77,*) input%groundstate%RDMFT%rdmxctype
      Write (77,*)
      Write (77, '("rdmmaxscl")')
      Write (77,*) input%groundstate%RDMFT%rdmmaxscl
      Write (77,*)
      Write (77, '("maxitn")')
      Write (77,*) input%groundstate%RDMFT%maxitn
      Write (77,*)
      Write (77, '("maxitc")')
      Write (77,*) input%groundstate%RDMFT%maxitc
      Write (77,*)
      Write (77, '("taurdmn")')
      Write (77,*) input%groundstate%RDMFT%taurdmn
      Write (77,*)
      Write (77, '("taurdmc")')
      Write (77,*) input%groundstate%RDMFT%taurdmc
      Write (77,*)
      Write (77, '("rdmalpha")')
      Write (77,*) input%groundstate%RDMFT%rdmalpha
      Write (77,*)
      Write (77, '("reducebf")')
      Write (77,*) input%groundstate%spin%reducebf
      Close (77)
End Subroutine dumpparams
!EOC
