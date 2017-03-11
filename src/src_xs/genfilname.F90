!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_genfilname
      Implicit None
Contains
!
!BOP
! !ROUTINE: genfilname
! !INTERFACE:
!
!
      Subroutine genfilname (nodotpar, basename, etype, asc, bzsampl, &
     & acont, nar, tord, nlf, fxctype, fxctypestr, scrtype, bsetype, markfxcbse, &
     & tq0, oc1, oc2, iq, iqmt, procs, rank, dotext, setfilext, &
     & revertfilext, appfilext, filnam, fileext, auxtype, lambda, dirname)
! !USES:
         use modmpi, only: terminate
         Use modmain, Only: filext
         Use modxs, Only: filextrevert,skipgnd
! !DESCRIPTION:
!   Generates file name and extension according to optional input parameters (
!   see routine).
!   Interpret bzsampl variable as default (Lorentzian) for 0, as
!   tetrahedron method for 1. Trilinear method to be followed.
!
! !REVISION HISTORY:
!   Created October 2007 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Optional, Intent (In) :: bzsampl, fxctype, oc1, oc2, &
        & iq, iqmt, lambda, procs, rank
         Integer, Optional, Intent (In) :: etype
         Logical, Optional, Intent (In) :: nodotpar, asc, acont, nar, &
        & tord, nlf, tq0, markfxcbse
         Logical, Optional, Intent (In) :: revertfilext, setfilext, &
        & appfilext
         Character (*), Optional, Intent (In) :: basename, dirname, dotext, fxctypestr, &
        & scrtype, bsetype, auxtype
         Character (256), Optional, Intent (Out) :: filnam, fileext
    ! local variables
         Logical :: nodot0, revert, setfxt, appfxt, dotxt, oct, lnar
         Character (*), Parameter :: thisnam = 'genfilname'
         Character (256) :: s, s1
    ! if file extension in "modmain" is to be reset to last value: reset
    ! else store current file extension

         revert = .False.
         setfxt = .False.
         If (present(revertfilext)) revert = revertfilext
         If (present(setfilext)) setfxt = setfilext
         If (revert) Then
            filext = trim (filextrevert)
         Else If (setfxt) Then
            filextrevert = filext
         End If
         appfxt = .False.
         If (present(appfilext)) appfxt = appfilext
         dotxt = .False.
         If (present(dotext)) dotxt = .True.
         If ((appfxt .And. setfxt) .Or. (appfxt .And. dotxt)) Then
            Write (*, '(a)') 'Error(' // trim (thisnam) // '): specifie&
           &d appfxt together with setfilext or dotext'
            Call terminate
         End If
         oct = present (oc1) .And. present (oc2)
    ! dot in front of filename in parallel output for rank eq. zero
         nodot0 = .False.
         If (present(nodotpar)) nodot0 = nodotpar
    ! start with empty string
         s = ''
    ! type of band combinations for plane wave matrix elements
         If (present(etype)) Then
            Select Case (etype)
            Case (0)
          ! all band combinations
               s = trim (s) // '_FULL'
            Case (1)
          ! v-c anc c-v combinations for response function
            Case (2)
          ! v-v and c-c combinations for screened interaction
               s = trim (s) // '_SCRI'
            Case Default
               Write (*, '(a)') 'Error(' // trim (thisnam) // '): unkno&
              &wn etype: ', etype
               Call terminate
            End Select
         End If
    ! ascii output identifier
         If (present(asc)) Then
            If (asc) Then
               s = trim (s) // '_ASC'
            End If
         End If
    ! sampling of Brillouine zone
         If (present(bzsampl)) Then
            Select Case (bzsampl)
            Case (0)
          ! do nothing (Lorentzian broadening)
            Case (1)
          ! tetrahedron method
               s = trim (s) // '_TET'
            Case Default
               Write (*, '(a)') 'Error(' // trim (thisnam) // '): unkno&
              &wn bzsampl: ', bzsampl
               Call terminate
            End Select
         End If
    ! analytic continuation
         If (present(acont)) Then
            If (acont) Then
               s = trim (s) // '_AC'
            End If
         End If
    ! exclusion of anti-resonant part
         lnar = .False.
         If (present(nar)) Then
            If (nar) Then
               lnar = .True.
               s = trim (s) // '_NAR'
            End If
         End If
    ! time-ordering
         If (present(tord)) Then
            If (tord .And. ( .Not. lnar)) Then
               s = trim (s) // '_TORD'
            End If
         End If
    ! no local field effects
         If (present(nlf)) Then
            If (nlf) Then
               s = trim (s) // '_NLF'
            End If
         End If
    ! xc-kernel type (numeric code)
         If (present(fxctype)) Then
            Write (s1, '("_FXC",i2.2)') fxctype
            s = trim (s) // trim (s1)
         End If
    ! xc-kernel type (string)
         If (present(fxctypestr)) Then
            s = trim (s) // '_FXC' // trim(adjustl(fxctypestr))
         End If
    ! BSE effective Hamiltonian type
         If (present(bsetype)) Then
            Write (s1, '("_BSE",a)') trim (adjustl(bsetype))
            s = trim (s) // trim (s1)
         End If
    ! screening type in screened Coulomb interaction
         If (present(scrtype)) Then  
            Write (s1, '("_SCR",a)') trim (adjustl(scrtype))
            s = trim (s) // trim (s1)
         End If
    ! optical components
         If (present(tq0) .And. oct) Then
            If (tq0) Then
               Write (s1, '("_OC",2i1.1)') oc1, oc2
               s = trim (s) // trim (s1)
            End If
         End If
    ! mark file if it is generated in combination with fxc-BSE
         If (present(markfxcbse)) Then
            If (markfxcbse) Then
               s = trim (s) // "_FXCBSE"
            End If
         End If
    ! Q-point (finite momentum transfer)
         If (present(iqmt)) Then
            Write (s1, '("_QMT",i3.3)') iqmt
            s = trim (s) // trim (s1)
         End If
    ! q-point
         If (present(iq)) Then
            Write (s1, '("_Q",i5.5)') iq
            s = trim (s) // trim (s1)
         End If
    ! auxilliary name 
         If (present(auxtype)) Then  
            Write (s1, '("_",a)') trim (adjustl(auxtype))
            s = trim (s) // trim (s1)
         End If
    ! auxilliary index
         If (present(lambda)) Then
            Write (s1, '("_LAMBDA",i6.6)') lambda
            s = trim (s) // trim (s1)
         End If
    ! parallelization
         If (present(rank) .And. present(procs)) Then
            If ((procs > 1) .And. ((nodot0 .And. (rank > 0)) .Or. ( &
           & .Not. nodot0))) Then
          ! tag for rank
               Write (s1, '("_par",i3.3)') rank + 1
               s = trim (s) // trim (s1)
            End If
         End If
    ! extension (including the dot)
         If (present(dotext)) Then
            If (skipgnd) Then
                s = trim (s) // '.OUT'
            Else 
                s = trim (s) // trim (dotext)
            End If
         Else If (appfxt) Then
            s = trim (s) // trim (filext)
         Else
            s = trim (s) // '.OUT'
         End If
    ! assign file extension if required
         If (present(fileext)) fileext = trim (s)
    ! assign global file extension if required
         If (setfxt) filext = trim (s)
    ! basename
         If (present(basename)) s = trim (basename) // trim (s)
    ! dirname
         If (present(dirname)) s = trim(dirname)//'/'// trim (s)
    ! dot in front of filename determined by procs, rank and nodotpar
         If (present(rank) .And. present(procs)) Then
            If (((procs > 1) .And. (rank > 0)) .Or. (( .Not. nodot0) &
           & .And. (procs > 1) .And. (rank == 0))) Then
               s = '.' // trim (s)
            End If
         End If
         If (present(filnam)) filnam = trim (s)
      End Subroutine genfilname
!EOC
!
End Module m_genfilname
