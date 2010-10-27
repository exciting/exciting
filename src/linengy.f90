!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: linengy
! !INTERFACE:
!
!
Subroutine linengy
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the new linearisation energies for both the APW and local-orbital
!   radial functions. See the routine {\tt findband}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ja, ias, jas
      Integer :: l, ilo, io1, io2
! automatic arrays
      Logical :: done (natmmax)
      Real (8) :: vr (nrmtmax)
      logical :: tfnd
! begin loops over atoms and species
      Do is = 1, nspecies
         done (:) = .False.
         Do ia = 1, natoms (is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
!-----------------------!
!     APW functions     !
!-----------------------!
               Do l = 0, input%groundstate%lmaxapw
                  Do io1 = 1, apword (l, is)
                     If (apwve(io1, l, is)) Then
! check if previous radial functions have same default energies
                        Do io2 = 1, io1 - 1
                           If (apwve(io2, l, is)) Then
                              If (Abs(apwe0(io1, l, is)-apwe0(io2, l, &
                             & is)) .Lt. 1.d-4) Then
                                 apwe (io1, l, ias) = apwe (io2, l, &
                                & ias)
                                 Go To 10
                              End If
                           End If
                        End Do
! find the band energy starting from default
                        apwe (io1, l, ias) = apwe0 (io1, l, is)
                        Call findband (input%groundstate%findlinentype, &
                         & l, 0, input%groundstate%nprad, nrmt(is), &
                         & spr(:, is), vr, input%groundstate%deband, input%groundstate%epsband, &
                         & apwe(io1, l, ias),tfnd)
                        if (.not.tfnd) then
                          write(*,*)
                          write(*,'("Warning(linengy): linearisation energy not found")')
                          write(*,'(" for species ",I4," and atom ",I4)') is, ia
                          write(*,'(" APW angular momentum ",I4)') l
                          write(*,'(" order ",I4)') io1
                          write(*,'(" and s.c. loop ",I5)') iscl
                        end if
                     else
                       if (input%groundstate%fermilinengy) &
                         apwe(io1,l,ias)=efermi + input%groundstate%dlinengyfermi
                     End If
10                   Continue
                  End Do
               End Do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
               Do ilo = 1, nlorb (is)
                  Do io1 = 1, lorbord (ilo, is)
                     If (lorbve(io1, ilo, is)) Then
! check if previous radial functions have same default energies
                        Do io2 = 1, io1 - 1
                           If (lorbve(io2, ilo, is)) Then
                              If (Abs(lorbe0(io1, ilo, is)-lorbe0(io2, &
                             & ilo, is)) .Lt. 1.d-4) Then
                                 lorbe (io1, ilo, ias) = lorbe (io2, &
                                & ilo, ias)
                                 Go To 20
                              End If
                           End If
                        End Do
                        l = lorbl (ilo, is)
! find the band energy starting from default
                        lorbe (io1, ilo, ias) = lorbe0 (io1, ilo, is)
                        Call findband (input%groundstate%findlinentype, &
                         & l, 0, input%groundstate%nprad, nrmt(is), &
                         & spr(:, is), vr, input%groundstate%deband, input%groundstate%epsband, &
                         & lorbe(io1, ilo, ias),tfnd)
                        if (.not.tfnd) then
                          write(*,*)
                          write(*,'("Warning(linengy): linearisation energy not found")')
                          write(*,'(" for species ",I4," and atom ",I4)') is, ia
                          write(*,'(" local-orbital ",I4)') ilo
                          write(*,'(" order ",I4)') io1
                          write(*,'(" and s.c. loop",I5)') iscl
                        end if
                     else
                       if (input%groundstate%fermilinengy) &
                         lorbe(io1,ilo,ias)=efermi + input%groundstate%dlinengyfermi
                     End If
20                   Continue
                  End Do
               End Do
               done (ia) = .True.
! copy to equivalent atoms
               Do ja = 1, natoms (is)
                  If (( .Not. done(ja)) .And. (eqatoms(ia, ja, is))) &
                 & Then
                     jas = idxas (ja, is)
                     Do l = 0, input%groundstate%lmaxapw
                        Do io1 = 1, apword (l, is)
                           apwe (io1, l, jas) = apwe (io1, l, ias)
                        End Do
                     End Do
                     Do ilo = 1, nlorb (is)
                        Do io1 = 1, lorbord (ilo, is)
                           lorbe (io1, ilo, jas) = lorbe (io1, ilo, &
                          & ias)
                        End Do
                     End Do
                     done (ja) = .True.
                  End If
               End Do
            End If
! end loops over atoms and species
         End Do
      End Do
      Return
End Subroutine
!EOC
