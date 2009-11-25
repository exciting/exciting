!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putevecfv (ik, evecfv)
      Use modmain
      Use modmpi
!
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
  ! local variables
!
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
      Integer :: recl, koffset
!
!
  ! find the record length
      Inquire (IoLength=Recl) vkl (:, ik), nmatmax, nstfv, nspnfv, &
     & evecfv
  !$OMP CRITICAL
      filetag = 'EVECFV'
      If (splittfile .Or. (rank .Eq. 0)) Then
         Open (70, File=outfilenamestring(filetag, ik), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         If (splittfile) Then
            koffset = ik - firstk (procofk(ik)) + 1
         Else
            koffset = ik
         End If
         Write (70, Rec=koffset) vkl (:, ik), nmatmax, nstfv, nspnfv, &
        & evecfv
         Close (70)
!
      End If
   !$OMP END CRITICAL
      Return
End Subroutine putevecfv
