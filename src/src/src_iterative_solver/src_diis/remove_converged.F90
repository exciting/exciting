
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine remove_converged (evecmap, iunconverged, rnorms, n, r, h, s, &
& eigenvector, eigenvalue, trialvecs)
      Use sclcontroll, Only: diismax, epsresid, maxdiisspace
      Use modmain, Only: nstfv
      Implicit None
      Integer, Intent (In) :: n
      Integer, Intent (Inout) :: evecmap (nstfv), iunconverged
      Complex (8), Intent (Inout) :: r (n, nstfv), h (n, nstfv, &
     & diismax), s (n, nstfv, diismax)
      Complex (8), Intent (Inout) :: eigenvector (n, nstfv), trialvecs &
     & (n, nstfv, diismax)
      Real (8), Intent (Inout) :: eigenvalue (nstfv, diismax), rnorms &
     & (nstfv)
      Integer :: i, skipp, idiis, oldindex, newindex
      skipp = 0
      Do i = 1, nstfv
         If (evecmap(i) .Gt. 0) Then
            If (rnorms(evecmap(i)) .Lt. epsresid) Then
               evecmap (i) = 0
               skipp = skipp + 1
            Else
               If (skipp .Gt. 0) Then
                  oldindex = evecmap (i)
                  newindex = evecmap (i) - skipp
                  Call zcopy (n, eigenvector(1, oldindex), 1, &
                 & eigenvector(1, newindex), 1)
                  Call zcopy (n, r(1, oldindex), 1, r(1, newindex), 1)
!
                  rnorms (newindex) = rnorms (oldindex)
                  Do idiis = 1, maxdiisspace
                     eigenvalue (newindex, idiis) = eigenvalue &
                    & (oldindex, idiis)
                     Call zcopy (n, h(1, oldindex, idiis), 1, h(1, &
                    & newindex, idiis), 1)
                     Call zcopy (n, s(1, oldindex, idiis), 1, s(1, &
                    & newindex, idiis), 1)
                     Call zcopy (n, trialvecs(1, oldindex, idiis), 1, &
                    & trialvecs(1, newindex, idiis), 1)
                  End Do
                  evecmap (i) = newindex
               End If
            End If
         End If
      End Do
      iunconverged = 0
      Do i = 1, nstfv
         If (evecmap(i) .Gt. 0) Then
            iunconverged = iunconverged + 1
         End If
      End Do
 ! write(*,*)iunconverged,"map",evecmap
End Subroutine remove_converged
