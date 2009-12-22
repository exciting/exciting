!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine readdveff (iq, is, ia, ip, dveffmt, dveffir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
      Complex (8), Intent (Out) :: dveffmt (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (Out) :: dveffir (ngrtot)
! local variables
      Integer :: js, iostat
      Integer :: version_ (3), nspecies_, lmmaxvr_
      Integer :: natoms_, nrcmt_, ngrid_ (3)
      Character (256) :: fext
      Call phfext (iq, is, ia, ip, fext)
      Open (50, File='DVEFF'//trim(fext), Action='READ', Form='UNFORMAT&
     &TED', Status='OLD', IoStat=IoStat)
      If (iostat .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(readdveff): error opening ", A)') 'STATE' &
        & // trim (filext)
         Write (*,*)
         Stop
      End If
      Read (50) version_
      If ((version(1) .Ne. version_(1)) .Or. (version(2) .Ne. &
     & version_(2)) .Or. (version(3) .Ne. version_(3))) Then
         Write (*,*)
         Write (*, '("Warning(readdveff): different versions")')
         Write (*, '(" current	 : ", I3.3, ".", I3.3, ".", I3.3)') &
        & version
         Write (*, '(" DVEFF.OUT : ", I3.3, ".", I3.3, ".", I3.3)') &
        & version_
      End If
      Read (50) nspecies_
      If (nspecies .Ne. nspecies_) Then
         Write (*,*)
         Write (*, '("Error(readdveff): differing nspecies")')
         Write (*, '(" current	 : ", I4)') nspecies
         Write (*, '(" DVEFF.OUT : ", I4)') nspecies_
         Write (*,*)
         Stop
      End If
      Read (50) lmmaxvr_
      If (lmmaxvr .Ne. lmmaxvr_) Then
         Write (*,*)
         Write (*, '("Error(readdveff): differing lmmaxvr")')
         Write (*, '(" current	 : ", I4)') lmmaxvr
         Write (*, '(" DVEFF.OUT : ", I4)') lmmaxvr_
         Write (*,*)
         Stop
      End If
      Do js = 1, nspecies
         Read (50) natoms_
         If (natoms(js) .Ne. natoms_) Then
            Write (*,*)
            Write (*, '("Error(readdveff): differing natoms for species&
           & ", I4)') js
            Write (*, '(" current   : ", I4)') natoms (js)
            Write (*, '(" DVEFF.OUT : ", I4)') natoms_
            Write (*,*)
            Stop
         End If
         Read (50) nrcmt_
         If (nrcmt(js) .Ne. nrcmt_) Then
            Write (*,*)
            Write (*, '("Error(readdveff): differing nrcmt for species &
           &", I4)') js
            Write (*, '(" current   : ", I6)') nrcmt (js)
            Write (*, '(" DVEFF.OUT : ", I6)') nrcmt_
            Write (*,*)
            Stop
         End If
      End Do
      Read (50) ngrid_
      If ((ngrid(1) .Ne. ngrid_(1)) .Or. (ngrid(2) .Ne. ngrid_(2)) .Or. &
     & (ngrid(3) .Ne. ngrid_(3))) Then
         Write (*,*)
         Write (*, '("Error(readdveff): differing ngrid")')
         Write (*, '(" current	 : ", 3I6)') ngrid
         Write (*, '(" DVEFF.OUT : ", 3I6)') ngrid_
         Write (*,*)
         Stop
      End If
      Read (50) dveffmt, dveffir
      Close (50)
      Return
End Subroutine
