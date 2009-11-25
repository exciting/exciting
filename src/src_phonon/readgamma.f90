!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine readgamma (gq)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (Out) :: gq (3*natmtot, nqpt)
! local variables
      Integer :: iq, i
      Integer :: natmtot_, nqpt_, iq_, i_
      Real (8) :: vql_ (3), vqc_ (3)
! external functions
      Real (8) :: r3dist
      External r3dist
      Open (50, File='GAMMAQ.OUT', Action='READ', Form='FORMATTED', &
     & Status='OLD')
      Read (50,*)
      Read (50,*) natmtot_
      If (natmtot .Ne. natmtot_) Then
         Write (*,*)
         Write (*, '("Error(readgamma): differing natmtot")')
         Write (*, '(" current	  : ", I4)') natmtot
         Write (*, '(" GAMMAQ.OUT : ", I4)') natmtot_
         Write (*,*)
         Stop
      End If
      Read (50,*) nqpt_
      If (nqpt .Ne. nqpt_) Then
         Write (*,*)
         Write (*, '("Error(readgamma): differing nqpt")')
         Write (*, '(" current	  : ", I6)') nqpt
         Write (*, '(" GAMMAQ.OUT : ", I6)') nqpt_
         Write (*,*)
         Stop
      End If
      Read (50,*)
      Do iq = 1, nqpt
         Read (50,*) iq_
         If (iq .Ne. iq_) Then
            Write (*,*)
            Write (*, '("Error(readgamma): incorrect q-point index in G&
           &AMMAQ.OUT for q-point ", I6)') iq
            Write (*,*)
            Stop
         End If
         Read (50,*) vql_
         If (r3dist(vql(:, iq), vql_) .Gt. input%structure%epslat) Then
            Write (*,*)
            Write (*, '("Error(readgamma): differing q-vectors in latti&
           &ce coordinates for q-point ", I6)') iq
            Write (*, '(" current    : ", 3G18.10)') vql (:, iq)
            Write (*, '(" GAMMAQ.OUT : ", 3G18.10)') vql_
            Write (*,*)
            Stop
         End If
         Read (50,*) vqc_
         If (r3dist(vqc(:, iq), vqc_) .Gt. input%structure%epslat) Then
            Write (*,*)
            Write (*, '("Error(readgamma): differing q-vectors in Carte&
           &sian coordinates for q-point ", I6)') iq
            Write (*, '(" current    : ", 3G18.10)') vqc (:, iq)
            Write (*, '(" GAMMAQ.OUT : ", 3G18.10)') vqc_
            Write (*,*)
            Stop
         End If
         Do i = 1, 3 * natmtot
            Read (50,*) i_, gq (i, iq)
            If (i .Ne. i_) Then
               Write (*,*)
               Write (*, '("Error(readgamma): incorrect mode index in G&
              &AMMAQ.OUT for q-point ", I6)') iq
               Write (*,*)
               Stop
            End If
         End Do
         Read (50,*)
      End Do
      Close (50)
      Return
      Stop
End Subroutine
