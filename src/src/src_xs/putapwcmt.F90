!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putapwcmt (fname, ik, vk, vq, apwcmt)
      Use modmain
      Use modinput
      Use m_getunit
      Implicit None
  ! arguments
      Character (*), Intent (In) :: fname
      Integer, Intent (In) :: ik
      Real (8), Intent (In) :: vk (3), vq (3)
      Complex (8), Intent (In) :: apwcmt (nstfv, apwordmax, lmmaxapw, &
     & natmtot)
  ! local variables
      Integer :: recl, un
      Call getunit (un)
      Inquire (IoLength=Recl) vq, vk, nstfv, apwordmax, &
     & input%groundstate%lmaxapw, apwcmt
      Open (un, File=trim(fname), Action='write', Form='unformatted', &
     & Access='direct', Recl=Recl)
      Write (un, Rec=ik) vq, vk, nstfv, apwordmax, &
     & input%groundstate%lmaxapw, apwcmt
      Close (un)
End Subroutine putapwcmt
