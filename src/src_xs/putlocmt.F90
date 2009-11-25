!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putlocmt (fname, ik, vk, vq, locmt)
      Use modmain
      Use m_getunit
      Implicit None
  ! arguments
      Character (*), Intent (In) :: fname
      Integer, Intent (In) :: ik
      Real (8), Intent (In) :: vk (3), vq (3)
      Complex (8), Intent (In) :: locmt (nstfv, nlomax,-lolmax:lolmax, &
     & natmtot)
  ! local variables
      Integer :: recl, un
      Call getunit (un)
      Inquire (IoLength=Recl) vq, vk, nstfv, nlomax, lolmax, locmt
      Open (un, File=trim(fname), Action='write', Form='unformatted', &
     & Access='direct', Recl=Recl)
      Write (un, Rec=ik) vq, vk, nstfv, nlomax, lolmax, locmt
      Close (un)
End Subroutine putlocmt
