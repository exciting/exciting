!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getematrad (iqr, iq)
      Use modmain
      Use modinput
      Use modxs
      Use m_genfilname
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iqr, iq
  ! local variables
      Integer :: lmax1, lmax2, lmax3, un
      Character (256) :: fname
      lmax1 = Max (input%xs%lmaxapwwf, lolmax)
      lmax2 = input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same!
      lmax3 = lmax1
      If (allocated(riaa)) deallocate (riaa)
      Allocate (riaa(0:lmax1, apwordmax, 0:lmax3, apwordmax, 0:lmax2, &
     & natmtot, ngq(iq)))
      If (allocated(riloa)) deallocate (riloa)
      Allocate (riloa(nlomax, 0:lmax3, apwordmax, 0:lmax2, natmtot, &
     & ngq(iq)))
      If (allocated(rilolo)) deallocate (rilolo)
      Allocate (rilolo(nlomax, nlomax, 0:lmax2, natmtot, ngq(iq)))
      Call genfilname (basename='EMATRAD', iq=iqr, filnam=fname)
      Call getunit (un)
      Open (un, File=trim(fname), Form='unformatted', Action='read', &
     & Status='old')
      Read (un) riaa, riloa, rilolo
      Close (un)
End Subroutine getematrad
