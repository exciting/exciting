! Copyright (C) 2015 C. Vorwerk and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
!
! !ROUTINE: xasinit
!
! !INTERFACE:
subroutine xasinit

! !DESCRIPTION:
!
! This is the main initialization subroutine of the BSE-XAS program.
! 
! !USES:
  Use modmain
  Use modmpi
  Use modinput
  Use modxs
  Use modxas
  use constants, only: pi
  Implicit none
  Integer :: is, ist, m, ic, ias, l, ir, ia, k
! !REVISION HISTORY:
!
! Created June 2015 by C. Vorwerk
! Revised August 2020 by C. Vorwerk
!
!EOP
!BOC
! Initialize core energy and wavefunction arrays and some other utilities
  call coreinit(input%xs%bse%xasspecies, input%xs%bse%xasatom)
! Obtain boundaries for different edges	
  if (input%xs%bse%xasedge .eq. 'K') then
    xasstart=1
    xasstop=2
    lxas=0
  elseif (input%xs%bse%xasedge .eq. 'L1') then
    xasstart=3
    xasstop=4
    lxas=0
  elseif (input%xs%bse%xasedge .eq. 'L2') then
    xasstart=5
    xasstop=6
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'L3') then
    xasstart=7
    xasstop=10
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'L23') then
    xasstart=5
    xasstop=10
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'M1') then
    xasstart=11
    xasstop=12
    lxas=0
  elseif (input%xs%bse%xasedge .eq. 'M2') then
    xasstart=13
    xasstop=14
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'M3') then
    xasstart=15
    xasstop=18
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'M23') then
    xasstart=13
    xasstop=18
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'M4') then
    xasstart=19
    xasstop=22
    lxas=2
  elseif (input%xs%bse%xasedge .eq. 'M5') then
    xasstart=23
    xasstop=28
    lxas=2
  elseif (input%xs%bse%xasedge .eq. 'M45') then
    xasstart=19
    xasstop=28
    lxas=2
  elseif (input%xs%bse%xasedge .eq. 'N1') then
    xasstart=29
    xasstop=30
    lxas=0
  elseif (input%xs%bse%xasedge .eq. 'N2') then
    xasstart=31
    xasstop=32
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'N3') then
    xasstart=33
    xasstop=36
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'N23') then
    xasstart=31
    xasstop=36
    lxas=1
  elseif (input%xs%bse%xasedge .eq. 'N4') then
    xasstart=37
    xasstop=40
    lxas=2
  elseif (input%xs%bse%xasedge .eq. 'N5') then
    xasstart=41
    xasstop=46
    lxas=2
  elseif (input%xs%bse%xasedge .eq. 'N45') then
    xasstart=37
    xasstop=46
    lxas=2
  end if
! define number of core states in XAS calculation
  nxas=xasstop-xasstart+1
!    additional arrays used for convenience
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
!         shortcut for atomic positions
      atposl(:,ia,is)=input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
    end do
!       calculate the muffin-tin volume
    vmt(is)=4.0d0*pi*rmt(is)*rmt(is)*rmt(is)/(3.0d0*omega)
  end do
end subroutine xasinit
! EOC
