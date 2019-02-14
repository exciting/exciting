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
  Implicit none
  Integer :: is, ist, m, ic, ias, l, ir, ia, k
! !REVISION HISTORY:
!
! Created June 2015 by C. Vorwerk
!
!EOP
!BOC
  call init0
  call init1
  call init2
  call gencore
  ncmax=0
  nclm=0
  ncg=0
  lcoremax=0
!   count all real core states and set radial wavefunction ucore=rwfcr/r
  ia=input%xs%bse%xasatom
  is=input%xs%bse%xasspecies
  ncore=0
  ic = 0
  do ist=1,spnst(is)
    if (spcore(ist,is)) then
      ncore=ncore+1
      l=spl(ist,is)
      lcoremax=max(lcoremax,l)
      do m=-spk(ist,is),spk(ist,is)-1
        ic=ic+1
      end do
    end if
  end do
        if (ic==0) then
                write(*,*)'No Core State in Species file. Check input and/or species file.'
                call terminate 
        end if
  ncmax=max(ncmax,ncore)
  nclm=max(nclm,ic)
  ncg=ic
!	Fill ucore
!	if (allocated(ucore)) deallocate(ucore)
  allocate(ucore(spnrmax,ncg))
  ias=idxas(ia,is)
  ic=0
  do ist=1,ncore
    do m=-spk(ist,is), spk(ist,is)-1
      ic=ic+1
      do ir=1,nrmt(is)
        ucore(ir,ic)=rwfcr(ir,1,ist,ias)/spr(ir,is)
      end do
    end do
  enddo ! ist
!	Set array with core state energies
  allocate(ecore(ncg))
  ic=0
  do ist=1,ncore
    do m=-spk(ist,is), spk(ist,is)-1
      ic=ic+1
      ecore(ic)=evalcr(ist, idxas(ia,is))
    end do
  end do
  
! Set j- and mj-quantum numbers
  allocate(spj(ncg))
  allocate(mj(ncg))
  ! K
  spj(1)=0.5d0  ! 1s(1/2,-1/2)
  mj(1)=-0.5d0
  spj(2)=0.5d0   ! 1s(1/2,1/2)
  mj(2)=0.5d0
  ! L1
  if (ncg .gt. 2) then
    spj(3)=0.5d0  ! 2s(1/2,-1/2)
    mj(3)=-0.5d0
    spj(4)=0.5d0  ! 2s(1/2,1/2)
    mj(4)=0.5d0
  end if
  ! L23
  if (ncg .gt. 4) then
    spj(5)=0.5d0  ! 2p(1/2,-1/2)
    mj(5)=-0.5d0
    spj(6)=0.5d0  ! 2p(1/2,1/2)
    mj(6)=0.5d0
    spj(7)=1.5d0  ! 2p(3/2,-3/2)
    mj(7)=-1.5d0 
    spj(8)=1.5d0  ! 2p(3/2,-1/2)
    mj(8)=-0.5d0
    spj(9)=1.5d0  ! 2p(3/2,1/2)
    mj(9)=0.5d0
    spj(10)=1.5d0 ! 2p(3/2,3/2)
    mj(10)=1.5d0
  end if
  ! M1
  if (ncg .gt. 10) then
    spj(11)=0.5d0 ! 3s(1/2,-1/2)
    mj(11)=-0.5d0
    spj(12)=0.5d0 ! 3s(1/2,1/2)
    mj(12)=0.5d0
  end if
  ! M23
  if (ncg .gt. 12) then
    spj(13)=0.5d0   ! 3p(1/2,-1/2)
    mj(13)=-0.5d0
    spj(14)=0.5d0   ! 3p(1/2,1/2)
    mj(14)=0.5d0
    spj(15)=1.5d0   ! 3p(3/2,-3/2)
    mj(15)=-1.5d0
    spj(16)=1.5d0   ! 3p(3/2,-1/2)
    mj(16)=-0.5d0
    spj(17)=1.5d0   ! 3p(3/2,1/2)
    mj(17)=0.5d0
    spj(18)=1.5d0   ! 3p(3/2,3/2)
    mj(18)=1.5d0
  end if
  ! M45
  if (ncg .gt. 12) then
    spj(19)=1.5d0   ! 3d(3/2,-3/2)
    mj(19)=-1.5d0
    spj(20)=1.5d0   ! 3d(3/2,-1/2)
    mj(20)=-0.5d0
    spj(21)=1.5d0   ! 3d(3/2,1/2)
    mj(21)=0.5d0
    spj(22)=1.5d0   ! 3d(3/2,3/2)
    mj(22)=1.5d0
    spj(23)=2.5d0   ! 3d(5/2,-5/2)
    mj(23)=-2.5d0
    spj(24)=2.5d0   ! 3d(5/2,-3/2)
    mj(24)=-1.5d0
    spj(25)=2.5d0   ! 3d(5/2,-1/2)
    mj(25)=-0.5d0
    spj(26)=2.5d0   ! 3d(5/2,1/2)
    mj(26)=0.5d0
    spj(27)=2.5d0   ! 3d(5/2,3/2)
    mj(27)=1.5d0
    spj(28)=2.5d0   ! 3d(5/2,5/2)
    mj(28)=2.5d0
  end if
  !N1
  if (ncg .gt. 28) then
    spj(29)=0.5d0 ! 4s(1/2,-1/2)
    mj(29)=-0.5d0
    spj(30)=0.5d0 ! 4s(1/2,1/2)
    mj(30)=0.5d0
  end if
  !N23
  if (ncg .gt. 30) then
    spj(31)=0.5d0   ! 4p(1/2,-1/2)
    mj(31)=-0.5d0
    spj(32)=0.5d0   ! 4p(1/2,1/2)
    mj(32)=0.5d0
    spj(33)=1.5d0   ! 4p(3/2,-3/2)
    mj(33)=-1.5d0
    spj(34)=1.5d0   ! 4p(3/2,-1/2)
    mj(34)=-0.5d0
    spj(35)=1.5d0   ! 4p(3/2,1/2)
    mj(35)=0.5d0
    spj(36)=1.5d0   ! 4p(3/2,3/2)
    mj(36)=1.5d0
  end if
  !N45
  if (ncg .gt. 36) then
    spj(37)=1.5d0   ! 4d(3/2,-3/2)
    mj(37)=-1.5d0
    spj(38)=1.5d0   ! 4d(3/2,-1/2)
    mj(38)=-0.5d0
    spj(39)=1.5d0   ! 4d(3/2,1/2)
    mj(39)=0.5d0
    spj(40)=1.5d0   ! 4d(3/2,3/2)
    mj(40)=1.5d0
    spj(41)=2.5d0   ! 4d(5/2,-5/2)
    mj(41)=-2.5d0
    spj(42)=2.5d0   ! 4d(5/2,-3/2)
    mj(42)=-1.5d0
    spj(43)=2.5d0   ! 4d(5/2,-1/2)
    mj(43)=-0.5d0
    spj(44)=2.5d0   ! 4d(5/2,1/2)
    mj(44)=0.5d0
    spj(45)=2.5d0   ! 4d(5/2,3/2)
    mj(45)=1.5d0
    spj(46)=2.5d0   ! 4d(5/2,5/2)
    mj(46)=2.5d0
  end if
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
