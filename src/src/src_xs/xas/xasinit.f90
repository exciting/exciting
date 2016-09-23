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
	
!	if(allocated(mj2ml)) deallocate(mj2ml)
	allocate(mj2ml(ncg,2))                      ! Not completely general, enough for K,L1,L2,L3,M1,M2,M3,M4-edge calculations
!	if(allocated(preml)) deallocate(preml)
	allocate(preml(ncg,2))                      ! Not completely general, enough for K,L1,L2,L3,M1,M2,M3,M4-edge calculations

! Set prefactor and mj2ml-pointers for the spin spherical harmonics
	preml=0.0
	mj2ml=0
	preml(1,1)=0.0d0	! 1s(1/2,-1/2)
	preml(1,2)=1.0d0
	preml(2,1)=1.0d0    ! 1s(1/2,1/2)
	preml(2,2)=0.0d0	
	if (ncg .gt. 2) then
		preml(3,1)=1.0d0	! 2s(1/2,-1/2)
		preml(4,2)=1.0d0	! 2s(1/2,1/2)
	end if
	if (ncg .gt. 4) then
		mj2ml(5,1)=-1   ! 2p(1/2,-1/2)
		mj2ml(5,2)=0
		preml(5,1)=-sqrt(2.0d0/3.0d0)
		preml(5,2)=sqrt(1.0d0/3.0d0)
		mj2ml(6,1)=0	! 2p(1/2,1/2)
		mj2ml(6,2)=1
		preml(6,1)=-sqrt(1.0d0/3.0d0)
		preml(6,2)=sqrt(2.0d0/3.0d0)
		mj2ml(7,1)=0	! 2p(3/2,-3/2)
		mj2ml(7,2)=-1
		preml(7,1)=-0.0d0
		preml(7,2)=1.0d0
		mj2ml(8,1)=-1	! 2p(3/2,-1/2)
		mj2ml(8,2)=0
		preml(8,1)=sqrt(1.0d0/3.0d0)
		preml(8,2)=sqrt(2.0d0/3.0d0)
		mj2ml(9,1)=0	! 2p(3/2,1/2)
		mj2ml(9,2)=1
		preml(9,1)=sqrt(2.0d0/3.0d0)
		preml(9,2)=sqrt(1.0d0/3.0d0)
		mj2ml(10,1)=1	! 2p(3/2,3/2)
		mj2ml(10,2)=0
		preml(10,1)=1.0d0
		preml(10,2)=0.0d0
	end if
	if (ncg .gt. 10) then
		preml(11,1)=1.0d0	! 3s(1/2,-1/2)
		preml(12,2)=1.0d0	! 3s(1/2,1/2)
	end if
	if (ncg .gt. 12) then
		mj2ml(13,1)=-1   ! 3p(1/2,-1/2)
		mj2ml(13,2)=0
		preml(13,1)=-sqrt(2.0d0/3.0d0)
		preml(13,2)=sqrt(1.0d0/3.0d0)
		mj2ml(14,1)=0	! 3p(1/2,1/2)
		mj2ml(14,2)=1
		preml(14,1)=-sqrt(1.0d0/3.0d0)
		preml(14,2)=sqrt(2.0d0/3.0d0)
		mj2ml(15,1)=0	! 3p(3/2,-3/2)
		mj2ml(15,2)=-1
		preml(15,1)=-0.0d0
		preml(15,2)=1.0d0
		mj2ml(16,1)=-1	! 3p(3/2,-1/2)
		mj2ml(16,2)=0
		preml(16,1)=sqrt(1.0d0/3.0d0)
		preml(16,2)=sqrt(2.0d0/3.0d0)
		mj2ml(17,1)=0	! 3p(3/2,1/2)
		mj2ml(17,2)=1
		preml(17,1)=sqrt(2.0d0/3.0d0)
		preml(17,2)=sqrt(1.0d0/3.0d0)
		mj2ml(18,1)=1	! 3p(3/2,3/2)
		mj2ml(18,2)=0
		preml(18,1)=1.0d0
		preml(18,2)=0.0d0
	end if
	if (ncg .gt. 12) then
		mj2ml(19,1)=-2	! 3d(3/2,-3/2)
		mj2ml(19,2)=-1
		preml(19,1)=-sqrt(4.0d0/5.0d0)
		preml(19,2)=sqrt(1.0d0/5.0d0)
		mj2ml(20,1)=-1	! 3d(3/2,-1/2)
		mj2ml(20,2)=0
		preml(20,1)=-sqrt(3.0d0/5.0d0)
		preml(20,2)=sqrt(2.0d0/5.0d0)
		mj2ml(21,1)=0	! 3d(3/2,1/2)
		mj2ml(21,2)=1
		preml(21,1)=-sqrt(2.0d0/5.0d0)
		preml(21,2)=sqrt(3.0d0/5.0d0)
		mj2ml(22,1)=1	! 3d(3/2,3/2)
		mj2ml(22,2)=2
		preml(22,1)=-sqrt(1.0d0/5.0d0)
		preml(22,2)=sqrt(4.0d0/5.0d0)
		mj2ml(23,1)=0	! d(5/2,-5/2)
		mj2ml(23,2)=-2	
		preml(23,1)=0.0d0
		preml(23,2)=1.0d0
		mj2ml(24,1)=-2	! d(5/2,-3/2)
		mj2ml(24,2)=-1	
		preml(24,1)=sqrt(1.0d0/5.0d0)
		preml(24,2)=sqrt(4.0d0/5.0d0)
		mj2ml(25,1)=-1	! d(5/2,-1/2)
		mj2ml(25,2)=0
		preml(25,1)=sqrt(2.0d0/5.0d0)
		preml(25,2)=sqrt(3.0d0/5.0d0)
		mj2ml(26,1)=0	! d(5/2,1/2)
		mj2ml(26,2)=1
		preml(26,1)=sqrt(3.0d0/5.0d0)
		preml(26,2)=sqrt(2.0d0/5.0d0)
		mj2ml(27,1)=1 	! d(5/2,3/2)
		mj2ml(27,2)=2
		preml(27,1)=sqrt(4.0d0/5.0d0)
		preml(27,2)=sqrt(1.0d0/5.0d0)
		mj2ml(28,1)=2 	! d(5/2,5/2)
		mj2ml(28,2)=0
		preml(28,1)=1.0d0
		preml(28,2)=0.0d0
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
