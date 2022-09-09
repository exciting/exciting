! Copyright (C) 2015 C. Vorwerk and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
!
! !ROUTINE: coreinit
!
! !INTERFACE:
subroutine coreinit (is, ia)

! !INPUT/OUTPUT PARAMETERS
! is : Species index for which the core states are initialized
! ia : Atom index for which the core states are initialized
!
! !DESCRIPTION:
!
! This is the routine to initialize arrays to store core eigenvalues,
! wavefunctions and other utilities. 
! 
! !USES:
  Use modmpi, only: terminate
  use mod_atoms, only: spnst, spcore, spl, spk, spnrmax, idxas, spr
  use mod_muffin_tin, only: nrmt
  use mod_corestate, only: rwfcr, evalcr
  Use modinput
  Use modxs
  Use modxas
  
  Implicit none
  !Arguments
  Integer, Intent(In) :: is ! Species Index
  Integer, Intent(In) :: ia ! Atom Index
  !Local Variables
  Integer :: ist, write_bandstr_hdf5m, ic, ias, l, m, ir, k
! !REVISION HISTORY:
!
! Separated from xasinit.f90 August 2020 by C. Vorwerk
!
!EOP
!BOC
  call gencore
  ncmax=0
  nclm=0
  ncg=0
  lcoremax=0
!   count all real core states and set radial wavefunction ucore=rwfcr/r
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
  
! Set j- and mj-quantum numbers for core levels 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d
  ! spj (   j quantum number of core level ) : J=L+-S
  ! mj  ( m_j quantum number of core level ) : M_J={-J,-J+1,...,J-1,J}
  allocate( spj( ncg ) )
  allocate(  mj( ncg ) )
  !
  ! Below, j and m_j quantum numbers are defined by core level, e.g. 1s, 2s, 2p ...
  ! Number pairs within the comment, e.g. (3/2, -1/2), indicate (j, m_j) quantum numbers
  !
  ! K (1s core level)
  spj(1)=0.5d0      ! 1s(1/2,-1/2)
  mj(1)=-0.5d0
  spj(2)=0.5d0      ! 1s(1/2,1/2)
  mj(2)=0.5d0
  !
  ! L1 (2s core level)
  if (ncg .gt. 2) then
    spj(3)=0.5d0    ! 2s(1/2,-1/2)
    mj(3)=-0.5d0
    spj(4)=0.5d0    ! 2s(1/2,1/2)
    mj(4)=0.5d0
  end if
  !
  ! L23 (2p core level)
  if (ncg .gt. 4) then
    spj(5)=0.5d0    ! 2p(1/2,-1/2)
    mj(5)=-0.5d0
    spj(6)=0.5d0    ! 2p(1/2,1/2)
    mj(6)=0.5d0
    spj(7)=1.5d0    ! 2p(3/2,-3/2)
    mj(7)=-1.5d0 
    spj(8)=1.5d0    ! 2p(3/2,-1/2)
    mj(8)=-0.5d0
    spj(9)=1.5d0    ! 2p(3/2,1/2)
    mj(9)=0.5d0
    spj(10)=1.5d0   ! 2p(3/2,3/2)
    mj(10)=1.5d0
  end if
  !
  ! M1 (3s core level)
  if (ncg .gt. 10) then
    spj(11)=0.5d0   ! 3s(1/2,-1/2)
    mj(11)=-0.5d0
    spj(12)=0.5d0   ! 3s(1/2,1/2)
    mj(12)=0.5d0
  end if
  !
  ! M23 (3p core level)
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
  !
  ! M45 (3d core level)
  if (ncg .gt. 18) then
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
  !
  ! N1 (4s core level)
  if (ncg .gt. 28) then
    spj(29)=0.5d0   ! 4s(1/2,-1/2)
    mj(29)=-0.5d0
    spj(30)=0.5d0   ! 4s(1/2,1/2)
    mj(30)=0.5d0
  end if
  !
  ! N23 (4p core level)
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
  !
  ! N45 (4d core level)
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

end subroutine coreinit
! EOC
