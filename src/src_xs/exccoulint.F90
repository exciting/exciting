
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exccoulint
  use modmain
  use modmpi
  use modxs
  use ioarray
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='exccoulint'
  integer, parameter :: iqmt=1
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,igq1,n
  integer :: iv(3),j1,j2
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp,nkkp
  real(8), allocatable :: potcl(:)
  complex(8), allocatable :: exclit(:,:),excli(:,:,:,:)
  complex(8), allocatable :: emat12(:,:),emat34(:,:)
  ! external functions
  integer, external :: iplocnr
  logical, external :: tqgamma
  complex(8), allocatable :: emat12k(:,:,:,:)
  !---------------!
  !   main part   !
  !---------------!
  emattype=1
  call init0
  call init1
  call init2
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call xssave0
  ! generate Gaunt coefficients
  call xsgauntgen(max(lmaxapw,lolmax),lmaxemat,max(lmaxapw,lolmax))
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(max(lmaxapwwf,lolmax),max(lmaxapwwf,lolmax),lmaxemat,xsgnt)
  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated &
       &within lmax values:', lmaxapw,lmaxemat,lmaxapw
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of q-points: ',nqpt
  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  ! only for systems with a gap in energy
  if (.not.ksgap) then
     write(*,*)
     write(*,'("Error(",a,"): exchange Coulomb interaction works only for &
          &systems with KS-gap.")') trim(thisnam)
     write(*,*)
     call terminate
  end if
  ! check number of empty states
  if (nemptyscr.lt.nempty) then
     write(*,*)
     write(*,'("Error(",a,"): too few empty states in screening eigenvector &
          &file - the screening should include many empty states &
          &(BSE/screening)",2i8)') trim(thisnam),nempty,nemptyscr
     write(*,*)
     call terminate
  end if
  call ematbdcmbs(emattype)
  nst12=nst1*nst2
  nst34=nst3*nst4
  nst13=nst1*nst3
  nst24=nst2*nst4
  call genfilname(dotext='_SCI.OUT',setfilext=.true.)
  if (rank.eq.0) then
     call writekpts
     call writeqpts
  end if
  n=ngq(iqmt)
  call ematrad(iqmt)
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  allocate(potcl(n))
  allocate(excli(nst1,nst2,nst1,nst2))
  allocate(exclit(nst12,nst34))
  allocate(emat12k(nst1,nst2,n,nkptnr))
  allocate(emat12(nst12,n),emat34(nst34,n))
  potcl(:)=0.d0
  excli(:,:,:,:)=zzero
  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
  call init1offs(qvkloff(1,iqmt))
  call ematqalloc
  do iknr=1,nkptnr
     call chkpt(3,(/task,1,iknr/),'task,sub,k-point; matrix elements of plane &
          &wave')
     ! matrix elements for k and q=0
     call ematqk1(iqmt,iknr)
     emat12k(:,:,:,iknr)=xiou(:,:,:)
     deallocate(xiou,xiuo)
  end do
  emattype=1
  call ematbdcmbs(emattype)
  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  nkkp=(nkptnr*(nkptnr+1))/2
  call genparidxran('p',nkkp)

  do ikkp=ppari,pparf
     call chkpt(3,(/task,2,ikkp/),'task,sub,(k,kp)-pair; exchange term of &
          &BSE-Hamiltonian')
     call kkpmap(ikkp,nkptnr,iknr,jknr)
     iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
     iv(:)=modulo(iv(:),ngridk(:))
     ! q-point (reduced)
     iqr=iqmapr(iv(1),iv(2),iv(3))
     ! q-point (non-reduced)
     iq=iqmap(iv(1),iv(2),iv(3))

     ! set G=0 term of Coulomb potential to zero [Ambegaoker-Kohn]
     potcl(1)=0.d0
     ! set up Coulomb potential
     do igq1=2,n
        call genwiq2xs(0,iqmt,igq1,igq1,potcl(igq1))
     end do

     call genfilname(dotext='_SCR.OUT',setfilext=.true.)
     j1=0
     do ist2=1,nst2
        do ist1=1,nst1
           j1=j1+1
           emat12(j1,:)=emat12k(ist1,ist2,:,iknr)
        end do
     end do
     j2=0
     do ist4=1,nst2
        do ist3=1,nst1
           j2=j2+1
           emat34(j2,:)=emat12k(ist3,ist4,:,jknr)*potcl(:)
        end do
     end do

     ! * calculate exchange matrix elements
!!$     exclit=matmul(conjg(emat12),transpose(emat34))/omega/nkptnr
     emat12=conjg(emat12)
     call zgemm('n','t', nst12, nst12, n, zone/omega/nkptnr, emat12, &
          nst12, emat34, nst12, zzero, exclit, nst12 )

     ! map back to individual band indices
     j2=0
     do ist4=1,nst2
        do ist3=1,nst1
           j2=j2+1
           j1=0
           do ist2=1,nst2
              do ist1=1,nst1
                 j1=j1+1
                 excli(ist1,ist2,ist3,ist4)=exclit(j1,j2)
              end do
           end do
        end do
     end do

     if (ikkp.le.100) then
        do ist1=1,nst1
           do ist2=1,nst2
              do ist3=1,nst1
                 do ist4=1,nst2
                    write(1200,'(i5,3x,3i4,2x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
                         ist2,jknr,ist3,ist4,excli(ist1,ist2,ist3,ist4),&
                         abs(excli(ist1,ist2,ist3,ist4))
                 end do
              end do
           end do
        end do
     end if

     ! parallel write
     call putbsemat('EXCLI.OUT',excli,ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst4,nst3)
     call genfilname(dotext='_SCI.OUT',setfilext=.true.)

     ! end loop over (k,kp) pairs
  end do
  call barrier
  call findgntn0_clear
  deallocate(emat12k,exclit,emat12,emat34)
  deallocate(potcl,excli)

  write(unitout,'(a)') "Info("//trim(thisnam)//"): Exchange Coulomb interaction&
       & finished"
end subroutine exccoulint
