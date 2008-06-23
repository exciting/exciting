
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exccoulint
  use modmain
  use modmpi
  use modxs
  use ioarray
  use m_tdgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  integer, parameter :: iqmt=1
  character(256) :: fname
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,igq1,n,iflg
  integer :: ngridkt(3),iv(3),j1,j2
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp,nkkp
  logical :: nosymt,reducekt
  real(8) :: vklofft(3)
  real(8), allocatable :: potcl(:)
  complex(8), allocatable :: exclit(:,:),excli(:,:,:,:)
  complex(8), allocatable :: emat12(:,:),emat34(:,:)
  ! external functions
  integer, external :: iplocnr
  logical, external :: tqgamma

  complex(8), allocatable :: emat12k(:,:,:,:)

  call genfilname(setfilext=.true.)

  !----------------!
  !   initialize   !
  !----------------!
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  ! map variables for screened Coulomb interaction
  call initbse
  nosym=nosymscr
  ! no symmetries implemented for screened Coulomb interaction
  reducek=.false.
  ! q-point set of screening corresponds to (k,kp)-pairs
  ngridk(:)=ngridq(:)
  vkloff(:)=vkloffbse(:)
  if (nemptyscr.eq.-1) nemptyscr=nempty

  !---------------!
  !   main part   !
  !---------------!
  emattype=1
  call init0
  call init1
  call init2xs
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call tdsave0
  ! generate Gaunt coefficients
  call tdgauntgen(max(lmaxapw,lolmax),lmaxemat,max(lmaxapw,lolmax))
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(max(lmaxapwtd,lolmax),max(lmaxapwtd,lolmax),lmaxemat,tdgnt)
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

  write(*,'(a,4i6)') 'nst1,2,3,4',nst1,nst2,nst3,nst4
  write(*,'(a,4i6)') 'nst12,34,13,24',nst12,nst34,nst13,nst24

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
  allocate(exclit(nst34,nst12))
  allocate(emat12k(nst1,nst2,n,nkptnr))
  allocate(emat12(nst12,n),emat34(nst34,n))

  potcl(:)=0.d0
  excli(:,:,:,:)=zzero
  ikkp=0

  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
  call init1xs(qvkloff(1,iqmt))
  call ematqalloc
  do iknr=1,nkptnr
     write(*,*) 'generation of matrix elements: k-point:',iknr
     ! matrix elements for k and q=0
     call ematqk1(iqmt,iknr)
     emat12k(:,:,:,iknr)=xiou(:,:,:)
     deallocate(xiou,xiuo)
  end do
  emattype=1
  call ematbdcmbs(emattype)

  write(*,'(a,4i6)') 'nst1,2,3,4',nst1,nst2,nst3,nst4
  write(*,'(a,4i6)') 'nst12,34,13,24',nst12,nst34,nst13,nst24

  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  nkkp=nkptnr*(nkptnr+1)/2
  call genparidxran('p',nkkp)

  do ikkp=ppari,pparf
     call kkpmap(ikkp,nkptnr,iknr,jknr)

     iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
     iv(:)=modulo(iv(:),ngridk(:))
     ! q-point (reduced)
     iqr=iqmapr(iv(1),iv(2),iv(3))
     ! q-point (non-reduced)
     iq=iqmap(iv(1),iv(2),iv(3))

     ! set G=0 term of Coulomb potential to zero [Ambegoaker-Kohn]
     potcl(1)=0.d0
     ! set up Coulomb potential
     do igq1=2,n
        call genwiq2xs(0,iqmt,igq1,igq1,potcl(igq1))
     end do

     call genfilname(dotext='_SCR.OUT',setfilext=.true.)

     write(*,'(a,i6,2x,2i5,2x,2i5,2x,i6)') 'ikkp,iknr,jknr,iq,iqr,n',&
          ikkp,iknr,jknr,iq,iqr,n

     ! help arrays h1(cc',G) = M_G(kcc'), h2(G',vv') = conjg(M_G'(kvv'))
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
     exclit=matmul(conjg(emat12),transpose(emat34))/omega/nkptnr

!     emat12=conjg(emat12)
!     call zgemm('n','t', nst12, nst12, n, zone, emat12, &
!          nst12, emat34, nst12, zzero, exclit, nst12 )

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

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !        excli(:,:,:,:)=zzero
     !        do igq1=1,n
     !           do ist1=1,nst1
     !              do ist2=1,nst2
     !                 do ist3=1,nst1
     !                    do ist4=1,nst2
     !                       excli(ist1,ist2,ist3,ist4)= &
     !                            excli(ist1,ist2,ist3,ist4)+ &
     !                            conjg(emat12k(ist1,ist2,igq1,iknr))* &
     !			    potcl(igq1)* &
     !                            (emat12k(ist3,ist4,igq1,jknr))/omega/nkptnr
     !                    end do
     !                 end do
     !              end do
     !           end do
     !        end do
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
     call putbsemat('EXCLI.OUT',excli,ikkp,iknr,jknr,iq,iqr,nst1,nst3,nst2,nst4)
     call genfilname(dotext='_SCI.OUT',setfilext=.true.)

     ! end loop over (k,kp) pairs
  end do
  call barrier

  call findgntn0_clear
  deallocate(emat12k,exclit,emat12,emat34)
  deallocate(potcl,excli)

  !--------------!
  !   finalize   !
  !--------------!
  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Exchange Coulomb interaction&
       & finished"
end subroutine exccoulint
