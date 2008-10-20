
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint
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
  character(*), parameter :: thisnam='scrcoulint'
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iqrnr,jsym,jsymi,igq1,igq2,n,iflg,recl
  integer :: ngridkt(3),iv(3),ivgsym(3),j1,j2,nkkp
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp
  logical :: nosymt,reducekt,tq0,nsc,tphf
  real(8) :: vklofft(3),vqr(3),vq(3),t1
  real(8), allocatable :: potcl(:,:)
  integer :: igqmap(maxsymcrys),sc(maxsymcrys),ivgsc(3,maxsymcrys)
  complex(8) :: zt1
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:)
  complex(8), allocatable :: scrni(:,:,:),tm(:,:),tmi(:,:),bsedt(:,:)
  complex(8), allocatable :: phf(:,:),emat12(:,:),emat34(:,:)
  ! external functions
  integer, external :: iplocnr,idxkkp
  logical, external :: tqgamma
!!$  real(8) :: cpu0,cpu1,cpu2,cpu3
!!$  real(8) :: cpu_init1offs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1,cpu_ematqdealloc
!!$  real(8) :: cpu_clph,cpu_suma,cpu_write

!!$  !----------------!
!!$  !   initialize   !
!!$  !----------------!
!!$  ! save global variables
!!$  nosymt=nosym
!!$  reducekt=reducek
!!$  ngridkt(:)=ngridk(:)
!!$  vklofft(:)=vkloff(:)
!!$  nosym=nosymscr
!!$  ! no symmetries implemented for screened Coulomb interaction
!!$  reducek=.false.
!!$  vkloff(:)=vkloffbse(:)
!!$  if (nemptyscr.eq.-1) nemptyscr=nempty

  !---------------!
  !   main part   !
  !---------------!
  ! q-point set of screening corresponds to (k,kp)-pairs
  ngridk(:)=ngridq(:)


  emattype=2
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
     write(*,'("Error(",a,"): screened Coulomb interaction works only for &
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

  ! local arrays
  allocate(phf(ngqmax,ngqmax),potcl(ngqmax,ngqmax))
  allocate(sccli(nst1,nst3,nst2,nst4))
  allocate(scrni(ngqmax,ngqmax,nqptr))
  potcl(:,:)=0.d0
  sccli(:,:,:,:)=zzero
  scrni(:,:,:)=zzero

  ! set file extension
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)

  !-----------------------------------!
  !     loop over reduced q-points    !
  !-----------------------------------!
  call genparidxran('q',nqptr)
  do iqr=qpari,qparf
     ! locate reduced q-point in non-reduced set
     iqrnr=iplocnr(ivqr(1,iqr),ngridq)
     n=ngq(iqrnr)
     ! obtain inverse of dielectric matrix
     call geniscreen(iqr,ngqmax,n,scrni(1,1,iqr))
  end do
  ! communicate array-parts wrt. q-points
  call zalltoallv(scrni,ngqmax**2,nqptr)

  !---------------------------------------!
  !     loop over non-reduced q-points    !
  !---------------------------------------!
  call genparidxran('q',nqpt)
  do iq=qpari,qparf
     write(*,*) 'radial integrals for q-point:',iq
     call putematrad(iq)
  end do
  call barrier

  ! information on size of output file
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4,sccli
  write(*,'(a,f12.3)') 'file size for screened Coulomb &
       &interaction (GB):',recl*nkptnr*(nkptnr+1)/2.d0/1024.d0**3

  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  nkkp=nkptnr*(nkptnr+1)/2
  call genparidxran('p',nkkp)
  allocate(bsedt(3,nkkp))
  bsedt(1,:)=1.d8
  bsedt(2,:)=-1.d8
  bsedt(3,:)=zzero

  do ikkp=ppari,pparf
     call kkpmap(ikkp,nkptnr,iknr,jknr)
     
     ! k-point difference
     iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
     iv(:)=modulo(iv(:),ngridk(:))
     ! q-point (reduced)
     iqr=iqmapr(iv(1),iv(2),iv(3))
     vqr(:)=vqlr(:,iqr)
     ! q-point (non-reduced)
     iq=iqmap(iv(1),iv(2),iv(3))
     vq(:)=vql(:,iq)
     ! locate reduced q-point in non-reduced set
     iqrnr=iplocnr(ivqr(1,iqr),ngridq)
     ! local field effects size
     tq0=tqgamma(iq)
     n=ngq(iq)
     
     write(*,'(i5,3x,2i5,2x,i5,2x,2i5)') ikkp,iknr,jknr,iq,iqr,iqrnr
     
     allocate(emat12(nst12,n),emat34(nst34,n))
     allocate(tm(n,n),tmi(n,n))
     allocate(scclit(nst12,nst34)) !!!****SAG
     
     call findsymeqiv(vq,vqr,nsc,sc,ivgsc)
     call findgqmap(iq,iqr,nsc,sc,ivgsc,ngqmax,n,jsym,jsymi,ivgsym,igqmap)
     
     ! generate phase factor for dielectric matrix due to non-primitive
     ! translations
     call genphasedm(iq,jsym,ngqmax,n,phf,tphf)
     
     ! rotate inverse of screening
     do igq1=1,n
        j1=igqmap(igq1)
        do igq2=1,n
           j2=igqmap(igq2)
           if (tphf) then
              tmi(igq1,igq2)=scrni(j2,j1,iqr)
           else
              tmi(igq1,igq2)=scrni(j1,j2,iqr)
           end if
        end do
     end do
     
     ! set up Coulomb potential
     do igq1=1,n
        do igq2=igq1,n
           ! calculate weights for Coulomb potential
           iflg=0
           if (tq0.and.((igq1.eq.1).or.(igq2.eq.1))) then
              ! consider only 1/q and 1/q^2 cases for q goint to zero
              iflg=bsediagweight
           end if
           call genwiq2xs(iflg,iq,igq1,igq2,potcl(igq1,igq2))
           potcl(igq2,igq1)=potcl(igq1,igq2)
        end do
     end do
     
     ! calculate matrix elements of the plane wave
     emattype=2
     call ematbdcmbs(emattype)
     call getematrad(iq)
     call ematqalloc
     call ematqk1(iq,iknr)
     emattype=2
     call ematbdcmbs(emattype)

     if (ikkp.eq.1) write(*,'(a,2g18.10)') 'W-diagonal long range term:',scrni(1,1,1)*potcl(1,1)/nkptnr/omega

     ! select screening level
     tm(:,:)=zzero
     select case (trim(screentype))
     case('longrange')
        ! constant screening (q=0 average tensor)
        forall(igq1=1:n)
           tm(igq1,igq1)=scrni(1,1,1)*potcl(igq1,igq1)
        end forall
     case('diag')
        ! only diagonal of screening
        forall(igq1=1:n)
           tm(igq1,igq1)=tmi(igq1,igq1)*potcl(igq1,igq1)
        end forall
     case('full')
        ! full screening
        tm(:,:)=tmi(:,:)*potcl(1:n,1:n)
     end select
     
     ! combine indices for matrix elements of plane wave
     j1=0
     do ist2=1,nst2
        do ist1=1,nst1
           j1=j1+1
           emat12(j1,:)=xiou(ist1,ist2,:)
        end do
     end do
     j2=0
     do ist4=1,nst4
        do ist3=1,nst3
           j2=j2+1
           emat34(j2,:)=xiuo(ist3,ist4,:)
        end do
     end do
     
     ! matrix elements of direct term (as in BSE-code of Peter and
     ! in the SELF-documentation of Andrea Marini)
!!$     scclit=matmul(emat34,matmul(transpose(tm),transpose(conjg(emat12)))) &
!!$          /omega/nkptnr

     scclit=matmul(conjg(emat12),matmul(tm,transpose(emat34)))/omega/nkptnr

	  
!!$     ! map back to individual band indices
!!$     j2=0
!!$     do ist4=1,nst4
!!$        do ist3=1,nst3
!!$           j2=j2+1
!!$           j1=0
!!$           do ist2=1,nst2
!!$              do ist1=1,nst1
!!$                 j1=j1+1
!!$                 sccli(ist1,ist3,ist2,ist4)=scclit(j2,j1)
!!$              end do
!!$           end do
!!$        end do
!!$     end do
     
     ! map back to individual band indices
     j2=0
     do ist4=1,nst4
        do ist3=1,nst3
           j2=j2+1
           j1=0
           do ist2=1,nst2
              do ist1=1,nst1
                 j1=j1+1
                 sccli(ist1,ist3,ist2,ist4)=scclit(j1,j2)
              end do
           end do
        end do
     end do
     
!<backup comment="slow version for checking">
!!$        ! set up the screened Coulomb interaction (diagonal screening)
!!$        sccli(:,:,:,:)=zzero
!!$        do igq1=1,n
!!$           do ist1=1,nst1
!!$              do ist3=1,nst3
!!$                 do ist2=1,nst2
!!$                    do ist4=1,nst4
!!$                       sccli(ist1,ist3,ist2,ist4)= &
!!$                            sccli(ist1,ist3,ist2,ist4)+ &
!!$                            tm(igq1,igq1)*  &
!!$                            conjg(xiou(ist1,ist2,igq1))* &
!!$                            xiuo(ist3,ist4,igq1) /omega/nkptnr
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
        
!!$        ! set up the screened Coulomb interaction (full screening)
!!$        sccli(:,:,:,:)=zzero
!!$        do igq1=1,n
!!$           do igq2=1,n
!!$              do ist1=1,nst1
!!$                 do ist3=1,nst3
!!$                    do ist2=1,nst2
!!$                       do ist4=1,nst4
!!$                          sccli(ist1,ist3,ist2,ist4)= &
!!$                               sccli(ist1,ist3,ist2,ist4)+ &
!!$                               tm(igq1,igq2)*  &
!!$                               conjg(xiou(ist1,ist2,igq1))* &
!!$                               xiuo(ist3,ist4,igq2) /omega/nkptnr
!!$                       end do
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!</backup>

     if (ikkp.le.100) then
        ! write to ASCII file
        do ist1=1,nst1
           do ist3=1,nst3
              do ist2=1,nst2
                 do ist4=1,nst4
		    zt1=sccli(ist1,ist3,ist2,ist4)
                    write(1100,'(i5,3x,3i4,2x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
                         ist3,jknr,ist2,ist4,zt1,abs(zt1)**2, &
			 atan2(aimag(zt1),dble(zt1))/pi
                 end do
              end do
           end do
        end do
     end if
     
!!$     ! processe write sequentially to one file
!!$     do j=0,procs-1
!!$        if (rank.eq.j) then
!!$           call putbsemat('SCCLI.OUT',sccli,ikkp,iknr,jknr,iq,iqr, &
!!$                nst1,nst3,nst4,nst2)
!!$        end if
!!$        call barrier
!!$     end do

     if (iknr.eq.jknr) then
        do ist1=1,nst1
           do ist3=1,nst3
              zt1=sccli(ist1,ist3,ist1,ist3)
              write(1107,'(i5,3x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
                   ist3,zt1,abs(zt1)**2, &
                   atan2(aimag(zt1),dble(zt1))/pi
              t1=dble(zt1)
              bsedt(1,ikkp)=min(dble(bsedt(1,ikkp)),t1)
              bsedt(2,ikkp)=max(dble(bsedt(2,ikkp)),t1)
              bsedt(3,ikkp)=bsedt(3,ikkp)+zt1/(nst1*nst3)
           end do
        end do
        write(*,*) 'kkp:bsedl,bsedu,bsed',bsedt(1,ikkp),bsedt(2,ikkp), &
             bsedt(3,ikkp)
     end if

     ! parallel write
     call putbsemat('SCCLI.OUT',sccli,ikkp,iknr,jknr,iq,iqr,nst1,nst3,nst2,nst4)

     deallocate(emat12,emat34)
     deallocate(tm,tmi)
     deallocate(scclit)

     ! end loop over (k,kp)-pairs
  end do

  call barrier
  
  ! communicate array-parts wrt. q-points
  call zalltoallv(bsedt,3,nkkp)
  ! BSE kernel diagonal parameters
  bsedl=minval(dble(bsedt(1,:)))
  bsedu=maxval(dble(bsedt(2,:)))
  bsedd=bsedu-bsedl
  bsed=sum(bsedt(3,:))/nkptnr
  deallocate(bsedt)
  ! write BSE kernel diagonal parameters
  call putbsediag

  call findgntn0_clear

!!$  !--------------!
!!$  !   finalize   !
!!$  !--------------!
!!$  ! restore global variables
!!$  nosym=nosymt
!!$  reducek=reducekt
!!$  ngridk(:)=ngridkt(:)
!!$  vkloff(:)=vklofft(:)
!!$

  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screened Coulomb interaction&
       & finished"

end subroutine scrcoulint
