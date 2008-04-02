
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint
  use modmain
  use modmpi
  use modxs
  use ioarray
  use m_tdgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  character(256) :: fname
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iqrnr,isym,jsym,jsymi,igq1,igq2,n,iflg,recl
  integer :: ngridkt(3),iv(3),ivgsym(3),un,j1,j2
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp
  logical :: nosymt,reducekt,tq0,nsc,tphf
  complex(8) :: zt1
  real(8) :: vklofft(3),vqr(3),vq(3),v2(3),s(3,3),si(3,3),t3
  real(8), allocatable :: potcl(:,:)
  integer :: igqmap(maxsymcrys),sc(maxsymcrys),ivgsc(3,maxsymcrys)
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:,:)
  complex(8), allocatable :: scrni(:,:,:),tm(:,:),tmi(:,:)
  complex(8), allocatable :: phf(:,:),emat12(:,:),emat34(:,:)
  logical, allocatable :: done(:)
  ! external functions
  integer, external :: iplocnr
  logical, external :: tqgamma

  real(8) :: cpu0,cpu1,cpu2,cpu3
  real(8) :: cpu_init1xs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1,cpu_ematqdealloc
  real(8) :: cpu_clph,cpu_suma,cpu_write
  complex(8), allocatable :: emat12k(:,:,:),emat12kp(:,:,:)

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
  emattype=2
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

  ! read dielectric matrix and invert for reduced q-point set
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call getunit(un)
  allocate(scrni(ngqmax,ngqmax,nqptr))
  scrni(:,:,:)=zzero
  !-----------------------------------!
  !     loop over reduced q-points    !
  !-----------------------------------!
  do iqr=1,nqptr
     ! locate reduced q-point in non-reduced set
     iqrnr=iplocnr(ivqr(1,iqr),ngridq)
     n=ngq(iqrnr)
     ! obtain inverse of dielectric matrix
     call geniscreen(iqr,ngqmax,n,scrni(1,1,iqr))
  end do

  !---------------------------------------!
  !     loop over non-reduced q-points    !
  !---------------------------------------!
  do iq=1,nqpt
     write(*,*) 'radial integrals for q-point:',iq
     call putematrad(iq)
  end do

  allocate(done(nqpt))
  allocate(phf(ngqmax,ngqmax),potcl(ngqmax,ngqmax))
  ! allocate array to keep values for all k-points in next loop
  allocate(sccli(nst1,nst3,nst2,nst4,nkptnr))
  allocate(emat12k(nst1,nst3,ngq(1)),emat12kp(nst1,nst3,ngq(1)))
  phf(:,:)=zzero
  potcl(:,:)=0.d0
  sccli(:,:,:,:,:)=zzero
  done(:)=.false.
  ikkp=0

  call genfilname(basename='SCCLI',dotext='.OUT',filnam=fname)
  call getunit(un)
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4, &
       sccli(:,:,:,:,1)
  open(un,file=trim(fname),form='unformatted',action='write', &
       status='replace',access='direct',recl=recl)


write(*,*) 'shape(sccli)',shape(sccli)
write(*,*) 'record length for SCI',recl

  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  ! first k-point
  do iknr=1,nkptnr

     emattype=2
     call ematbdcmbs(emattype)

     ! second k-point
     do jknr=iknr,nkptnr

        call cpu_time(cpu2)
        cpu_init1xs=0.d0
        cpu_ematrad=0.d0
        cpu_ematqalloc=0.d0
        cpu_ematqk1=0.d0
        cpu_ematqdealloc=0.d0
        cpu_clph=0.d0
        cpu_suma=0.d0
        cpu_write=0.d0

        ikkp=ikkp+1
        iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
        iv(:)=modulo(iv(:),ngridk(:))
        ! q-point (reduced)
        iqr=iqmapr(iv(1),iv(2),iv(3))
        vqr(:)=vqlr(:,iqr)
        ! q-point (non-reduced)
        iq=iqmap(iv(1),iv(2),iv(3))
        tq0=tqgamma(iq)
        vq(:)=vql(:,iq)
        ! locate reduced q-point in non-reduced set
        iqrnr=iplocnr(ivqr(1,iqr),ngridq)
        ! local field effects size
        n=ngq(iq)

        call findsymeqiv(vq,vqr,nsc,sc,ivgsc)
        call findgqmap(iq,iqr,nsc,sc,ivgsc,ngqmax,n,jsym,jsymi,ivgsym,igqmap)

        ! store symmetry
        s(:,:)=dble(symlat(:,:,lsplsymc(jsym)))
        ! store inverse of symmetry
        si(:,:)=dble(symlat(:,:,lsplsymc(jsymi)))

        ! cross check symmetry relation (q1 = s^-1 * q + G_s)
        if (sum(abs(vq-(matmul(transpose(s),vqr)+dble(ivgsym)))).gt.epslat) then
           write(*,*) 'deviation:',iknr,jknr,iqr,iq,v2
           write(*,*)
           write(*,'("Error(",a,"): cross checking of symmetry reduction for &
                &q-vector failed:")') trim(thisnam)
           write(*,'(" non-reduced q-point                    :",i8)') iq
           write(*,'(" reduced q-point                        :",i8)') iqr
           write(*,'(" reduced q-point in non-reduced set     :",i8)') &
                iqrnr
           write(*,'(" crystal symmetry operation             :",i8)') isym
           write(*,'(" umklapp G-vector                       :",3i8)') ivgsym
           write(*,*)
           call terminate
        end if

        ! temporary arrays
        allocate(tm(n,n),tmi(n,n),emat12(nst12,n),emat34(nst34,n))
        allocate(scclit(nst34,nst12))

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
!!$write(9000+iq,'(4i6,3g18.10)') iqr,igq1,igq2,j2,tmi(igq1,igq2),abs(tmi(igq1,igq2))
           end do
        end do

        call cpu_time(cpu0)

        ! set up Coulomb potential
        do igq1=1,n
           do igq2=igq1,n
              ! calculate weights for Coulomb potential
              iflg=0
              if (tq0.and.((igq1.eq.1).or.(igq2.eq.1))) then
                 ! consider only 1/q and 1/q^2 cases for q goint to zero
                 iflg=bsediagweight
              else if ((igq1.eq.1).and.(igq2.eq.1)) then
                 ! consider only 1/q^2 cases for non-zero q-point
                 iflg=bsediagweight
              end if
              call genwiq2xs(iflg,iq,igq1,igq2,potcl(igq1,igq2))
              potcl(igq2,igq1)=potcl(igq1,igq2)
              
              !if (iflg.ne.0) &
              !write(50,'(a,6i8,2g18.10)') 'ik,jk,q,bsediagweight,g,gp,potcl',iknr,jknr,iq,iflg,igq1,igq2,potcl(igq1,igq2),fourpi/(gqc(igq1,iq)*gqc(igq2,iq))


              ! end loop over (G,Gp)-vectors
           end do
        end do
	call cpu_time(cpu1)
	cpu_clph=cpu_clph+cpu1-cpu0

        call genfilname(iq=iq,dotext='_SCI.OUT',setfilext=.true.)
        if (.not.done(iq)) call writegqpts(iq)
        call genfilname(dotext='_SCR.OUT',setfilext=.true.)

        emattype=2
        call ematbdcmbs(emattype)

        ! calculate matrix elements of the plane wave
        call cpu_time(cpu0)
        call init1xs(qvkloff(1,iq))
        call cpu_time(cpu1)
        cpu_init1xs=cpu_init1xs+cpu1-cpu0

        write(*,'(a,i6,2x,2i5,2x,2i5,2x,i6)') 'ikkp,iknr,jknr,iq,iqr,ngq(iq)',&
             ikkp,iknr,jknr,iq,iqr,ngq(iq)

        call cpu_time(cpu0)
        call getematrad(iq)
        call cpu_time(cpu1)
        cpu_ematrad=cpu_ematrad+cpu1-cpu0

        call cpu_time(cpu0)
        call ematqalloc
        call cpu_time(cpu1)
        cpu_ematqalloc=cpu_ematqalloc+cpu1-cpu0

        call cpu_time(cpu0)
        emattype=2
        ! calculate matrix elements of the plane wave
        call ematqk1(iq,iknr)
	call ematbdcmbs(emattype) !!! ***
        call cpu_time(cpu1)
        cpu_ematqk1=cpu_ematqk1+cpu1-cpu0

        call cpu_time(cpu0)
        call ematqdealloc
        call cpu_time(cpu1)
        cpu_ematqdealloc=cpu_ematqdealloc+cpu1-cpu0
        call cpu_time(cpu0)

        ! help arrays h1(cc',G) = M_G(kcc'), h2(G',vv') = conjg(M_G'(kvv'))
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


        select case (trim(screentype))
        case('longrange')
           ! long range screening
	   tm(:,:)=zzero
           ! keep (0,0)-element
	   tm(1,1)=tmi(1,1)*potcl(1,1)
           do igq1=2,n
              tm(igq1,igq1)=fourpi*tmi(1,1)/gqc(igq1,iq)**2
           end do
        case('diag')
           ! only diagonal of screening
           tm(:,:)=potcl(1:n,1:n)*tmi(:,:)
           forall(igq1=1:n,igq2=1:n,igq1.ne.igq2)
              tm(igq1,igq2)=zzero
           end forall
	case('full')
	   ! full screening
	   !!tm(:,:)=phf(:,:)*potcl(:,:)*tmi(:,:)
           tm(:,:)=potcl(1:n,1:n)*tmi(:,:)
        end select
        
!write(1234,*) 'ikkp',ikkp,tm
!write(1234,*)


        ! * version 1
!!	scclit=matmul(conjg(emat34),matmul(tm,transpose(emat12)))/omega/nkptnr

        ! * version 2 : like in calkWD.frc and SELF documentation
scclit=matmul(emat34,matmul(transpose(tm),transpose(conjg(emat12))))/omega/nkptnr

write (*,'(a,3g18.10)') 'upper limits:',maxval(abs(tm)),maxval(abs(emat12)),maxval(abs(emat34))

        ! * version 3 like in pep-thesis
!!!	scclit=matmul(emat34,matmul(tm,transpose(conjg(emat12))))/omega/nkptnr

        sccli(:,:,:,:,jknr)=zzero
!!$        do igq1=1,n
!!$           do igq2=1,n
!!$              do ist1=1,nst1
!!$                 do ist3=1,nst3
!!$                    do ist2=1,nst2
!!$                       do ist4=1,nst4
!!$                          sccli(ist1,ist3,ist2,ist4,jknr)= &
!!$                               sccli(ist1,ist3,ist2,ist4,jknr)+ &
!!$                               tm(igq1,igq2)* &
!!$                               conjg(xiou(ist1,ist2,igq2))* &
!!$                               xiuo(ist3,ist4,igq1) /omega/nkptnr
!!$                       end do
!!$                    end do
!!$                 end do
!!$              end do
!!$              ! end loop over (G,Gp) pairs
!!$           end do
!!$        end do

        ! map back to individual band indices
        j2=0
        do ist4=1,nst4
           do ist3=1,nst3
              j2=j2+1
              j1=0
              do ist2=1,nst2
                 do ist1=1,nst1
                    j1=j1+1
                    sccli(ist1,ist3,ist2,ist4,jknr)=scclit(j2,j1)
                 end do
              end do
           end do
        end do

        call cpu_time(cpu1)
        cpu_suma=cpu_suma+cpu1-cpu0
        call cpu_time(cpu0)

	do ist1=1,nst1
	   do ist3=1,nst3
	      do ist2=1,nst2
		 do ist4=1,nst4
		    write(1100,'(i5,3x,3i4,2x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
			 ist3,jknr,ist2,ist4,sccli(ist1,ist3,ist2,ist4,jknr),&
			 abs(sccli(ist1,ist3,ist2,ist4,jknr))
		 end do
	      end do
	   end do
	end do

        ! write screened Coulomb interaction to direct-access file
        write(un,rec=ikkp) ikkp,iknr,jknr,iq,iqr,nst1,nst3,nst4,nst2, &
             sccli(:,:,:,:,jknr)

        call cpu_time(cpu1)
        cpu_write=cpu_write+cpu1-cpu0

        call genfilname(dotext='_SCI.OUT',setfilext=.true.)


        done(iq)=.true.
        deallocate(tm,tmi,emat12,emat34,scclit)

        call cpu_time(cpu3)
        t3=cpu_ematqdealloc+cpu_ematqk1+cpu_ematqalloc+cpu_ematrad+cpu_init1xs+cpu_clph+cpu_suma+cpu_write
        write(*,'(a,f12.3)') 'init1xs     :',cpu_init1xs
        write(*,'(a,f12.3)') 'ematrad     :',cpu_ematrad
        write(*,'(a,f12.3)') 'ematqalloc  :',cpu_ematqalloc
        write(*,'(a,f12.3)') 'ematqk1     :',cpu_ematqk1
        write(*,'(a,f12.3)') 'ematqdealloc:',cpu_ematqdealloc
        write(*,'(a,f12.3)') 'ph+cl       :',cpu_clph
        write(*,'(a,f12.3)') 'summation   :',cpu_suma
        write(*,'(a,f12.3)') 'write       :',cpu_write
        write(*,'(a,f12.3)') '*** sum     :',t3
        write(*,'(a,f12.3)') '*** rest    :',cpu3-cpu2-t3
        write(*,'(a,f12.3)') '*** overall :',cpu3-cpu2
        write(*,*)

        ! end loop over (k,kp) pairs
     end do     
  end do
  close(un)


  call findgntn0_clear
  deallocate(emat12k,emat12kp)
  deallocate(done,scrni,phf,potcl,sccli)

  !--------------!
  !   finalize   !
  !--------------!
  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screened Coulomb interaction&
       & finished"
end subroutine scrcoulint
