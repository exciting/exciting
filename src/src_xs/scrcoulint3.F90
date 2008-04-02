
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint3
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
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:)
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

  !---------------------------------------!
  !     loop over non-reduced q-points    !
  !---------------------------------------!
  do iq=1,nqpt
     write(*,*) 'radial integrals for q-point:',iq
     call putematrad(iq)
  end do

  allocate(phf(ngqmax,ngqmax),potcl(ngqmax,ngqmax))
  allocate(sccli(nst1,nst3,nst2,nst4))
  allocate(emat12k(nst1,nst3,ngq(1)),emat12kp(nst1,nst3,ngq(1)))
  potcl(:,:)=0.d0
  sccli(:,:,:,:)=zzero
  

  call genfilname(basename='SCCLI',dotext='.OUT',filnam=fname)
  call getunit(un)
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4, &
       sccli(:,:,:,:)
  open(un,file=trim(fname),form='unformatted',action='write', &
       status='replace',access='direct',recl=recl)


write(*,*) 'shape(sccli)',shape(sccli)
write(*,*) 'record length for SCI',recl
  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  ikkp=0
  ! first k-point
  do iknr=1,nkptnr

     ! second k-point
     do jknr=iknr,nkptnr

        ikkp=ikkp+1

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

        allocate(emat12(nst12,n),emat34(nst34,n))
        allocate(scclit(nst34,nst12))

        write(*,'(i5,3x,2i5,2x,i5,2x,2i5)') ikkp,iknr,jknr,iq,iqr,iqrnr


        ! Coulomb potential
        do igq1=1,n
           if ((gqc(igq1,iq).gt.epslat)) then
              potcl(igq1,igq1)=fourpi/gqc(igq1,iq)**2
           else
              call genwiq2xs(1,iq,igq1,igq1,potcl(igq1,igq1))
           end if
        end do

        emattype=2
        call ematbdcmbs(emattype)
        call getematrad(iq)
        call ematqalloc
        ! calculate matrix elements of the plane wave
        call ematqk1(iq,iknr)
        emattype=2
        call ematbdcmbs(emattype)

        sccli(:,:,:,:)=zzero
        do igq1=1,n
!           do igq2=1,n
              do ist1=1,nst1
                 do ist3=1,nst3
                    do ist2=1,nst2
                       do ist4=1,nst4
                          sccli(ist1,ist3,ist2,ist4)= &
                               sccli(ist1,ist3,ist2,ist4)+ &
                               potcl(igq1,igq1)/14.d0*  &
                               conjg(xiou(ist1,ist2,igq1))* &
                               xiuo(ist3,ist4,igq1) /omega/nkptnr
                       end do
                    end do
                 end do
              end do
              ! end loop over (G,Gp) pairs
 !          end do
        end do

	do ist1=1,nst1
	   do ist3=1,nst3
	      do ist2=1,nst2
		 do ist4=1,nst4
		    write(1100,'(i5,3x,3i4,2x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
			 ist3,jknr,ist2,ist4,sccli(ist1,ist3,ist2,ist4),&
			 abs(sccli(ist1,ist3,ist2,ist4))
		 end do
	      end do
	   end do
	end do

        ! write screened Coulomb interaction to direct-access file
        write(un,rec=ikkp) ikkp,iknr,jknr,iq,iqr,nst1,nst3,nst4,nst2, &
             sccli(:,:,:,:)

        deallocate(emat12,emat34)
        deallocate(scclit)


        ! end loop over k,kp
     end do
  end do

  close(un)


  call findgntn0_clear

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

end subroutine scrcoulint3
