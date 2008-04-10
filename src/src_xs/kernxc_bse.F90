
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine kernxc_bse(oct)
  !
  ! BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
  !
  use modmain
  use modmpi
  use modtetra
  use modxs
  use modfxcifc
  use invert
  use m_tdgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genwgrid
  use m_dyson
  use m_tdzoutpr3
  use m_getpemat
  use m_getx0
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: oct
  ! local variables
  character(*), parameter :: thisnam = 'kernxs_bse'
  integer, parameter :: iqmt=1
  real(8), parameter :: delt=1.d-5
  character(256) :: filnam,filnam2,filnam3,filnam4
  complex(8),allocatable :: fxc(:,:,:), idf(:,:), mdf1(:),w(:), chi0hd(:)
  complex(8),allocatable :: chi0h(:),chi0wg(:,:,:),chi0(:,:),chi0i(:,:)
  complex(8),allocatable :: chi0h2(:),chi0wg2(:,:,:),chi02(:,:),chi02i(:,:)
  integer :: n,m,recl,j,iw,wi,wf,nwdfp,nc,octh,oct1,oct2,igmt
  integer, external :: l2int
  complex(8) :: zt1,bsediagshift,bsediagshiftc
  integer :: sh(2),ig
  ! ***    
  character(256) :: fname
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iknrq,jknrq,iqr,iq,iqrnr,isym,jsym,jsymi,igq1,igq2,iflg
  integer :: ngridkt(3),iv(3),ivgsym(3),un,un2,un3,j1,j2,iws1,iws2,a1,a2,a3
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp,ikkph
  integer :: ikkp_,iknr_,jknr_,iq_,iqr_,nst1_,nst2_,nst3_,nst4_
  logical :: nosymt,reducekt,tq0,nsc,tphf,ksintp,tks1,tks2
  real(8) :: vklofft(3),vqr(3),vq(3),v2(3),s(3,3),si(3,3),t1,t2,t3,ws,ws1,ws2,tp
  real(8), allocatable :: potcl(:,:),dek(:,:),dekp(:,:),dok(:,:),dde(:,:)
  real(8), allocatable :: dokp(:,:),scisk(:,:),sciskp(:,:)
  real(8), allocatable :: zmr(:,:),zmq(:,:),deval(:,:,:),docc(:,:,:),scis(:,:,:)
  integer :: igqmap(maxsymcrys),sc(maxsymcrys),ivgsc(3,maxsymcrys)
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:),scclih(:,:,:,:)
  complex(8), allocatable :: scrni(:,:,:),tm(:,:),tmi(:,:),den1(:),den2(:)
  complex(8), allocatable :: phf(:,:),emat12(:,:),emat12p(:,:),emat(:,:,:,:)
  complex(8), allocatable :: residr(:,:),residq(:,:),osca(:,:),oscb(:,:),hdg(:,:,:)
  logical, allocatable :: done(:)
  ! external functions
  integer, external :: iplocnr,idxkkp
  logical, external :: tqgamma

  real(8) :: cpu0,cpu1,cpu2,cpu3
  real(8) :: cpu_init1xs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1
  real(8) :: cpu_ematqdealloc,cpu_clph,cpu_suma,cpu_write
  complex(8), allocatable :: emat12k(:,:,:),emat12kp(:,:,:)

  call genfilname(setfilext=.true.)

t3=1.d0
bsediagshift=zzero
bsediagshiftc=zzero

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
  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated&
       & within lmax values:', lmaxapw,lmaxemat,lmaxapw
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
  call getunit(un)

write(*,*) 'nst1,2,3,4',nst1,nst2,nst3,nst4

  call genparidxran('w')
  ! sampling type for Brillouin zone sampling
  bzsampl=l2int(tetra)
  ! limits for w-points
  wi=wpari
  wf=wparf
  nwdfp=wparf-wpari+1
  ! matrix size for local field effects (first q-point is Gamma-point)
  n=ngq(iqmt)

  ! allocate global arrays
  if (allocated(xiou)) deallocate(xiou)
  allocate(xiou(nst1,nst3,n))
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3,nst1,nst3))
  if (allocated(deou)) deallocate(deou)
  allocate(deou(nst1,nst3))
  if (allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1,nst3))
  ! allocate local arrays
  allocate(emat12(nst13,n),emat12p(nst13,n),zmr(nst13,nst13),zmq(nst13,nst13))
  allocate(dek(nst1,nst3),dekp(nst1,nst3),dde(nst1,nst3))
  allocate(dok(nst1,nst3),dokp(nst1,nst3))
  allocate(scisk(nst1,nst3),sciskp(nst1,nst3))
  allocate(fxc(n,n,nwdf))
  allocate(sccli(nst1,nst3,nst2,nst4))
  allocate(scclih(nst1,nst3,nst2,nst4))
  allocate(scclit(nst13,nst13))
  allocate(emat12k(n,nst1,nst3),emat12kp(nst1,nst3,n))
  allocate(residr(nst13,n),residq(nst13,n))
  allocate(w(nwdf),osca(n,n),oscb(n,n),den1(nwdf),den2(nwdf))
  allocate(chi0i(n,n),chi0(n,n),chi0wg(n,2,3),chi0h(9))
  allocate(chi02i(n,n),chi02(n,n),chi0wg2(n,2,3),chi0h2(9))
  fxc(:,:,:)=zzero
  sccli(:,:,:,:)=zzero
  allocate(emat(nst1,nst3,n,nkptnr))
  allocate(deval(nst1,nst3,nkptnr))
  allocate(docc(nst1,nst3,nkptnr))
  allocate(scis(nst1,nst3,nkptnr))

allocate(hdg(nst1,nst3,nkptnr))
hdg=zzero

  ! generate energy grid
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)

  ! open file for screened Coulomb interaction
  call genfilname(basename='SCCLI',dotext='.OUT',filnam=fname)
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4, &
       sccli(:,:,:,:)
  open(un,file=trim(fname),form='unformatted',action='read', &
       status='old',access='direct',recl=recl)

  ! precalculate matrix elements
  emattype=1
  call ematbdcmbs(emattype)
  call ematrad(iqmt)
  call ematqalloc
  do iknr=1,nkptnr
     iknrq=ikmapikq(iknr,iqmt)
     write(*,*) 'generation of matrix elements: k-point:',iknr
     ! matrix elements for k and q=0
     call ematqk1(iqmt,iknr)
     emat(:,:,:,iknr)=xiou(:,:,:)
     deallocate(xiou,xiuo)
     call getdevaldoccsv(iqmt,iknr,iknrq,istlo1,isthi1,istlo2,isthi2,deou, &
          docc12,scisk)
     deval(:,:,iknr)=deou(:,:)
     docc(:,:,iknr)=docc12(:,:)
     scis(:,:,iknr)=scisk(:,:)
  end do
  emattype=2
  call ematbdcmbs(emattype)


  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  ikkp=0
  ! first k-point
  do iknr=1,nkptnr
     iknrq=ikmapikq(iknr,iqmt)

!<new>
     emattype=1
     call ematbdcmbs(emattype)
     allocate(xiou(nst1,nst3,n))
     xiou(:,:,:)=emat(:,:,:,iknr)
     deou(:,:)=deval(:,:,iknr)
     docc12(:,:)=docc(:,:,iknr)

 write(7712,*) xiou(2,3,1),xiou(2,3,7),xiou(1,2,1),xiou(1,2,6)
 write(8712,*) deou(2,3),deou(3,2),deou(1,2),deou(2,1)

     ! apply gauge wrt. symmetrized Coulomb potential
     call getpemat(iqmt,iknr,'PMAT_SCR.OUT','',m12=xiou,p12=pmou)
 write(6712,*) pmou(oct,2,3),pmou(oct,1,2)
     dek(:,:)=deou(:,:)
     dok(:,:)=docc12(:,:)
     scisk(:,:)=scis(:,:,iknr)
     ! assign optical component
     xiou(:,:,1)=pmou(oct,:,:)
     do igq1=1,n
        emat12k(igq1,:,:)=xiou(:,:,igq1)
     end do
!</new>


!!$     ! matrix elements for k and q=0
!!$     emattype=1
!!$     call ematbdcmbs(emattype)
!!$     call init1xs(qvkloff(1,iqmt))
!!$     call getdevaldoccsv(iqmt,iknr,iknrq,istlo1,isthi1,istlo2,isthi2,deou, &
!!$          docc12,scisk)
!!$     dek(:,:)=deou(:,:)
!!$     dok(:,:)=docc12(:,:)
!!$     call getematrad(iqmt)
!!$     call ematqalloc
!!$!     call cpu_time(cpu0)
!!$     call ematqk1(iqmt,iknr)
!!$!     call cpu_time(cpu1)
!!$!     cpu_ematqk1=cpu_ematqk1+cpu1-cpu0
!!$     ! apply gauge wrt. symmetrized Coulomb potential
!!$     call getpemat(iqmt,iknr,'PMAT_SCR.OUT','',m12=xiou,p12=pmou)
!!$     ! assign optical component
!!$     xiou(:,:,1)=pmou(oct,:,:)
!!$     do igq1=1,n
!!$        emat12k(igq1,:,:)=xiou(:,:,igq1)
!!$     end do



     deallocate(xiou)
     emattype=2
     call ematbdcmbs(emattype)

     residr(:,:)=zzero
     residq(:,:)=zzero
     ! second k-point
     do jknr=1,nkptnr
        jknrq=ikmapikq(jknr,iqmt)

        call cpu_time(cpu2)
        cpu_init1xs=0.d0
        cpu_ematrad=0.d0
        cpu_ematqalloc=0.d0
        cpu_ematqk1=0.d0
        cpu_ematqdealloc=0.d0
        cpu_clph=0.d0
        cpu_suma=0.d0
        cpu_write=0.d0

        if (iknr.le.jknr) then
           ! index for upper triangle
           ikkp=idxkkp(iknr,jknr,nkptnr)
        else
           ! swapped index for lower triangle
           ikkp=idxkkp(jknr,iknr,nkptnr)
        end if

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

   write(*,'(a,i6,2x,2i5,2x,2i5)') 'ikkp,iknr,jknr,iq,iqr',&
             ikkp,iknr,jknr,iq,iqr



!<new>
        emattype=1
        call ematbdcmbs(emattype)
        allocate(xiou(nst1,nst3,n))
        xiou(:,:,:)=emat(:,:,:,jknr)
        deou(:,:)=deval(:,:,jknr)
        docc12(:,:)=docc(:,:,jknr)

 write(7713,*) xiou(2,3,1),xiou(2,3,7),xiou(1,2,1),xiou(1,2,6)
 write(8713,*) deou(2,3),deou(3,2),deou(1,2),deou(2,1)

        ! apply gauge wrt. symmetrized Coulomb potential
        call getpemat(iqmt,jknr,'PMAT_SCR.OUT','',m12=xiou,p12=pmou)
 write(6713,*) pmou(oct,2,3),pmou(oct,1,2)
        dekp(:,:)=deou(:,:)
        dokp(:,:)=docc12(:,:)
        sciskp(:,:)=scis(:,:,jknr)
        ! assign optical component
        xiou(:,:,1)=pmou(oct,:,:)
        emat12kp(:,:,:)=xiou(:,:,:)
!</new>




!!$        ! matrix elements for kp and q=0
!!$        emattype=1
!!$        call ematbdcmbs(emattype)
!!$        call init1xs(qvkloff(1,iqmt))
!!$        call getdevaldoccsv(iqmt,jknr,jknrq,istlo1,isthi1,istlo2,isthi2,deou, &
!!$             docc12,sciskp)
!!$        dekp(:,:)=deou(:,:)
!!$        dokp(:,:)=docc12(:,:)
!!$        call getematrad(iqmt)
!!$        call ematqalloc
!!$        call cpu_time(cpu0)
!!$        call ematqk1(iqmt,jknr)
!!$        call cpu_time(cpu1)
!!$        cpu_ematqk1=cpu_ematqk1+cpu1-cpu0
!!$        ! apply gauge wrt. symmetrized Coulomb potential
!!$        call getpemat(iqmt,jknr,'PMAT_SCR.OUT','',m12=xiou,p12=pmou)
!!$        ! assign optical component
!!$        xiou(:,:,1)=pmou(oct,:,:)
!!$        emat12kp(:,:,:)=xiou(:,:,:)

        deallocate(xiou)
        emattype=2
        call ematbdcmbs(emattype)

!!!!!!!!!!!!!!!!!!!!!!!! CHECK
        write(7711,*) emat12kp(2,3,1),emat12kp(2,3,7),emat12kp(1,2,1),emat12kp(1,2,6)


!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        do igq1=1,n
!!$           do ist1=1,nst1
!!$              do ist2=1,nst3
!!$                 write(4000,'(4i5,3g18.10)') iknr,igq1,ist1,ist2, &
!!$                      emat12k(igq1,ist1,ist2),abs(emat12k(igq1,ist1,ist2))**2
!!$              end do
!!$           end do
!!$        end do
!!$        do igq1=1,n
!!$           do ist1=1,nst1
!!$              do ist2=1,nst3
!!$                 write(5000,'(4i5,3g18.10)') jknr,igq1,ist1,ist2, &
!!$                      emat12kp(ist1,ist2,igq1),abs(emat12kp(ist1,ist2,igq1))**2
!!$              end do
!!$           end do
!!$        end do

        ! get screened Coulomb interaction
        if (iknr.le.jknr) then
           ! read from file
           read(un,rec=ikkp) ikkp_,iknr_,jknr_,iq_,iqr_,nst1_,nst2_,nst3_, &
                nst4_,sccli


!!$write(*,*) '====nst1,2,3,4',nst1,nst2,nst3,nst4
!!$
!!$
!!$           if ((ikkp.ne.ikkp_).or.(iknr.ne.iknr_).or.(jknr.ne.jknr_).or. &
!!$                (iq.ne.iq_).or.(iqr.ne.iqr_).or.(nst1.ne.nst1_).or. &
!!$                (nst2.ne.nst2_).or.(nst3.ne.nst3_).or.(nst4.ne.nst4_)) then
!!$              write(*,*)
!!$              write(*,'("Error(kernxc_bse): wrong indices for screened Coulomb&
!!$                   & interaction")')
!!$              write(*,'(" indices (ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,&
!!$                   &nst4)")')
!!$              write(*,'(" current:",i6,3x,2i4,2x,2i4,2x,4i4)') ikkp,iknr,jknr,&
!!$                   iq,iqr,nst1,nst2,nst3,nst4
!!$              write(*,'(" file   :",i6,3x,2i4,2x,2i4,2x,4i4)') ikkp_,iknr_,&
!!$                   jknr_,iq_,iqr_,nst1_,nst2_,nst3_,nst4_
!!$              write(*,*)
!!$              call terminate
!!$           end if
        else
           write(*,*) 'using Hermitian transposed scr. Coul. int.'
           ! use Hermitian property for lower triangle
           read(un,rec=ikkp) ikkp_,iknr_,jknr_,iq_,iqr_,nst1_,nst2_,nst3_, &
                nst4_,scclih
           do ist1=1,nst1
              do ist2=1,nst3
                 sccli(ist1,ist2,:,:)=conjg(scclih(:,:,ist1,ist2))
              end do
           end do
        end if
        ! diagonal of BSE-kernel (approximate by first value in matrix)
        ! *** improve later
!        if (ikkp.eq.1) then
!           bsediagshift=dble(sccli(1,1,1,1))
!        end if

        ! rescale
        sccli = - sccli

        ! set diagonal of Bethe-Salpeter kernel to zero (cf. A. Marini, PRL 2003)
        if (iknr.eq.jknr) then
           do ist3=1,nst3
              do ist1=1,nst1
	      
	      zt1=sccli(ist1,ist3,ist1,ist3)
!!!	      hdg(ist1,ist3,iknr)=-zt1
	      t1=dble(zt1)
	      t2=aimag(zt1)
    tp=atan2(t2,t1)
    if (tp.lt.0.d0) tp=tp+twopi
    write(1107,'(3i5,2g18.10,2x,2g18.10)') iknr,ist1,ist3,zt1,abs(zt1),tp*180.d0/pi
	      
                 sccli(ist1,ist3,ist1,ist3)= zzero
              end do
           end do
        end if
	
!!$if (iknr.le.jknr) then
!!$        	! * write out screened Coulomb interaction
!!$	do ist1=1,nst1
!!$	   do ist3=1,nst3
!!$	      do ist2=1,nst2
!!$		 do ist4=1,nst4
!!$		    write(1101,'(i5,3x,3i4,2x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
!!$			 ist3,jknr,ist2,ist4,sccli(ist1,ist3,ist2,ist4),&
!!$			 abs(sccli(ist1,ist3,ist2,ist4))
!!$		 end do
!!$	      end do
!!$	   end do
!!$	end do
!!$end if
        j1=0
        do ist2=1,nst3
           do ist1=1,nst1
              j1=j1+1
              emat12p(j1,:)=conjg(emat12kp(ist1,ist2,:))
           end do
        end do
        ! map 
        j2=0
        do ist3=1,nst3
           do ist1=1,nst1
              j2=j2+1
              j1=0
              do ist4=1,nst3
                 do ist2=1,nst1
                    j1=j1+1
                    zt1=sccli(ist1,ist3,ist2,ist4)
                    ! four point energy difference
                    t1=dekp(ist2,ist4)-dek(ist1,ist3)
                    ! arrays for R- and Q-residuals
                    if (abs(t1).ge.delt) then
                       zmr(j2,j1)=zt1/t1
                       zmq(j2,j1)=0.d0
                    else
                       zmr(j2,j1)=0.d0
                       zmq(j2,j1)=zt1
                    end if
                 end do
              end do
           end do
        end do

        ! calculate residual "R" 
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))
        residr=residr+matmul(zmr,emat12p)

        ! calculate residual "Q" 
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))
        residq=residq+matmul(zmq,emat12p)

!!$        call cpu_time(cpu3)
!!$        t3=cpu_ematqdealloc+cpu_ematqk1+cpu_ematqalloc+cpu_ematrad+cpu_init1xs+cpu_clph+cpu_suma+cpu_write
!!$        write(*,'(a,f12.3)') 'init1xs     :',cpu_init1xs
!!$        write(*,'(a,f12.3)') 'ematrad     :',cpu_ematrad
!!$        write(*,'(a,f12.3)') 'ematqalloc  :',cpu_ematqalloc
!!$        write(*,'(a,f12.3)') 'ematqk1     :',cpu_ematqk1
!!$        write(*,'(a,f12.3)') 'ematqdealloc:',cpu_ematqdealloc
!!$        write(*,'(a,f12.3)') 'summation   :',cpu_suma
!!$        write(*,'(a,f12.3)') 'write       :',cpu_write
!!$        write(*,'(a,f12.3)') '*** sum     :',t3
!!$        write(*,'(a,f12.3)') '*** rest    :',cpu3-cpu2-t3
!!$        write(*,'(a,f12.3)') '*** overall :',cpu3-cpu2
        write(*,*)

        ! end inner loop over k-points
     end do

     !--------------------------!
     !     set up BSE-kernel    !
     !--------------------------!

     do ist3=1,nst3
        do ist1=1,nst1
           osca(:,:)=zzero
           oscb(:,:)=zzero
           j1=ist1+(ist3-1)*nst1
           ! set up inner part of kernel
           
           ! generate oscillators
           call tdzoutpr3(n,n,zone,emat12k(:,ist1,ist3),residr(j1,:),osca)
           ! add Hermitian transpose
           osca=osca+conjg(transpose(osca))
           call tdzoutpr3(n,n,zone,emat12k(:,ist1,ist3),residq(j1,:),oscb)

!           ! set up energy denominators
!           den1(:)=2.d0/(w(:)+bsediagshift+dek(ist1,ist3)+zi*brdtd)/nkptnr/&
!                omega          
!           den2(:)=2.d0/(w(:)+bsediagshift+dek(ist1,ist3)+zi*brdtd)**2/nkptnr/&
!                omega

!           ! set up energy denominators
!           den1(:)=2.d0/(w(:)+bsediagshift+dek(ist1,ist3)+zi*brdtd) + &
!                2.d0/(-w(:)-bsediagshift+dek(ist1,ist3)-zi*brdtd)
!           den2(:)=2.d0/(w(:)+bsediagshift+dek(ist1,ist3)+zi*brdtd)**2 + &
!                2.d0/(-w(:)-bsediagshift+dek(ist1,ist3)-zi*brdtd)**2
!           den1=den1/nkpt/omega
!           den2=den2/nkpt/omega

!!$! *** this part is working for Si_lapw and Si_APW+lo ***
!!$           ! set up energy denominators
!!$           den1(:)=2.d0/(w(:)+hdg(ist1,ist3,iknr)+dek(ist1,ist3)+zi*brdtd)
!!$           den2(:)=2.d0/(w(:)+hdg(ist1,ist3,iknr)+dek(ist1,ist3)+zi*brdtd)**2
!!$           den1=den1/nkpt/omega
!!$           den2=den2/nkpt/omega
!!$! *** end

           ! set up energy denominators
           den1(:)=2.d0/(w(:)+hdg(ist1,ist3,iknr)+dek(ist1,ist3)+zi*brdtd) + &
                2.d0/(-w(:)+hdg(ist1,ist3,iknr)+dek(ist1,ist3)-zi*brdtd)
           den2(:)=2.d0/(w(:)+hdg(ist1,ist3,iknr)+dek(ist1,ist3)+zi*brdtd)**2 +&
                2.d0/(-w(:)+hdg(ist1,ist3,iknr)+dek(ist1,ist3)-zi*brdtd)**2
           den1=den1/nkpt/omega
           den2=den2/nkpt/omega
           ! update kernel
           do iw=1,nwdf
              fxc(:,:,iw)=fxc(:,:,iw)+osca(:,:)*den1(iw)+oscb(:,:)*den2(iw)
           end do

           ! end loop over states #1
        end do
        ! end loop over states #3
     end do
     
     ! end outer loop over k-points
  end do
  close(un)


!!!! write BSE diagonal to file
write(1108) hdg

  do iw=1,nwdf
     do igq1=1,n
        do igq2=1,n
           write(777,'(i6,2x,2i5,2g18.10)') iw,igq1,igq2,fxc(igq1,igq2,iw)
        end do
     end do
  end do


  ! filename for response function file
  call genfilname(basename='X0',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam)

  ! filename for xc-kernel (ASCII)
  call genfilname(basename='FXC_BSE',asc=.true.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam2)
  open(un,file=trim(filnam2),form='formatted',action='write',status='replace')

  call getunit(un2)
  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam3)
  inquire(iolength=recl) fxc(:,:,1)
  open(un2,file=trim(filnam3),form='unformatted',action='write', &
       status='replace',access='direct',recl=recl)
  
  call getunit(un3)
  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE_HEAD',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam4)
  open(un3,file=trim(filnam4),form='formatted',action='write',status='replace')
  
  ! set up kernel
  do iw=1,nwdf
!!$     ! locate shifted energy on grid
!!$     t1=dble(w(iw))-bsediagshift
!!$     ws=1.d0+(t1-wdos(1))*dble(nwdf)/(wdos(2)-wdos(1))
!!$     iws1=floor(ws)
!!$     iws2=ceiling(ws)
!!$     ws1=dble(w(iws1))
!!$     ws2=dble(w(iws2))
!!$
!!$t3=dble(w(iw))
!!$write(*,*) iw,t3,t1,ws,iws1,iws2
!!$
!!$!     if (iws1.lt.1) then
!!$!        write(*,*)
!!$!        write(*,'("Error(",a,"): negative shifted w-point")') &
!!$!             trim(thisnam)
!!$!        write(*,*)
!!$!        call terminate
!!$!     end if
!!$     ksintp=.true.
!!$     tks1=.true.
!!$     tks2=.true.
!!$     if (iws1.eq.iws2) ksintp=.false.
!!$     if (iws1.gt.nwdf) tks1=.false.
!!$     if (iws1.lt.1) then
!!$     	tks1=.false.
!!$	tks2=.false.
!!$     end if
!!$     if ((.not.ksintp).or.(iws2.gt.nwdf)) tks2=.false.
!!$     ! zero in case frequency above the interval is required
!!$     chi0i(:,:)=zzero
!!$     oct1=oct
!!$     oct2=oct
!!$     octh=1+(oct-1)*4
!!$
!!$
!!$write(*,*) 'octs',oct,oct1,oct2,octh
!!$     if (tks1) then
!!$        ! get KS/QP response function for lower w-point
!!$        call getx0(.true.,iqmt,iws1,trim(filnam),'',chi0,chi0wg,chi0h)
!!$        ! head
!!$        chi0(1,1)=chi0h(octh)
!!$        ! wings
!!$        if (n.gt.1) then
!!$           chi0(1,2:)=chi0wg(2:,1,oct1)
!!$           chi0(2:,1)=chi0wg(2:,2,oct2)
!!$        end if
!!$
!!$     do igq1=1,n
!!$        do igq2=1,n
!!$           write(772,'(i6,2x,2i5,2g18.10)') iw,igq1,igq2,chi0(igq1,igq2)
!!$        end do
!!$     end do
!!$
!!$        ! invert
!!$        call zinvert_hermitian(0,chi0,chi0i)
!!$     end if
!!$     if (tks2) then
!!$        ! get KS/QP response function for lower w-point
!!$        call getx0(.true.,iqmt,iws2,trim(filnam),'',chi02,chi0wg2,chi0h2)
!!$        ! head
!!$        chi02(1,1)=chi0h2(octh)
!!$        ! wings
!!$        if (n.gt.1) then
!!$           chi02(1,2:)=chi0wg2(2:,1,oct1)
!!$           chi02(2:,1)=chi0wg2(2:,2,oct2)
!!$        end if
!!$        ! invert
!!$        call zinvert_hermitian(0,chi02,chi02i)
!!$     end if
!!$     if (ksintp) then
!!$        t3=(t1-ws1)/(ws2-ws1)
!!$        ! linear interpolation
!!$        chi0i(:,:)=(1.d0-t3)*chi0i(:,:) + t3*chi02i(:,:)
!!$     end if

     ! simple version without interpolation
     oct1=oct
     oct2=oct
     octh=1+(oct-1)*4
     ! get KS/QP response function for lower w-point
     call getx0(.true.,iqmt,iw,trim(filnam),'',chi0,chi0wg,chi0h)
     ! head
     chi0(1,1)=chi0h(octh)
     ! wings
     if (n.gt.1) then
        chi0(1,2:)=chi0wg(2:,1,oct1)
        chi0(2:,1)=chi0wg(2:,2,oct2)
     end if
     ! invert
     call zinvert_hermitian(0,chi0,chi0i)

     ! multiply inner part of kernel with inverse QP-response function from
     ! both sides
     fxc(:,:,iw)=matmul(chi0i,matmul(fxc(:,:,iw),chi0i))
     ! write kernel to file for each frequency
     write(un2,rec=iw) fxc(:,:,iw)
     do igq1=1,n
        do igq2=1,n
           write(un,'(i6,2x,2i5,2g18.10)') iw,igq1,igq2,fxc(igq1,igq2,iw)
           write(773,'(i6,2x,2i5,2g18.10)') iw,igq1,igq2,chi0i(igq1,igq2)
        end do
     end do
     write(un3,'(i6,2x,g18.10,2x,2g18.10)') iw,dble(w(iw)),fxc(1,1,iw)
  end do
  close(un)
  close(un2)
  close(un3)

  ! deallocate
  !deallocte(..............................)
end subroutine kernxc_bse


integer function idxkkp(ik,ikp,n)
  implicit none
  ! arguments
  integer, intent(in) :: ik,ikp,n
  ! local variables
  integer :: a,s
  if ((ik.le.0).or.(ikp.le.0).or.(n.le.0)) then
     write(*,*)
     write(*,'("Error(idxkkp): negative indices or number of points")')
     write(*,*)
     call terminate
  end if
  if (ik.gt.ikp) then
     write(*,*)
     write(*,'("Error(idxkkp): ik > ikp")')
     write(*,*)
     call terminate
  end if
  s=0
  do a=1,ik-1
     s=s+n-a+1
  end do
  idxkkp=s+ikp-ik+1
end function idxkkp


