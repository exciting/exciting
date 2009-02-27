
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: kernxc_bse
! !INTERFACE:
subroutine kernxc_bse
! !USES:
  use modmain
  use modtetra
  use modxs
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genwgrid
  use m_xszoutpr3
  use m_getpemat
  use m_getunit
  use m_genfilname
! !INPUT/OUTPUT PARAMETERS:
!   oct   : optical diagonal tensor component (in,integer)
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='kernxc_bse'
  integer, parameter :: iqmt=1,noptc=3
  character(256) :: filnam2,filnam3,filnam4
  logical :: tq0
  integer :: iv(3),iw,wi,wf,nwdfp,n,recl,un,un2,un3,j1,j2,oct
  integer :: ikkp,iknr,jknr,iknrq,jknrq,iqr,iq,iqrnr,igq1,igq2
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24
  real(8) :: vqr(3),vq(3),t1,brd
  real(8) :: cpu_init1offs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1
  real(8) :: cpu_ematqdealloc,cpu_clph,cpu_suma,cpu_write
  complex(8) :: zt1
  ! allocatable arrays
  real(8), allocatable :: dek(:,:),dok(:,:),scisk(:,:)
  real(8), allocatable :: dekp(:,:),dokp(:,:),sciskp(:,:)
  real(8), allocatable :: deval(:,:,:),docc(:,:,:),scis(:,:,:)
  real(8), allocatable :: dde(:,:)
  complex(8), allocatable :: zmr(:,:),zmq(:,:),zmra(:,:),zmqa(:,:)
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:),scclih(:,:,:,:)
  complex(8), allocatable :: emat(:,:,:,:),emata(:,:,:,:)
  complex(8), allocatable :: den1(:),den2(:),den1a(:),den2a(:)
  complex(8), allocatable :: emat12p(:,:),emat12pa(:,:)
  complex(8), allocatable :: emat12k(:,:,:),emat12kp(:,:,:)
  complex(8), allocatable :: emat12ka(:,:,:),emat12kpa(:,:,:)
  complex(8), allocatable :: residr(:,:),residq(:,:),osca(:,:),oscb(:,:)
  complex(8), allocatable :: residra(:,:),residqa(:,:),oscaa(:,:),oscba(:,:)
  complex(8), allocatable :: fxc(:,:,:),w(:),bsedg(:,:),bufou(:,:,:),bufuo(:,:,:),pufou(:,:,:),pufuo(:,:,:)
  ! external functions
  integer, external :: idxkkp,l2int
  logical, external :: tqgamma
  brd=broad
  emattype=2
write(unitout,*) 'enty...'; call flushifc(unitout)
  call init0
  call init1
  call init2
write(unitout,*) 'done init2...'; call flushifc(unitout)
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call xssave0
  ! generate Gaunt coefficients
  call xsgauntgen(max(lmaxapw,lolmax),lmaxemat,max(lmaxapw,lolmax))
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(max(lmaxapwwf,lolmax),max(lmaxapwwf,lolmax),lmaxemat,xsgnt)
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

  call genparidxran('w',nwdf)
  ! sampling type for Brillouin zone sampling
  bzsampl=l2int(tetradf)
  ! limits for w-points
  wi=wpari
  wf=wparf
  nwdfp=wparf-wpari+1
  ! matrix size for local field effects (first q-point is Gamma-point)
  n=ngq(iqmt)

  allocate(bufou(nst1,nst3,n))
  allocate(bufuo(nst3,nst1,n))
  allocate(pufou(3,nst1,nst3))
  allocate(pufuo(3,nst3,nst1))
  ! allocate global arrays
  if (allocated(xiou)) deallocate(xiou)
  allocate(xiou(nst1,nst3,n))
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiuo(nst3,nst1,n))
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3,nst1,nst3))
  if (allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3,nst3,nst1))
  if (allocated(deou)) deallocate(deou)
  allocate(deou(nst1,nst3))
  if (allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst3,nst1))
  if (allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1,nst3))
  if (allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst3,nst1))
  ! allocate local arrays
  allocate(emat12p(nst13,-3:n),zmr(nst13,nst13), &
       zmq(nst13,nst13))
  allocate(emat12pa(nst13,-3:n),zmra(nst13,nst13), &
       zmqa(nst13,nst13))
  allocate(dek(nst1,nst3),dekp(nst1,nst3),dde(nst1,nst3))
  allocate(dok(nst1,nst3),dokp(nst1,nst3))
  allocate(scisk(nst1,nst3),sciskp(nst1,nst3))
  allocate(fxc(-3:n,-3:n,nwdf))
  allocate(sccli(nst1,nst3,nst2,nst4))
  allocate(scclih(nst1,nst3,nst2,nst4))
  allocate(scclit(nst13,nst13))
  allocate(emat12k(-3:n,nst1,nst3),emat12kp(nst1,nst3,-3:n))
  allocate(residr(nst13,-3:n),residq(nst13,-3:n))
  allocate(emat12ka(-3:n,nst3,nst1),emat12kpa(nst3,nst1,-3:n))
  allocate(residra(nst13,-3:n),residqa(nst13,-3:n))
  allocate(w(nwdf))
  allocate(osca(-3:n,-3:n),oscb(-3:n,-3:n))
  allocate(oscaa(-3:n,-3:n),oscba(-3:n,-3:n))
  allocate(den1(nwdf),den2(nwdf),den1a(nwdf),den2a(nwdf))
  fxc(:,:,:)=zzero
  sccli(:,:,:,:)=zzero
  allocate(emat(nst1,nst3,n,nkptnr))
  allocate(emata(nst3,nst1,n,nkptnr))
  allocate(deval(nst1,nst3,nkptnr))
  allocate(docc(nst1,nst3,nkptnr))
  allocate(scis(nst1,nst3,nkptnr))
  allocate(bsedg(nst1,nst3))


  if ((fxctype.eq.7).or.(fxctype.eq.8)) then
     call getbsediag
     write(unitout,'("Info(",a,"): read diagonal of BSE kernel")') trim(thisnam)
     write(unitout,'(" mean value : ",2g18.10)') bsed
  end if

  ! generate energy grid
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)

  ! precalculate matrix elements
  emattype=1
  call ematbdcmbs(emattype)
  call ematrad(iqmt)
  call ematqalloc

  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
  do iknr=1,nkptnr
     call chkpt(3,(/task,1,iknr/),'task,sub,k-point; generate matrix elements &
          &of plane wave')    
     iknrq=ikmapikq(iknr,iqmt)
     ! matrix elements for k and q=0
     call ematqk1(iqmt,iknr)
     emat(:,:,:,iknr)=xiou(:,:,:)
     emata(:,:,:,iknr)=xiuo(:,:,:)
     deallocate(xiou,xiuo)
     call getdevaldoccsv(iqmt,iknr,iknrq,istl1,istu1,istl2,istu2,deou, &
          docc12,scisk)
     call getdevaldoccsv(iqmt,iknr,iknrq,istl2,istu2,istl1,istu1,deuo, &
     	  docc21,sciskp)
     deval(:,:,iknr)=deou(:,:)
     docc(:,:,iknr)=docc12(:,:)
     scis(:,:,iknr)=scisk(:,:)
  end do
  emattype=2
  call ematbdcmbs(emattype)

  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  if (allocated(xiou)) deallocate(xiou)
  if (allocated(xiuo)) deallocate(xiuo)
  ikkp=0
  ! first k-point
  do iknr=1,nkptnr
     call chkpt(3,(/task,3,iknr/),'task,sub,k-point; BSE-fxc-kernel')
     iknrq=ikmapikq(iknr,iqmt)

     call getbsedg('BSED.OUT',iknr,nst1,nst3,bsedg)
     !@@@@@@@@@@@@@@@@@@@@@bsedg(:,:)=bsed ! REVERTED TO CONSTANT SHIFT
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     bsedg(:,:)=0.d0

     emattype=1
     call ematbdcmbs(emattype)
     allocate(xiou(nst1,nst3,n))
     allocate(xiuo(nst3,nst1,n))
     xiou(:,:,:)=emat(:,:,:,iknr)
     xiuo(:,:,:)=emata(:,:,:,iknr)
     deou(:,:)=deval(:,:,iknr)
     docc12(:,:)=docc(:,:,iknr)
     ! apply gauge wrt. symmetrized Coulomb potential
     call getpemat(iqmt,iknr,'PMAT_SCR.OUT','',m12=bufou,p12=pufou,m34=bufuo, &
     	p34=pufuo)
     dek(:,:)=deou(:,:)
     dok(:,:)=docc12(:,:)
     ! add BSE diagonal
     scisk(:,:)=scis(:,:,iknr)+bsedg(:,:)

     ! assign optical components
     do oct=1,noptc
       emat12k(-oct,:,:)=pufou(oct,:,:)
       emat12ka(-oct,:,:)=pufuo(oct,:,:)
     end do
     do igq1=1,n
       emat12k(igq1,:,:)=bufou(:,:,igq1)
       emat12ka(igq1,:,:)=bufuo(:,:,igq1)
     end do

     deallocate(xiou,xiuo)
     emattype=2
     call ematbdcmbs(emattype)

     residr(:,:)=zzero
     residq(:,:)=zzero
     residra(:,:)=zzero
     residqa(:,:)=zzero
     ! second k-point
    do jknr=1,nkptnr
        jknrq=ikmapikq(jknr,iqmt)

        cpu_init1offs=0.d0
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

        emattype=1
        call ematbdcmbs(emattype)
        allocate(xiou(nst1,nst3,n))
        allocate(xiuo(nst3,nst1,n))
        xiou(:,:,:)=emat(:,:,:,jknr)
        xiuo(:,:,:)=emata(:,:,:,jknr)
        deou(:,:)=deval(:,:,jknr)
        docc12(:,:)=docc(:,:,jknr)

        ! apply gauge wrt. symmetrized Coulomb potential
        call getpemat(iqmt,jknr,'PMAT_SCR.OUT','',m12=bufou,p12=pufou,m34=bufuo, &
         p34=pufuo)
        dekp(:,:)=deou(:,:)
        dokp(:,:)=docc12(:,:)
        sciskp(:,:)=scis(:,:,jknr)
        ! assign optical component
	do oct=1,noptc
	  emat12kp(:,:,-oct)=pufou(oct,:,:)
	  emat12kpa(:,:,-oct)=pufuo(oct,:,:)
	end do
        emat12kp(:,:,1:)=bufou(:,:,:)
        emat12kpa(:,:,1:)=bufuo(:,:,:)

        deallocate(xiou,xiuo)
        emattype=2
        call ematbdcmbs(emattype)

        ! get screened Coulomb interaction
        if (iknr.le.jknr) then
           call getbsemat('SCCLI.OUT',ikkp,nst1,nst3,sccli)
        else
           call getbsemat('SCCLI.OUT',ikkp,nst1,nst3,scclih)
           ! use Hermitian property for lower triangle
           do ist1=1,nst1
              do ist3=1,nst3
                 do ist2=1,nst1
                    do ist4=1,nst3
                       sccli(ist1,ist3,ist2,ist4)= &
                            conjg(scclih(ist2,ist4,ist1,ist3))
                    end do
                 end do
              end do
           end do
        end if
        ! proper sign of screened Coulomb interaction
        sccli=-sccli
        ! set diagonal of Bethe-Salpeter kernel to zero
        ! (cf. A. Marini, PRL 2003)
        if (iknr.eq.jknr) then
           do ist3=1,nst3
              do ist1=1,nst1
                 sccli(ist1,ist3,ist1,ist3)=zzero
              end do
           end do
        end if
        j1=0
        do ist2=1,nst3
           do ist1=1,nst1
              j1=j1+1
              emat12p(j1,:)=conjg(emat12kp(ist1,ist2,:))
              emat12pa(j1,:)=conjg(emat12kpa(ist2,ist1,:))
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
                    if (abs(t1).ge.fxcbsesplit) then
                       zmr(j2,j1)=zt1/t1
                       zmq(j2,j1)=zzero
                       zmra(j2,j1)=conjg(zt1)/t1
                       zmqa(j2,j1)=zzero
                    else
                       zmr(j2,j1)=zzero
                       zmq(j2,j1)=zt1
                       zmra(j2,j1)=zzero
                       zmqa(j2,j1)=conjg(zt1)
                    end if
		    
if ((iknr.eq.1).and.(jknr.eq.1).and.(ist1.eq.1).and.(ist3.eq.1).and.(ist2.eq.2).and. &
 (ist4.eq.1)) then
write(unitout,*) '1-1,2-1',j1,j2,sccli(ist1,ist3,ist2,ist4),zmr(j2,j1),dekp(2,1),dek(1,1),dekp(2,1)-dek(1,1)
end if			    
		    
                 end do
              end do
           end do
        end do
        ! calculate residual "R"; partial fraction decomposition without
	! double poles
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))
        residr=residr+matmul(zmr,emat12p)
        residra=residra+matmul(zmra,emat12pa)
        ! calculate residual "Q"; double poles part
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))
        residq=residq+matmul(zmq,emat12p)
        residqa=residqa+matmul(zmqa,emat12pa)
	
	
	
        ! end inner loop over k-points
     end do

!@@@@@@@@@@@@@@@@@@@@@@@@
write(unitout,*) 'iknr',iknr
write(unitout,*) '1-1',residr(1,-3:2)
write(unitout,*) '1-2',residr(2,-3:2)

     !--------------------------!
     !     set up BSE-kernel    !
     !--------------------------!
     t1=1.d0/(nkptnr*omega)
     do ist3=1,nst3
        do ist1=1,nst1
           osca(:,:)=zzero
           oscb(:,:)=zzero
           oscaa(:,:)=zzero
           oscba(:,:)=zzero
           j1=ist1+(ist3-1)*nst1
           ! set up inner part of kernel           
           ! generate oscillators
           call xszoutpr3(n+noptc+1,n+noptc+1,zone,emat12k(:,ist1,ist3), &
	   	residr(j1,:),osca)
          call xszoutpr3(n+noptc+1,n+noptc+1,zone,emat12ka(:,ist3,ist1), &
	   	residra(j1,:),oscaa)
		
           ! add Hermitian transpose 
	   forall(igq1=-3:n,igq2=-3:n)
	      osca(igq1,igq2)=osca(igq1,igq2)+conjg(osca(igq2,igq1))
	      oscaa(igq1,igq2)=oscaa(igq1,igq2)+conjg(oscaa(igq2,igq1))
	   end forall
	   
           call xszoutpr3(n+noptc+1,n+noptc+1,zone,emat12k(:,ist1,ist3), &
	   	residq(j1,:),oscb)
           call xszoutpr3(n+noptc+1,n+noptc+1,zone,emat12ka(:,ist3,ist1), &
	   	residqa(j1,:),oscba)
          ! set up energy denominators
           den1(:)=2.d0*t1/(w(:)+scisk(ist1,ist3)+dek(ist1,ist3)+zi*brd)
           den2(:)=2.d0*t1/(w(:)+scisk(ist1,ist3)+dek(ist1,ist3)+zi*brd)**2
           den1a(:)=2.d0*t1/(w(:)+scisk(ist1,ist3)-dek(ist1,ist3)+ &
	    	torfxc*zi*brd)
           den2a(:)=-2.d0*t1/(w(:)+scisk(ist1,ist3)-dek(ist1,ist3)+ &
	    	torfxc*zi*brd)**2
           ! update kernel
           do iw=1,nwdf
              ! resonant and antiresonant contributions
              fxc(:,:,iw)=fxc(:,:,iw)+osca(:,:)*den1(iw) !@@@@@@@@@+oscb(:,:)*den2(iw)
!@@@@@@              if (aresfxc) fxc(:,:,iw)=fxc(:,:,iw)+oscaa(:,:)*den1a(iw)+ &
!@@@@@@                   oscba(:,:)*den2a(iw)
           end do
           ! end loop over states #1
        end do
        ! end loop over states #3
     end do
     ! end outer loop over k-points
  end do

  ! filename for xc-kernel (ASCII)
  call genfilname(basename='FXC_BSE',asc=.true.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam2)
  call getunit(un)
  open(un,file=trim(filnam2),form='formatted',action='write',status='replace')

  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam3)
!  inquire(iolength=recl) (/(fxc(-oct,-oct,1),oct=1,noptc)/), &
!     	(/(fxc(-oct,1:,1),oct=1,noptc)/), &
!	(/(fxc(1:,-oct,1),oct=1,noptc)/),fxc(1:,1:,1)
  inquire(iolength=recl) fxc(-3:-1,-3:-1,1), fxc(-3:-1,1:,1), fxc(1:,-3:-1,1), &
	fxc(1:,1:,1)
  call getunit(un2)
  open(un2,file=trim(filnam3),form='unformatted',action='write', &
       status='replace',access='direct',recl=recl)

  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE_HEAD',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam4)
  call getunit(un3)
  open(un3,file=trim(filnam4),form='formatted',action='write',status='replace')
  
  do iw=1,nwdf
!     write(un2,rec=iw) (/(fxc(-oct,-oct,iw),oct=1,noptc)/), &
!     	(/(fxc(-oct,1:,iw),oct=1,noptc)/),(/(fxc(1:,-oct,iw),oct=1,noptc)/), &
!	fxc(1:,1:,iw)     
     write(un2,rec=iw) fxc(-3:-1,-3:-1,iw), fxc(-3:-1,1:,iw), fxc(1:,-3:-1,iw), &
	fxc(1:,1:,iw)     
     write(un3,'(i6,2x,g18.10,2x,6g18.10)') iw,dble(w(iw)),(fxc(-oct,-oct,iw),&
     	oct=1,noptc)
  end do
  do iw=1,nwdf,10
     do igq1=-noptc,n
        do igq2=-noptc,n
           write(un,'(3i6,3g18.10)') iw,igq1,igq2,fxc(igq1,igq2,iw), &
                abs(fxc(igq1,igq2,iw))
        end do
     end do
  end do  
  close(un)
  close(un2)
  close(un3)

  ! deallocate
  deallocate(den1,den2,den1a,den2a)
  deallocate(emat12p,zmr,zmq,dek,dekp,dde,dok,dokp,scisk,fxc)
  deallocate(sccli,scclih,scclit,emat12k,emat12kp,residr,residq,w,osca,oscb)
  deallocate(emat,deval,docc,scis)
  ! deallocate antiresonant parts
  deallocate(emata,emat12pa,emat12ka,emat12kpa,residra,residqa,zmra,zmqa)
  deallocate(oscaa,oscba)

  deallocate(bsedg)
  deallocate(bufou,bufuo,pufou,pufuo)
  
end subroutine kernxc_bse
!EOC
