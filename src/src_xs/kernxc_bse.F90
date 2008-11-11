
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
  !******************************************************************
  integer, parameter :: oct=1
  !******************************************************************
  character(*), parameter :: thisnam='kernxs_bse'
  integer, parameter :: iqmt=1
  real(8), parameter :: delt=1.d-6
  real(8), parameter :: epsortho=1.d-12
  character(256) :: fname,filnam,filnam2,filnam3,filnam4
  logical :: tq0
  integer :: iv(3),iw,wi,wf,nwdfp,n,recl,un,un2,un3,j1,j2
  integer :: ikkp,iknr,jknr,iknrq,jknrq,iqr,iq,iqrnr,igq1,igq2
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24
  integer :: ikkp_,iknr_,jknr_,iq_,iqr_
  integer :: nst1_,nst2_,nst3_,nst4_
  real(8) :: vqr(3),vq(3),t1,brd
  real(8) :: cpu_init1offs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1
  real(8) :: cpu_ematqdealloc,cpu_clph,cpu_suma,cpu_write
  complex(8) :: zt1
  real(8), allocatable :: dek(:,:),dok(:,:),scisk(:,:)
  real(8), allocatable :: dekp(:,:),dokp(:,:),sciskp(:,:)
  real(8), allocatable :: deval(:,:,:),docc(:,:,:),scis(:,:,:)
  real(8), allocatable :: dde(:,:)
  complex(8), allocatable :: zmr(:,:),zmq(:,:)
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:),scclih(:,:,:,:)
  complex(8), allocatable :: emat(:,:,:,:),den1(:),den2(:)
  complex(8), allocatable :: emat12(:,:),emat12p(:,:)
  complex(8), allocatable :: emat12k(:,:,:),emat12kp(:,:,:)
  complex(8), allocatable :: residr(:,:),residq(:,:),osca(:,:),oscb(:,:)
  complex(8), allocatable :: fxc(:,:,:),w(:)
  ! external functions
  integer, external :: idxkkp,l2int
  logical, external :: tqgamma
  brd=broad
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
  fxc(:,:,:)=zzero
  sccli(:,:,:,:)=zzero
  allocate(emat(nst1,nst3,n,nkptnr))
  allocate(deval(nst1,nst3,nkptnr))
  allocate(docc(nst1,nst3,nkptnr))
  allocate(scis(nst1,nst3,nkptnr))

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
          & of plane wave')    
     iknrq=ikmapikq(iknr,iqmt)
     ! matrix elements for k and q=0
     call ematqk1(iqmt,iknr)
     emat(:,:,:,iknr)=xiou(:,:,:)
     deallocate(xiou,xiuo)
     call getdevaldoccsv(iqmt,iknr,iknrq,istl1,istu1,istl2,istu2,deou, &
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
     call chkpt(3,(/task,3,ikkp/),'task,sub,(k,kp)-pair; BSE-fxc-kernel')
     iknrq=ikmapikq(iknr,iqmt)

     emattype=1
     call ematbdcmbs(emattype)
     allocate(xiou(nst1,nst3,n))
     xiou(:,:,:)=emat(:,:,:,iknr)
     deou(:,:)=deval(:,:,iknr)
     docc12(:,:)=docc(:,:,iknr)

     ! apply gauge wrt. symmetrized Coulomb potential
     call getpemat(iqmt,iknr,'PMAT_SCR.OUT','',m12=xiou,p12=pmou)
     dek(:,:)=deou(:,:)
     dok(:,:)=docc12(:,:)
     ! add BSE diagonal
     scisk(:,:)=scis(:,:,iknr)+bsed
     ! assign optical component
     xiou(:,:,1)=pmou(oct,:,:)
     do igq1=1,n
        emat12k(igq1,:,:)=xiou(:,:,igq1)
     end do

     deallocate(xiou)
     emattype=2
     call ematbdcmbs(emattype)

     residr(:,:)=zzero
     residq(:,:)=zzero
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
        iqrnr=iqmap(ivqr(1,iqr),ivqr(2,iqr),ivqr(3,iqr))

        emattype=1
        call ematbdcmbs(emattype)
        allocate(xiou(nst1,nst3,n))
        xiou(:,:,:)=emat(:,:,:,jknr)
        deou(:,:)=deval(:,:,jknr)
        docc12(:,:)=docc(:,:,jknr)

        ! apply gauge wrt. symmetrized Coulomb potential
        call getpemat(iqmt,jknr,'PMAT_SCR.OUT','',m12=xiou,p12=pmou)
        dekp(:,:)=deou(:,:)
        dokp(:,:)=docc12(:,:)
        sciskp(:,:)=scis(:,:,jknr)
        ! assign optical component
        xiou(:,:,1)=pmou(oct,:,:)
        emat12kp(:,:,:)=xiou(:,:,:)

        deallocate(xiou)
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
!!!                    zt1=sccli(ist2,ist4,ist1,ist3)
                    ! four point energy difference
                    t1=dekp(ist2,ist4)-dek(ist1,ist3)
                    ! arrays for R- and Q-residuals
                    if (abs(t1).ge.delt) then
                       zmr(j2,j1)=zt1/t1
                       zmq(j2,j1)=zzero
                    else
                       zmr(j2,j1)=zzero
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
           call xszoutpr3(n,n,zone,emat12k(:,ist1,ist3),residr(j1,:),osca)
           ! add Hermitian transpose
           osca=osca+conjg(transpose(osca))
           call xszoutpr3(n,n,zone,emat12k(:,ist1,ist3),residq(j1,:),oscb)
           ! *** this part is working for Si_lapw and Si_APW+lo ***
           ! set up energy denominators
           den1(:)=2.d0/(w(:)+scisk(ist1,ist3)+dek(ist1,ist3)+zi*brd)
           den2(:)=2.d0/(w(:)+scisk(ist1,ist3)+dek(ist1,ist3)+zi*brd)**2
           den1=den1/nkpt/omega
           den2=den2/nkpt/omega
           ! *** end
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

  ! filename for xc-kernel (ASCII)
  call genfilname(basename='FXC_BSE',asc=.true.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam2)
  call getunit(un)
  open(un,file=trim(filnam2),form='formatted',action='write',status='replace')

  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam3)
  inquire(iolength=recl) fxc(:,:,1)
  call getunit(un2)
  open(un2,file=trim(filnam3),form='unformatted',action='write', &
       status='replace',access='direct',recl=recl)

  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE_HEAD',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam4)
  call getunit(un3)
  open(un3,file=trim(filnam4),form='formatted',action='write',status='replace')
  
  do iw=1,nwdf
     write(un2,rec=iw) fxc(:,:,iw)     
     write(un3,'(i6,2x,g18.10,2x,2g18.10)') iw,dble(w(iw)),fxc(1,1,iw)
  end do
  do iw=1,nwdf,10
     do igq1=1,n
        do igq2=1,n
           write(un,'(3i6,3g18.10)') iw,igq1,igq2,fxc(igq1,igq2,iw), &
                abs(fxc(igq1,igq2,iw))
	end do
     end do
  end do  
  close(un)
  close(un2)
  close(un3)

  ! deallocate
  deallocate(emat12,emat12p,zmr,zmq,dek,dekp,dde,dok,dokp,scisk,sciskp,fxc)
  deallocate(sccli,scclih,scclit,emat12k,emat12kp,residr,residq,w,osca,oscb)
  deallocate(emat,deval,docc,scis)
  
end subroutine kernxc_bse
!EOC
