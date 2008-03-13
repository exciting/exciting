
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine kernxc_bse
  !
  ! BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
  !
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
  use m_getx0
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(256) :: filnam,filnam2
  complex(8),allocatable :: chi0(:,:), fxc(:,:), idf(:,:), mdf1(:),w(:)
  complex(8),allocatable :: chi0hd(:),chi0wg(:,:,:),chi0h(:)
  integer :: n,m,recl,j,iw,wi,wf,nwdfp,nc,oct,oct1,oct2,igmt
  integer, external :: l2int
  ! local variables
  character(*), parameter :: thisnam = 'kernxs_bse'
  complex(8) :: zt1
  integer :: sh(2),ig
  ! ***    
  character(256) :: fname
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iqrnr,isym,jsym,jsymi,igq1,igq2,iflg
  integer :: ngridkt(3),iv(3),ivgsym(3),un,j1,j2
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp
  logical :: nosymt,reducekt,tq0,nsc,tphf
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
  real(8) :: cpu_init1xs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1
  real(8) :: cpu_ematqdealloc,cpu_clph,cpu_suma,cpu_write
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


  call genparidxran('w')
  ! sampling type for Brillouin zone sampling
  bzsampl=l2int(tetra)
  ! limits for w-points
  wi=wpari
  wf=wparf
  nwdfp=wparf-wpari+1
  ! matrix size for local field effects (first q-point is Gamma-point)
  n=ngq(1)
  allocate(fxc(n,n))
  fxc(:,:)=zzero


  allocate(sccli(nst1,nst3,nst2,nst4))
  allocate(emat12k(nst1,nst3,ngq(1)),emat12kp(nst1,nst3,ngq(1)))
  sccli(:,:,:,:)=zzero

  call genfilname(basename='SCCLI',dotext='.OUT',filnam=fname)
  call getunit(un)
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4, &
       sccli(:,:,:,:)
  open(un,file=trim(fname),form='unformatted',action='read', &
       status='old',access='direct',recl=recl)

  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  ikkp=0
  ! first k-point
  do iknr=1,nkptnr

     ! matrix elements for k and q=0
     emattype=1
     call init1xs(qvkloff(1,1))
     call getematrad(1)
     call ematqalloc
     call cpu_time(cpu0)
     call ematqk1(1,iknr)
     call cpu_time(cpu1)
     cpu_ematqk1=cpu_ematqk1+cpu1-cpu0
     emat12k(:,:,:)=xiou(:,:,:)
     deallocate(xiou,xiuo)
     emattype=2
     call ematbdcmbs(emattype)

     ! get eigenvalue differences for k-point
     ! ******


     ! second k-point
     do jknr=1,nkptnr

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

        write(*,'(a,i6,2x,2i5,2x,2i5)') 'ikkp,iknr,jknr,iq,iqr',&
             ikkp,iknr,jknr,iq,iqr

        ! matrix elements for kp and q=0
        emattype=1
        call init1xs(qvkloff(1,1))
        call getematrad(1)
        call ematqalloc
        call cpu_time(cpu0)
        call ematqk1(1,jknr)
        call cpu_time(cpu1)
        cpu_ematqk1=cpu_ematqk1+cpu1-cpu0
        emat12kp(:,:,:)=xiou(:,:,:)
        deallocate(xiou,xiuo)
        emattype=2
        call ematbdcmbs(emattype)


!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        do igq1=1,n
!!$           do ist1=1,nst1
!!$              do ist2=1,nst3
!!$                 write(4000,'(4i5,3g18.10)') iknr,igq1,ist1,ist2, &
!!$                      emat12k(ist1,ist2,igq1),abs(emat12k(ist1,ist2,igq1))**2
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

        ! get eigenvalue differences for kp-point
        ! ******

        ! get screened Coulomb interaction

        ! assign lower triangular part
        if (jknr.gt.iknr) then
           ! read from file

        else
           ! use Hermitian property

        end if


        ! calculate residual "R" 
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))

        ! calculate residual "Q" 
        ! (cf. A. Marini, Phys. Rev. Lett. 91, 256402 (2003))


        ! set up inner part of kernel

        ! multiply with inverse KS-response function from both sides
        ! with energies shifted by BSE-diagonal


        call cpu_time(cpu3)
        t3=cpu_ematqdealloc+cpu_ematqk1+cpu_ematqalloc+cpu_ematrad+cpu_init1xs+cpu_clph+cpu_suma+cpu_write
        write(*,'(a,f12.3)') 'init1xs     :',cpu_init1xs
        write(*,'(a,f12.3)') 'ematrad     :',cpu_ematrad
        write(*,'(a,f12.3)') 'ematqalloc  :',cpu_ematqalloc
        write(*,'(a,f12.3)') 'ematqk1     :',cpu_ematqk1
        write(*,'(a,f12.3)') 'ematqdealloc:',cpu_ematqdealloc
        write(*,'(a,f12.3)') 'summation   :',cpu_suma
        write(*,'(a,f12.3)') 'write       :',cpu_write
        write(*,'(a,f12.3)') '*** sum     :',t3
        write(*,'(a,f12.3)') '*** rest    :',cpu3-cpu2-t3
        write(*,'(a,f12.3)') '*** overall :',cpu3-cpu2
        write(*,*)
        ! end loop over (k,kp) pairs
     end do

     if (iknr.eq.1) stop 'stop in kernxc_bse'

  end do



  ! deallocate
  deallocate(fxc)
end subroutine kernxc_bse
