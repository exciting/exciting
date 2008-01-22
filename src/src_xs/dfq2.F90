
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dfq2(iq)
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_genwgrid
  use m_gensymdf
  use m_getpemat2
  use m_getdevalsv2
  use m_dfqoschd
  use m_dfqoscwg
  use m_dfqoscbo
  use m_dftim
  use m_gettetcw
  use m_chi0upd
  use m_putx0
  use m_getunit
  use m_filedel
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='dfq2'
  real(8), parameter :: epstetra=1.d-8
  complex(8), allocatable :: w(:)
  complex(8), allocatable :: chi0(:,:,:),hou(:,:),huo(:,:)
  complex(8), allocatable :: chi0w(:,:,:,:),chi0h(:,:)
  complex(8), allocatable :: xou(:),xouc(:),xuo(:),xuoc(:),wou(:),wuo(:)
  complex(8) :: wout
  real(8), allocatable :: wreal(:),cw(:),cwa(:),cwsurf(:)
  real(8) :: brd,cpu0,cpu1,cpuread,cpuosc,cpuupd,cputot
  integer :: n,j,ik,iw,wi,wf,ist1,ist2,nwdfp,ikt,oct,un
  logical :: tq0
  logical, external :: tqgamma
  ! sampling of Brillouin zone
  bzsampl=0
  if (tetra) bzsampl=1
  ! initial and final w-point
  wi=wpari
  wf=wparf
  nwdfp=wf-wi+1
  ! matrix size for response function
  n=ngq(iq)
  ! zero broadening for analytic contiunation
  brd=brdtd
  if (acont) brd=zzero
  ! filenames for input
  call genfilname(basename='TETW',iq=iq,filnam=fnwtet)
  call genfilname(basename='EMAT',iq=iq,filnam=fnemat)
  call genfilname(basename='DEVALSV',iq=iq,filnam=fndevalsv)
  ! filenames for output
  call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,nar=.not.aresdf,&
       iq=iq,filnam=fnchi0)
  call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,nar=.not.aresdf,&
       iq=iq,procs=procs,rank=rank,filnam=fnchi0_t)
  call genfilname(nodotpar=.true.,basename='X0_TIMING',iq=iq,&
       procs=procs,rank=rank,filnam=fnxtim)
  ! file extension for q-point
  call genfilname(iq=iq,setfilext=.true.)
  ! remove timing files from previous runs
  call filedel(trim(fnxtim))
  ! calculate k+q and G+k+q related variables
  call init1xs(qvkloff(1,iq))
  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  ! find limits for band combinations
  call ematbdcmbs(emattype)
  ! check if q-point is Gamma point
  tq0=tqgamma(iq)
  if (tq0) then
     write(unitout,'(a)') 'Info('//trim(thisnam)//'): Gamma q-point: using &
          &momentum matrix elements for dielectric function'
  end if
  ! write out matrix size of response function
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors &
       &(local field effects):',ngq(iq)
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): lowest (partially) &
       & unoccupied state: ',istunocc0
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): highest (partially)&
       & occupied state  : ',istocc0
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): limits for band combinations&
       & nst1,nst2,nst3,nst4:',nst1,nst2,nst3,nst4
  ! allocate arrays for eigenvalue and occupation number differences
  if (allocated(deou)) deallocate(deou)
  allocate(deou(nst1,nst2))
  if(allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst3,nst4))
  if (allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1,nst2))
  if (allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst3,nst4))
  ! allocate matrix elements arrays
  if (allocated(xiou)) deallocate(xiou)
  allocate(xiou(nst1,nst2,n))
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiuo(nst3,nst4,n))
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3,nst1,nst2))
  if (allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3,nst3,nst4))
  ! allocate arrays
  allocate(w(nwdf))
  allocate(wreal(nwdfp))
  allocate(chi0h(3,nwdfp))
  allocate(chi0w(n,2,3,nwdfp))
  allocate(chi0(n,n,nwdfp))
  allocate(wou(nwdf))
  allocate(wuo(nwdf))
  allocate(xou(n))
  allocate(xouc(n))
  allocate(xuo(n))
  allocate(xuoc(n))
  allocate(hou(n,n))
  allocate(huo(n,n))
  if (tetra) allocate(cw(nwdf),cwa(nwdf),cwsurf(nwdf))
  ! generate complex energy grid
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
  wreal(:)=w(wi:wf)
  if (wreal(1).lt.epstetra) wreal(1)=epstetra
  ! initializations
  chi0(:,:,:)=zzero
  chi0w(:,:,:,:)=zzero
  chi0h(:,:)=zzero
  ! loop over k-points
  call getunit(un)
  ikt=0
  do ik=1,nkpt
     cpuosc=0.d0
     cpuupd=0.d0
     call cpu_time(cpu0)
     ! read Kohn-Sham energy differences
     call getdevalsv2(iq,ik,.true.,trim(fndevalsv),deou,docc12,deuo,docc21)
     ! read Kohn-Sham energy differences (random k-point set)
     ! get matrix elements (exp. expr. or momentum)
     call getpemat2(iq,ik,trim(fnpmat),trim(fnemat),m12=xiou,m34=xiuo, &
          p12=pmou,p34=pmuo)
     ! turn off antiresonant terms (type 21 band combiantions) for Kohn-Sham
     ! response function
     if (.not.aresdf) then
        xiuo(:,:,:)=zzero
        pmuo(:,:,:)=zzero
     end if
     ! avoid double counting *** zero lower/upper triangle
     ! of xiou,xiou because of scissors shift *** scissors: sign relevant:
     ! calculate explicitly the scissors shift, depending on sign of band
     ! energy difference
     if (istunocc0.le.istocc0) then
        xiuo(:istocc0-istunocc0+1,istunocc0:,:)=zzero
        pmuo(:,istocc0-istunocc0+1,istunocc0:)=zzero
!!$        xiou(istunocc0:,:istocc0-istunocc0+1,:)=zzero
!!$        pmou(:,istunocc0:,:istocc0-istunocc0+1)=zzero
        ! take into account intraband contributions
        if (tq0.or.(.not.intraband)) then
           do ist1=1,nst1
              do ist2=1,nst2
                 if (ist1.eq.istunocc0+ist2-1) then
                    xiou(ist1,ist2,:)=zzero
                    pmou(:,ist1,ist2)=zzero
                 end if
                 if (istunocc0+ist1-1.eq.ist2) then
                    xiuo(ist2,ist1,:)=zzero
                    pmuo(:,ist2,ist1)=zzero
                 end if
              end do
           end do
        end if
     end if
     
     call cpu_time(cpu1)
     cpuread=cpu1-cpu0
     do ist1=1,nst1
        do ist2=1,nst2
           !---------------------!
           !     denominator     !
           !---------------------!
           call cpu_time(cpu0)
           ! user request termination
           call terminate_inqr('dfq')
           if (tetra) then
              ! read weights for tetrahedron method
              call gettetcw(iq,ik,ist1,ist2,nwdf,trim(fnwtet),cw,cwa, &
                   cwsurf)
              ! include occupation number differences
              wou(wi:wf)=docc12(ist1,ist2)*cmplx(cw(wi:wf),cwsurf(wi:wf),8)/ &
                   omega
              wuo(wi:wf)=-docc21(ist2,ist1)*cmplx(cwa(wi:wf),0.d0,8)/omega
           else
              ! include occupation number differences
              wou(:)=docc12(ist1,ist2)*wkpt(ik)/omega/(w(:)+deou(ist1,ist2)- &
                   scissor+zi*brd)
              wuo(:)=docc21(ist2,ist1)*wkpt(ik)/omega/(w(:)+deuo(ist2,ist1)+ &
                   scissor+zi*brd)
           end if
           hou(:,:)=zzero
           huo(:,:)=zzero
           !---------------------!
           !     oscillators     !
           !---------------------!
           ! calculate oscillators
           if (.not.tq0) then
              ! set up body, head and wings in one
              call dfqoscbo(n,xiou(ist1,ist2,:),xiuo(ist2,ist1,:),hou,huo)
           end if

if (ik.eq.7) then
   write(3000,*) ik,ist1,ist2,hou(1,1)/fourpi 
   write(4000,*) ik,ist1,ist2,deou(ist1,ist2)
   write(5000,*) ik,ist1,ist2,docc12(ist1,ist2)
end if

           if (tq0.and.(n.gt.1)) then
              ! set up body
              call dfqoscbo(n-1,xiou(ist1,ist2,2:),xiuo(ist2,ist1,2:), &
                   hou(2:,2:),huo(2:,2:))
           end if
           ! loop over longitudinal Cartesian (diagonal) components of
           ! response function
           do oct=1,3
              ! symmetrization matrix for dielectric function
              call gensymdf(oct,oct)
              optcomp(1,1)=oct
              optcomp(2,1)=oct
              ! Gamma q-point
              if (tq0) then
                 ! head
                 call dfqoschd(pmou(:,ist1,ist2),pmuo(:,ist2,ist1),hou(1,1), &
                      huo(1,1))
                 do iw=wi,wf
                    wout=wou(iw)
                    ! be careful with gauge in the w-variable
                    ! one has to subtract the scissor's shift
                    if (tetra) wout=cmplx(dble(wou(iw)),aimag(wou(iw))*&
                         deou(ist1,ist2)**2/(wreal(iw-wi+1)-scissor)**2)
                    chi0h(oct,iw-wi+1)=chi0h(oct,iw-wi+1)+ &
                         wout*hou(1,1)+wuo(iw)*huo(1,1)
                 end do
              end if
              ! Gamma q-point
              if (tq0.and.(n.gt.1)) then
                 ! wings
                 call dfqoscwg(1,pmou(:,ist1,ist2),pmuo(:,ist2,ist1), &
                      xiou(ist1,ist2,2:),xiuo(ist2,ist1,2:),hou(1,2:), &
                      huo(1,2:))
                 call dfqoscwg(2,pmou(:,ist1,ist2),pmuo(:,ist2,ist1), &
                      xiou(ist1,ist2,2:),xiuo(ist2,ist1,2:),hou(2:,1), &
                      huo(2:,1))
                 do iw=wi,wf
                    wout=wou(iw)
                    ! be careful with gauge in the w-variable
                    ! one has to subtract the scissor's shift
                    if (tetra) wout=cmplx(dble(wou(iw)),aimag(wou(iw))*&
                         deou(ist1,ist2)/(-wreal(iw-wi+1)+scissor))
                    chi0w(2:,1,oct,iw-wi+1)=chi0w(2:,1,oct,iw-wi+1)+&
                         wout*hou(1,2:)+wuo(iw)*huo(1,2:)
                    chi0w(2:,2,oct,iw-wi+1)=chi0w(2:,2,oct,iw-wi+1)+&
                         wout*hou(2:,1)+wuo(iw)*huo(2:,1)
                 end do
              end if
              call cpu_time(cpu1)
              cpuosc=cpuosc+cpu1-cpu0
           end do !oct
           !----------------------------------!
           !     update response function     !
           !----------------------------------!
           do iw=wi,wf
              ! * most time-consuming part of routine *
              call chi0upd(n,wou(iw),wuo(iw),hou,huo,&
                   chi0(:,:,iw-wi+1))
           end do
           call cpu_time(cpu0)
           cpuupd=cpuupd+cpu0-cpu1
           ! end loop over states combinations
        end do
     end do
     cputot=cpuread+cpuosc+cpuupd
     ! timing information
     call dftim(iq,ik,trim(fnxtim),cpuread,cpuosc,cpuupd, &
          cputot)
     ! synchronize
     call barrier
     ! end loop over k-points
  end do
  ! write response function to file
  do j=0,procs-1
     if (rank.eq.j) then
        do iw=wi,wf
           call putx0(tq0,iq,iw-wi+1,trim(fnchi0_t),'',&
                chi0(:,:,iw-wi+1),chi0w(:,:,:,iw-wi+1),chi0h(:,iw-wi+1))
        end do
     end if
     call barrier
  end do
  deallocate(chi0h)
  deallocate(chi0w)
  deallocate(deou,deuo,wou,wuo)
  deallocate(xiou,xiuo,pmou,pmuo)
  deallocate(w,wreal,chi0)
  deallocate(xou,xouc,xuo,xuoc,hou,huo)
  if (tetra) deallocate(cw,cwa,cwsurf)
end subroutine dfq2
