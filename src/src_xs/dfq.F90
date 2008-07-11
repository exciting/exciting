
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dfq(iq)
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_genwgrid
  use m_gensymdf
  use m_getpemat
  use m_dfqoschd
  use m_dfqoscwg
  use m_dfqoscbo
  use m_dftim
  use m_gettetcw
  use m_chi0upd
  use m_putx0
  use m_getunit
  use m_writevars
  use m_filedel
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='dfq'
  character(256) :: fnscreen
  real(8), parameter :: epstetra=1.d-8
  complex(8), allocatable :: w(:)
  complex(8), allocatable :: chi0(:,:,:),hou(:,:),huo(:,:),hdg(:,:,:)
  complex(8), allocatable :: chi0w(:,:,:,:),chi0h(:,:)
  complex(8), allocatable :: xou(:),xouc(:),xuo(:),xuoc(:),wou(:),wuo(:)
  complex(8) :: wout
  real(8), allocatable :: wreal(:),cw(:),cwa(:),cwsurf(:)
  real(8), allocatable :: scis12(:,:),scis21(:,:)
  real(8) :: brd,cpu0,cpu1,cpuread,cpuosc,cpuupd,cputot,rv1(9),r1
  integer :: n,j,i1,i2,j1,j2,ik,ikq,igq,iw,wi,wf,ist1,ist2,nwdfp
  integer :: oct,oct1,oct2,un,ig1,ig2
  logical :: tq0
  integer, external :: octmap
  logical, external :: tqgamma
  if (acont.and.tscreen) then
     write(*,*)
     write(*,'("Error(",a,"): analytic continuation does not work for &
          &screening")')
     write(*,*)
     call terminate
  end if
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
  if (acont) brd=0.d0
  ! *** experimental *** zero broadening for dielectric matrix (w=0)
  ! for band-gap systems
  if (task.eq.430) brd=0.d0
  ! file extension for q-point
  if (.not.tscreen) call genfilname(iqmt=iq,setfilext=.true.)
  ! filenames for input
  ! filenames for output
  if (tscreen) then
     call genfilname(basename='TETW',iq=iq,appfilext=.true.,filnam=fnwtet)
     call genfilname(basename='PMAT',appfilext=.true.,filnam=fnpmat)
     call genfilname(basename='SCREEN',bzsampl=bzsampl,nar=.not.aresdf,&
          iq=iq,filnam=fnscreen)
     call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
          etype=emattype,procs=procs,rank=rank,appfilext=.true.,filnam=fnetim)
     call genfilname(nodotpar=.true.,basename='X0_TIMING',iq=iq,&
          bzsampl=bzsampl,acont=acont,nar=.not.aresdf,procs=procs,rank=rank, &
          appfilext=.true.,filnam=fnxtim)
  else
     call genfilname(basename='TETW',iqmt=iq,filnam=fnwtet)
     call genfilname(basename='PMAT_XS',filnam=fnpmat)
     call genfilname(basename='EMAT',iqmt=iq,filnam=fnemat)
     call genfilname(nodotpar=.true.,basename='X0_TIMING',bzsampl=bzsampl,&
          acont=acont,nar=.not.aresdf,iqmt=iq,procs=procs,rank=rank, &
          filnam=fnxtim)
     call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,nar=.not.aresdf,&
          iqmt=iq,filnam=fnchi0)
     call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,nar=.not.aresdf,&
          iqmt=iq,procs=procs,rank=rank,filnam=fnchi0_t)
  end if
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
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): limits for band combinations&
       & istlo1,isthi1,istlo2,isthi2:',istlo1,isthi1,istlo2,isthi2
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): limits for band combinations&
       & istlo3,isthi3,istlo4,isthi4:',istlo3,isthi3,istlo4,isthi4
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
  
  allocate(hdg(nst1,nst2,nkpt))
  
  allocate(scis12(nst1,nst2))
  allocate(scis21(nst2,nst1))
  allocate(w(nwdf))
  allocate(wreal(nwdfp))
  allocate(chi0h(9,nwdfp))
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
  scis12(:,:)=0.d0
  scis21(:,:)=0.d0
  if (tetra) allocate(cw(nwdf),cwa(nwdf),cwsurf(nwdf))
  ! generate complex energy grid
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
  wreal(:)=w(wi:wf)
  if (wreal(1).lt.epstetra) wreal(1)=epstetra
  ! initializations
  chi0(:,:,:)=zzero
  chi0w(:,:,:,:)=zzero
  chi0h(:,:)=zzero
  if (tscreen) then
     ! generate radial integrals wrt. sph. Bessel functions
     call ematrad(iq)
     ! delete timing information of previous runs
     call filedel(trim(fnetim))
     ! write information
     write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors:', &
          ngq(iq)
     call ematqalloc
  end if

!*******************************************************************************
hdg=zzero
write(*,*) 'dfq, shape(hdg)',shape(hdg)
!	read(1108) hdg
!*******************************************************************************

  ! loop over k-points
  do ik=1,nkpt
     write(*,'(a,i5,3x,2i6)') 'dfq: q-point/k-point/k+q-point:',iq,ik, &
          ikmapikq(ik,iq)
     cpuosc=0.d0
     cpuupd=0.d0
     call cpu_time(cpu0)
     ikq=ikmapikq(ik,iq)
     call getdevaldoccsv(iq,ik,ikq,istlo1,isthi1,istlo2,isthi2,deou,docc12, &
          scis12)
     call getdevaldoccsv(iq,ik,ikq,istlo2,isthi2,istlo1,isthi1,deuo,docc21, &
          scis21)
     if (tscreen) then
        ! for screening calculate matrix elements of plane wave on the fly
        call ematqk1(iq,ik)
        if (.not.allocated(xiuo)) allocate(xiuo(nst3,nst4,n))
        if (.not.allocated(pmuo)) allocate(pmuo(3,nst3,nst4))
     end if

!*******************************************************************************
! *** this is working for Si_lapw/apw+lo
	scis12=scis12+hdg(:,:,ik)
	scis21=scis21-hdg(:,:,ik)
!*******************************************************************************

     ! get matrix elements (exp. expr. or momentum op.)
     call getpemat(iq,ik,trim(fnpmat),trim(fnemat),m12=xiou,m34=xiuo, &
          p12=pmou,p34=pmuo)
     if (tscreen) then
        ! we don't need anti-resonant parts here, assign them the same
        ! value as for resonant parts, resulting in a factor of two.
        do igq=1,n
           xiuo(:,:,igq)=transpose(xiou(:,:,igq))
        end do
        do j=1,3
           pmuo(j,:,:)=transpose(pmou(j,:,:))
        end do
        deuo(:,:)=transpose(deou(:,:))
        docc21(:,:)=transpose(docc12(:,:))
        scis21(:,:)=transpose(scis12(:,:))
     end if
     ! turn off antiresonant terms (type 2-1 band combiantions) for Kohn-Sham
     ! response function
     if (.not.aresdf) then
        xiuo(:,:,:)=zzero
        pmuo(:,:,:)=zzero
     end if
     do ist1=1,istocc0-istunocc0+1
        do ist2=1,istocc0-istunocc0+1
           j=ist1+istunocc0-1
           ! set lower triangle of first block to zero
           if (ist1.gt.ist2) then
              xiou(j,ist2,:)=zzero
              pmou(:,j,ist2)=zzero
           end if
           ! set diagonal to zero (project out intraband contributions)
           if ((.not.intraband).and.(ist1.eq.ist2)) then
              xiou(j,ist2,:)=zzero
              pmou(:,j,ist2)=zzero
           end if
           ! set upper triangle of second block to zero
           ! also set diagonal to zero to avoid double counting
           if (ist1.ge.ist2) then
              xiuo(ist2,j,:)=zzero
              pmuo(:,ist2,j)=zzero
           end if
        end do
     end do
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
              ! absolute band indices
              i1=ist1; i2=istunocc0+ist2-1
              ! mirror index pair on diagonal if necessary
              if (i1.gt.i2) then
                 j1=ist2
                 j2=ist1-istunocc0+1
              else
                 j1=ist1
                 j2=ist2
              end if
              ! read weights for tetrahedron method
              call gettetcw(iq,ik,j1,j2,nst1,nst2,nwdf,trim(fnwtet),cw,cwa, &
                   cwsurf)
              ! include occupation number differences
              wou(wi:wf)=docc12(ist1,ist2)*cmplx(cw(wi:wf),cwsurf(wi:wf),8)/ &
                   omega
              wuo(wi:wf)=-docc21(ist2,ist1)*cmplx(cwa(wi:wf),0.d0,8)/omega
           else
              ! include occupation number differences
              wou(:)=docc12(ist1,ist2)*wkpt(ik)/omega/(w(:)+deou(ist1,ist2) &
                   +scis12(ist1,ist2)+zi*brd)
              wuo(:)=docc21(ist2,ist1)*wkpt(ik)/omega/(w(:)+deuo(ist2,ist1) &
                   +scis21(ist2,ist1)+zi*brd)
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
           if (tq0.and.(n.gt.1)) then
              ! set up body
              call dfqoscbo(n-1,xiou(ist1,ist2,2:),xiuo(ist2,ist1,2:), &
                   hou(2:,2:),huo(2:,2:))
           end if
           ! loop over longitudinal Cartesian (diagonal) components of
           ! response function
           if (tq0) then
              do oct1=1,3
                 optcomp(1,1)=oct1
                 optcomp(2,1)=oct1
                 if (n.gt.1) then
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
                            deou(ist1,ist2)/(-wreal(iw-wi+1)-scis12(ist1,ist2)))
                       chi0w(2:,1,oct1,iw-wi+1)=chi0w(2:,1,oct1,iw-wi+1)+&
                            wout*hou(1,2:)+wuo(iw)*huo(1,2:)
                       chi0w(2:,2,oct1,iw-wi+1)=chi0w(2:,2,oct1,iw-wi+1)+&
                            wout*hou(2:,1)+wuo(iw)*huo(2:,1)
                    end do
                 end if
                 do oct2=1,3
                    oct=octmap(oct1,oct2)
                    ! symmetrization matrix for dielectric function
                    call gensymdf(oct1,oct2)
                    optcomp(1,1)=oct1
                    optcomp(2,1)=oct2
                    ! head
                    call dfqoschd(pmou(:,ist1,ist2),pmuo(:,ist2,ist1),hou(1,1),&
                         huo(1,1))
                    do iw=wi,wf
                       wout=wou(iw)
                       ! be careful with gauge in the w-variable
                       ! one has to subtract the scissor's shift
                       if (tetra) wout=cmplx(dble(wou(iw)),aimag(wou(iw))*&
                            deou(ist1,ist2)**2 &
                            /(wreal(iw-wi+1)+scis12(ist1,ist2))**2) !SAG
                       chi0h(oct,iw-wi+1)=chi0h(oct,iw-wi+1)+ &
                            wout*hou(1,1)+wuo(iw)*huo(1,1)
                    end do
                 end do !oct2
              end do !oct1
           end if
           call cpu_time(cpu1)
           cpuosc=cpuosc+cpu1-cpu0	   
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
     if (.not.tscreen) call barrier
     ! end loop over k-points
  end do
  if (tscreen) call ematqdealloc
  ! write response function to file
  if (tscreen) then
     ! write out screening
     call getunit(un)
     open(un,file=trim(fnscreen),form='formatted',action='write', &
          status='replace')
     rv1(:)=0.d0
     rv1(1::4)=1.d0
     do ig1=1,n
        do ig2=1,n
           r1=0.d0
           if (ig1.eq.ig2) r1=1.d0
           if (tq0) then
              if ((ig1.eq.1).and.(ig2.eq.1)) then
                 write(un,'(3i8,2g18.10)') (ig1,ig2,j,rv1(j)-chi0h(j,1),j=1,9)
              end if
              if ((ig1.eq.1).and.(ig2.ne.1)) then
                 write(un,'(3i8,2g18.10)') (ig1,ig2,j,-chi0w(ig2,1,j,1),j=1,3)
              end if
              if ((ig1.ne.1).and.(ig2.eq.1)) then
                 write(un,'(3i8,2g18.10)') (ig1,ig2,j,-chi0w(ig1,2,j,1),j=1,3)
              end if
              if ((ig1.ne.1).and.(ig2.ne.1)) write(un,'(3i8,2g18.10)') ig1,ig2,&
                   0,r1-chi0(ig1,ig2,1)
           else
              write(un,'(3i8,2g18.10)') ig1,ig2,0,r1-chi0(ig1,ig2,1)
           end if
        end do
     end do
     call writevars(un,iq)
     close(un)
  else
     do j=0,procs-1
        if (rank.eq.j) then
           do iw=wi,wf
              call putx0(tq0,iq,iw-wi+1,trim(fnchi0_t),'',&
                   chi0(:,:,iw-wi+1),chi0w(:,:,:,iw-wi+1),chi0h(:,iw-wi+1))
           end do
        end if
        call barrier
     end do
  end if

  deallocate(docc12,docc21,scis12,scis21)
  deallocate(chi0h)
  deallocate(chi0w)
  deallocate(deou,deuo,wou,wuo)
  deallocate(xiou,xiuo,pmou,pmuo)
  deallocate(w,wreal,chi0)
  deallocate(xou,xouc,xuo,xuoc,hou,huo)
  if (tetra) deallocate(cw,cwa,cwsurf)
end subroutine dfq
