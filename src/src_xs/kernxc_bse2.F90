
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: kernxc_bse2
! !INTERFACE:
subroutine kernxc_bse2
! !USES:
  use modmain
  use modxs
  use m_getevalsvr
  use m_getpemat
  use m_xsgauntgen
  use m_findgntn0
  use m_genwgrid
  use m_xszoutpr3
  use m_genfilname
  use m_getunit
  use m_xszoutpr3
! !INPUT/OUTPUT PARAMETERS:
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
!
! !REVISION HISTORY:
!   Created February 2009 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='kernxc_bse2'
  real(8), parameter :: eps=1.d-5
  integer, parameter :: iqmt=1,noptc=3,iop=3,ig0=1
  character(256) :: filnam
  integer :: ikkp,iknr,jknr,iv,ic,jv,jc,si,sj,nv,nc,wsiz,n,i,j,iw,un,recl
  real(8) :: t1,dei,dej
  complex(8) :: zt1
  complex(8), allocatable :: w(:),wmat(:,:),wmatq(:,:),wm(:,:,:,:)
  complex(8), allocatable :: resr(:,:),resq(:,:),oscr(:,:),oscq(:,:)
  complex(8), allocatable :: denr(:),denq(:),fxc(:,:,:)
  complex(8), allocatable :: fxch(:),fxcw1(:,:),fxcw2(:,:),me(:,:),me2(:,:)
  complex(8), allocatable :: xiout(:,:,:), pmout(:,:,:),resr2(:,:)
  real(8), allocatable :: ev(:),de(:),scisk(:,:)
  integer, allocatable :: widx(:,:,:)
  integer, external :: idxkkp

  write(*,*) 'initializing...'

  emattype=2
  call init0
  call init1
  call init2

  write(*,*) 'preparing...'

  call readfermi
  call xssave0
  call xsgauntgen(max(lmaxapw,lolmax),lmaxemat,max(lmaxapw,lolmax))
  call findgntn0(max(lmaxapwwf,lolmax),max(lmaxapwwf,lolmax),lmaxemat,xsgnt)
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  call ematbdcmbs(emattype)
  allocate(w(nwdf))
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
  
  write(*,*) 'done.'

  n=ngq(iqmt)
  nv=nst1
  nc=nst3
  wsiz=nv*nc*nkptnr
  
  write(*,*) 'n,nv,nc,wsiz',n,nv,nc,wsiz
  
  allocate(wm(nv,nc,nv,nc),widx(nv,nc,nkptnr))
  allocate(de(wsiz),wmat(wsiz,wsiz),wmatq(wsiz,wsiz),me(n,wsiz),me2(n,wsiz))
  allocate(resr(n,wsiz),resq(n,wsiz),oscr(n,n),oscq(n,n),denr(nwdf),denq(nwdf))
  allocate(ev(nstsv))
  allocate(fxc(n,n,nwdf),fxch(3),fxcw1(3,n),fxcw2(n,3))
  allocate(resr2(n,wsiz))


  ! set up indices
  si=0
  do iknr=1,nkptnr
  do iv=1,nv
  do ic=1,nc
    si=si+1
    widx(iv,ic,iknr)=si
write(50,*) 'iknr,ic,iv;s',iknr,ic,iv,si
  end do
  end do
  end do



  ! set up energies and their differences
  do iknr=1,nkptnr
    call getevalsvr(1,nstsv,vkl(:,iknr),ev)
    do iv=1,nv
    do ic=1,nc
      si=widx(iv,ic,iknr)
      de(si)=ev(istocc+ic)-ev(iv)
    end do
    end do
  end do
  
  write(*,*) 'calculating matrix elements....'
  
  ! calculate matrix elements
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3,nv,nc))
  if (allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3,nc,nv))
  if (allocated(deou)) deallocate(deou)
  allocate(deou(nst1,nst3))
  if (allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst3,nst1))
  if (allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1,nst3))
  if (allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst3,nst1))
  allocate(xiout(nv,nc,n))
  allocate(pmout(3,nv,nc))
  allocate(scisk(nst1,nst3))
!!$  if ((fxctype.eq.7).or.(fxctype.eq.8)) then
!!$     call getbsediag
!!$     write(unitout,'("Info(): read diagonal of BSE kernel")')
!!$     write(unitout,'(" mean value : ",2g18.10)') bsed
!!$  end if
  emattype=1
  call ematbdcmbs(emattype)
  call ematrad(iqmt)
  call ematqalloc

  write(*,*) 'calculating matrix elements....'

  do iknr=1,nkptnr
write(unitout,*) 'matrix elements iknr=',iknr
call flushifc(unitout)
    call ematqk1(iqmt,iknr)
    call getdevaldoccsv(iqmt,iknr,iknr,istl1,istu1,istl2,istu2,deou, &
          docc12,scisk)
    call getpemat(iqmt,iknr,'PMAT_SCR.OUT','',m12=xiout,p12=pmout)    
    do iv=1,nv
      do ic=1,nc
        si=widx(iv,ic,iknr)
        me(ig0,si)=pmout(iop,iv,ic)
        if (n.gt.1) then
	  do i=2,n
	    me(i,si)=xiout(iv,ic,i)
	  end do
	end if
      end do
    end do
  end do
  emattype=2
  call ematbdcmbs(emattype)
  
  
  write(*,*) 'done.'

  ! set up W-matrix
  wmat(:,:)=zzero
  wmatq(:,:)=zzero
  do iknr=1,nkptnr; do jknr=iknr,nkptnr
    call getbsemat('SCCLI.OUT',idxkkp(iknr,jknr,nkptnr),nv,nc,wm)
    do iv=1,nv;  do ic=1,nc
      si=widx(iv,ic,iknr)
      dei=de(si)
      do jv=1,nv; do jc=1,nc
        sj=widx(jv,jc,jknr)
	if (si.ne.sj) then
          dej=de(sj)
	  zt1=-wm(iv,ic,jv,jc)
	  if (abs(dei-dej).lt.eps) then
            wmatq(si,sj)=zt1
	  else
            wmat(si,sj)=zt1/(dei-dej)
	  end if
	end if
	
if ((iknr.eq.2).and.(jknr.eq.3).and.(iv.eq.1).and.(ic.eq.1).and.(jv.eq.2).and. &
 (jc.eq.1)) then
   write(unitout,*) '1-1,2-1',wm(iv,ic,jv,jc),de(si),de(sj),de(si)-de(sj)
   write(unitout,*) 'indices:',si,sj
end if	
	
	
      end do; end do
    end do; end do
  end do; end do
  
  do si=1,wsiz
    do sj=si+1,wsiz
      wmat(sj,si)=-conjg(wmat(si,sj))
      wmatq(sj,si)=conjg(wmatq(si,sj))
    end do
  end do  
  
  write(*,*) 'setting up residuals...'
  
  ! set up residuals
  resr=transpose(matmul(wmat,conjg(transpose(me))))
  resq=transpose(matmul(wmatq,conjg(transpose(me))))

!@@@@@@@@@@@@@@@@@@@@@@@@

si=widx(1,1,3); sj=widx(2,1,2)
write(unitout,*) 'wmat:k=1,1-1;k=1,2-1:',si,sj,wmat(widx(1,1,1),widx(2,1,1)),&
     de(widx(1,1,1)),de(widx(2,1,1)),de(widx(1,1,1))-de(widx(2,1,1))
write(unitout,*) 'wmat:k=3,1-1;k=2,2-1',wmat(si,sj),de(si),de(sj),de(si)-de(sj)
write(unitout,*) 'resr:k=3,1-1',resr(1:2,si)
write(unitout,*) 'resr:k=2,1-2',resr(1:2,sj)

  ! deallocate the wmat arrays
  deallocate(wmat,wmatq)
  
  write(*,*) 'setting up kernel...'

  t1=2.d0/(nkptnr*omega)
  fxc(:,:,:)=zzero
  do si=1,wsiz

write(unitout,*) 'si=',si
call flushifc(unitout)

    oscr(:,:)=zzero
    call xszoutpr3(n,n,zone,me(:,si),resr(:,si),oscr)
    !!!call ZGERU(n,n,zone,me(:,si),1,resr(:,si),1,oscr,n)
    ! add Hermitian transpose
    oscr=oscr+conjg(transpose(oscr))
    denr(:)=1.d0/(w(:)-de(si)+zi*broad)
    do iw=1,nwdf
      fxc(:,:,iw)=fxc(:,:,iw)+t1*denr(iw)*oscr(:,:)
    end do
  end do
  
  
  write(*,*) 'writing out kernel...'


  ! write out kernel
  inquire(iolength=recl) fxch,fxcw1,fxcw2,fxc(:,:,1)  
  call genfilname(basename='FXC_BSE',asc=.false.,bzsampl=0,&
       acont=acont,nar=.not.aresdf,iqmt=iqmt,filnam=filnam)
  call getunit(un)
  open(un,file=trim(filnam),form='unformatted',action='write', &
       status='replace',access='direct',recl=recl)
  do iw=1,nwdf
     fxch(:)=fxc(1,1,iw)
     if (n.gt.1) then
       do i=1,noptc
         fxcw1(i,:)=fxc(1,:,iw)
         fxcw2(:,i)=fxc(:,1,iw)
       end do
     end if
     write(un,rec=iw) fxch,fxcw1,fxcw2,fxc(:,:,iw)
     write(8888,'(i8,g18.10,2x,3g18.10)') iw,dble(w(iw)),fxc(1,1,iw),aimag(w(iw))  
  end do
  close(un)
  
  
  do iw=1,nwdf,10
     do i=-noptc,n
        do j=-noptc,n
	   if ((i.gt.0).and.(j.gt.0)) then
             write(9999,'(3i6,3g18.10)') iw,i,j,fxc(i,j,iw), &
                abs(fxc(i,j,iw))
	   else
	      write(9999,'(3i6,3g18.10)') iw,i,j,fxc(1,1,iw),abs(fxc(1,1,iw))
	   end if
        end do
     end do
  end do    
end subroutine kernxc_bse2
!EOC
