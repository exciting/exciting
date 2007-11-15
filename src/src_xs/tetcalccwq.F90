
module m_tetcalccwq
  implicit none
contains

  subroutine tetcalccwq(iq)
    use modmain
    use modxs
    use modtetra
    use modmpi
    use m_genwgrid
    use m_getunit
    use m_filedel
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'tetcalccwq'
    character(256) :: filnam,filnamt
    complex(8), allocatable :: w(:)
    real(8), parameter :: epstetra=1.d-8
    real(8), allocatable :: eb(:,:)
    real(8), allocatable :: wreal(:)
    real(8), allocatable :: cwsurft2(:,:),cwt2(:,:),cwat2(:,:)
    real(8), allocatable :: cwsurft1(:),cwt1(:),cwat1(:)
    real(8), allocatable :: cwsurf(:,:,:),cw(:,:,:),cwa(:,:,:)
    real(8) :: wt, vkloff_save(3)
    integer :: ik,ist,iv,ic
    integer :: iw,wi,wf,nwdfp,un,un2,recl,recl2,irec,irec2
    logical :: exis,tq0,tetrat

!!$    ! debug output in tetrahedron integration library
!!$    call tetrasetdbglv(1)
!!$    ! safer pointer handling in tetrahedron integration library
!!$    call tetrasetpointerhandling(1)

    ! init1 should be called for settings in libbzint


    tq0 = tq1gamma.and.(iq.eq.1)
    ! save k-point offset
    vkloff_save = vkloff
    ! shift k-mesh by q-point    
    vkloff(:)=qvkloff(:,iq)

    ! initial and final w-point
    wi=wpari
    wf=wparf
    nwdfp=wf-wi+1

    ! allocations
    allocate(w(nwdf))
    allocate(wreal(nwdfp))
    ! generate complex energy grid
    call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
    wreal(:)=w(wi:wf)
    if (wreal(1).lt.epstetra) wreal(1)=epstetra

    ! set q-dependent file extension
    call genfilname(iq=iq,setfilext=.true.)

    ! generate filenames
    call genfilname(basename='TETW',iq=iq,rank=rank,procs=procs,&
         filnam=filnam)
    call genfilname(basename='TETWT',iq=iq,rank=rank,procs=procs,&
         filnam=filnamt)
    
    ! calculate k+q and G+k+q related variables
    call init1td

    ! tetrahedron method
    if (.not.tq0) then
       write(unitout,'(a)') 'Error('//trim(thisnam)//'): non-Gamma q-point &
            &and tetrahedron method chosen'
       tetrat=.false.
       call terminate
    end if

    ! read Fermi energy
    call genfilname(iq=0,setfilext=.true.)
    call readfermi
    call genfilname(iq=iq,setfilext=.true.)

    ! allocate arrays
    allocate(eb(nstsv,nkpt))

    ! get the eigenvalues from file
    do ik=1,nkpt
       call getevalsv(vkl(1,ik),evalsv(1,ik))
    end do
    eb(:,:)=evalsv(:,:)
    ! scissors shift
    do ik=1,nkpt
       do ist=nstval+1,nstsv
          eb(ist,ik)=eb(ist,ik)+scissor
       end do
    end do
    deallocate(evalsv)
    allocate(cw(nstsv,nstsv,nkpt))
    allocate(cwa(nstsv,nstsv,nkpt))
    allocate(cwsurf(nstsv,nstsv,nkpt))
    allocate(cwsurft2(nstval,nstcon),cwt2(nstval,nstcon),cwat2(nstval,nstcon))

    call getunit(un)
    inquire(iolength=recl) cw(:nstval,nstval+1:,1),cwa(:nstval,nstval+1:,1),&
         cwsurf(:nstval,nstval+1:,1)
    open(un,file=trim(filnamt),form='unformatted',&
         action='write',status='replace',access='direct',recl=recl)
    ! calculate weights
    do iw=1,nwdfp
       wt=wreal(iw)
       if (abs(wt).lt.epstetra) wt=epstetra
       ! switch 2 below in tetcw defines bulk integration for real part
       ! resonant contribution
       call tetcw(nkpt,ntet,nstsv,wtet,eb,tnodes,link,tvol,efermi, &
            wt,2,cw)
       ! anti-resonant contribution
       call tetcw(nkpt,ntet,nstsv,wtet,eb,tnodes,link,tvol,efermi, &
            -wt,2,cwa)
       ! switch 4 below in tetcw defines surface integration for imag. part
       call tetcw(nkpt,ntet,nstsv,wtet,eb,tnodes,link,tvol,efermi, &
            wt,4,cwsurf)
       do ik=1,nkpt
          irec=(ik-1)*nwdfp+iw
          cwsurft2(:,:)=cwsurf(:nstval,nstval+1:,ik)
          cwt2(:,:)=cw(:nstval,nstval+1:,ik)
          cwat2(:,:)=cwa(:nstval,nstval+1:,ik)
          write(un,rec=irec) cwt2,cwat2,cwsurft2
       end do
       if (iw <= nwdf/procs) then
          ! synchronize for common number of w-points to all processes
          call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')
       end if
    end do
    close(un)
    deallocate(cw,cwa,cwsurf)
    allocate(cw(nwdfp,nstval,nstcon))
    allocate(cwa(nwdfp,nstval,nstcon))
    allocate(cwsurf(nwdfp,nstval,nstcon))
    allocate(cwsurft1(nwdfp),cwt1(nwdfp),cwat1(nwdfp))

    open(un,file=trim(filnamt),form='unformatted',action='read',&
         status='old',access='direct',recl=recl)
    call getunit(un2)
    inquire(iolength=recl2) cw(:,1,1),cwa(:,1,1),cwsurf(:,1,1)
    open(un2,file=trim(filnam),form='unformatted',&
         action='write',status='replace',access='direct',recl=recl2)
    irec=0
    irec2=0
    do ik=1,nkpt
       do iw=1,nwdfp
          irec=irec+1
          read(un,rec=irec) cwt2,cwat2,cwsurft2
          cw(iw,:,:)=cwt2(:,:)
          cwa(iw,:,:)=cwat2(:,:)
          cwsurf(iw,:,:)=cwsurft2(:,:)
       end do
       do iv=1,nstval
          do ic=1,nstcon
             irec2=irec2+1
             cwsurft1(:)=cwsurf(:,iv,ic)
             cwt1(:)=cw(:,iv,ic)
             cwat1(:)=cwa(:,iv,ic)
             write(un2,rec=irec2) cwt1,cwat1,cwsurft1
          end do
       end do
    end do
    close(un)
    call filedel(trim(filnamt))
    close(un2)

    deallocate(cwt2,cwat2,cwsurft2)
    deallocate(cw,cwa,cwsurf,eb)
    deallocate(cwt1,cwat1,cwsurft1)
    deallocate(w,wreal)
    write(unitout,'(a)') 'Info('//trim(thisnam)//'): weights for tetrahedron &
         &method finished.'

    ! restore offset
    vkloff(:) = vkloff_save(:)
    ! restore file extension
    call genfilname(setfilext=.true.)

  end subroutine tetcalccwq

end module m_tetcalccwq
