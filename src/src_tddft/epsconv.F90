
subroutine epsconv
  use modmain
  use modtddft
  use modtetra
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'epsconv'
  character(256) :: filnam,str
  integer, parameter :: numlines_top=48
  integer :: iq,iw,iwp,j,m,n,oct,nc,un
  logical :: exis,tq0
  real(8), parameter :: epsc=1.d-8
  real(8) :: t0
  real(8), allocatable :: w(:), epst(:,:),lor(:),f(:),f1(:),g(:),g1(:),cf(:,:)
  complex(8), allocatable :: eps(:)

  call init0
  call init1
  call init2td

  allocate(w(nwdf),epst(nwdf,2),eps(nwdf),lor(nwdf))
  allocate(f(nwdf),f1(nwdf),g(nwdf),g1(nwdf),cf(3,nwdf))
  t0=1.d0
  if (tevout) t0=h2ev

  ! loop over q-points
  do iq=1,nqpt
     ! matrix size for local field effects
     n=ngq(iq)
     tq0 = tq1gamma.and.(iq.eq.1)
     nc=1
     if (tq0) nc=3
     ! neglect/include local field effects
     do m=1,n,max(n-1,1)
        ! loop over longitudinal components for optics
        do oct=1,nc
           ! file extension
           write(filext,'("_Q",i5.5,".OUT")') iq
           ! string for xc-kernel type
           write(str,'(i2.2)') fxctype
           if (tq0) write(str,'(i2.2,"_OC",i2.2)') fxctype,11*oct
           str='_FXC'//trim(str)
           if (m.eq.1) str='_NLF'//trim(str)
           if (.not.aresdf) str='_NAR'//trim(str)
           str='_TET'//trim(str)

           filnam='EPSILON'//trim(str)
           inquire(file=trim(filnam)//trim(filext),exist=exis)
           if (.not.exis) then
              write(*,*) 'Error('//trim(thisnam)//'): file does not exist: '&
                   //trim(filnam)//trim(filext)
              call terminate
           end if
           call getunit(un)
           open(unit=un,file=trim(filnam)//trim(filext),form='formatted',&
                action='read',status='old')

           ! read comments at the top of the file (check the number of lines
           ! from version to version!)
           do j=1,numlines_top
              read(un,*)
           end do

           ! read energies and Re and Im of eps
           do iw=1,nwdf
              read(un,*) w(iw),epst(iw,1),epst(iw,2)
           end do
           close(un)
           open(unit=un,file=trim(filnam)//'_sL'//trim(filext),&
                form='formatted',action='write',status='replace')

           w(:)=w(:)/t0
           ! complex dielectric function
           eps(:)=epst(:,1)+zi*epst(:,2)
           ! convolution with Lorentzian
           do iw=1,nwdf
              do iwp=1,nwdf
                 ! standard Lorentzian with peak at w(iw)
!!$                 lor(iwp)=(1/pi)*brdtd/((w(iw)-w(iwp))**2+brdtd**2)
                 ! antisymmetric Lorentzian at w(iw) and -w(iw)
                 ! with norm arctan(w/brdtd) to assure zero crossing
                 lor(iwp)=(1.d0/(2.d0*atan(w(iw)/brdtd)))*( &
                      brdtd/((w(iw)-w(iwp))**2+brdtd**2) - &
                      brdtd/((-w(iw)-w(iwp))**2+brdtd**2) )
                 if (w(iw) < epsc) lor(iwp) = 0.d0
                 f(iwp)=lor(iwp)*aimag(eps(iwp))
                 f1(iwp)=lor(iwp)*dble(eps(iwp))
              end do
              ! 
              ! do the convolution
              call fderiv(-1,nwdf,w,f,g,cf)
              call fderiv(-1,nwdf,w,f1,g1,cf)

              write(un,'(4g18.10)') w(iw)*t0,g1(nwdf),g(nwdf),&
                   (pi/brdtd)*brdtd**2/((w(iw)-w(iwp))**2+brdtd**2)
           end do ! iw

!!$           !/////////////////////////////////////////////////////////
!!$           g=aimag(eps)
!!$           g1=dble(eps)
!!$           call fsmooth(nsmdos,nwdf,1,g)
!!$           call fsmooth(nsmdos,nwdf,1,g1)
!!$           write(un,'(3g18.10)') (w(iw)*t0,g1(iw),g(iw),iw=1,nwdf)
!!$           !/////////////////////////////////////////////////////////


           close(un)
        end do ! oct
     end do
  end do

  deallocate(w,epst,eps,lor,f,f1,g,g1,cf)

end subroutine epsconv
