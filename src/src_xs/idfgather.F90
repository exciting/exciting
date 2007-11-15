
subroutine idfgather
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_filedel
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'idfgather'
  character(256) :: filnam,filnam2
  integer :: n,m,iq,iw,iproc,recl,nc,oct
  logical :: tq0
  real(8) :: vkloff_save(3)
  complex(8), allocatable :: mdf1(:)

  ! save k-point offset
  vkloff_save = vkloff

  allocate(mdf1(nwdf))
  inquire(iolength=recl) mdf1(1)
  call getunit(unit1)

  ! loop over q-points
  do iq = 1, nqpt
     tq0 = tq1gamma.and.(iq.eq.1)
     ! number of components (3 for q=0)
     nc=1
     if (tq0) nc=3

     ! matrix size for local field effects
     n=ngq(iq)
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td

     do m=1,n,max(n-1,1)

        ! loop over longitudinal components for optics
        do oct=1,nc
           do iproc=0,procs-1
              wpari=firstofset(iproc,nwdf)
              wparf=lastofset(iproc,nwdf)
              ! filename for proc
              call genfilname(basename='IDF',bzsampl=bzsampl,acont=acont,&
                   nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,tq0=tq0,oc=oct,&
                   iq=iq,procs=procs,rank=iproc,filnam=filnam2)

              open(unit1,file=trim(filnam2),form='unformatted',&
                   action='read',status='old',access='direct',recl=recl)
              do iw=wpari,wparf
                 read(unit1,rec=iw-wpari+1) mdf1(iw)
              end do
              close(unit1)
           end do ! iproc

           ! write to file
           call genfilname(basename='IDF',bzsampl=bzsampl,&
                acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
                tq0=tq0,oc=oct,iq=iq,filnam=filnam)
           open(unit1,file=trim(filnam),form='unformatted', &
                action='write',status='replace',access='direct',recl=recl)
           do iw=1,nwdf
              write(unit1,rec=iw) mdf1(iw)
           end do
           close(unit1)

           ! remove partial files
           do iproc=0,procs-1
              call genfilname(basename='IDF',bzsampl=bzsampl,acont=acont,&
                   nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,tq0=tq0,oc=oct,&
                   iq=iq,procs=procs,rank=iproc,filnam=filnam2)
              call filedel(trim(filnam2))
           end do

        end do ! oct

     end do ! m

     write(unitout,'(a,i8)') 'Info('//thisnam//'): inverse dielectric &
          &function gathered for q-point:',iq
  end do

  deallocate(mdf1)

  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)

end subroutine idfgather
