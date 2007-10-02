
subroutine idfgather
  use modmain
  use modtddft
  use modtetra
  use modpar
  use m_filedel
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'idfgather'
  character(256) :: filnam2,filnam3,str
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

     ! file extension for q-point
     write(filext,'("_Q",i5.5,".OUT")') iq

     do m=1,n,max(n-1,1)

        ! loop over longitudinal components for optics
        do oct=1,nc

           ! string for xc-kernel type
           write(str,'(i2.2)') fxctype
           if (tq0) write(str,'(i2.2,"_OC",i2.2)') fxctype,11*oct
           str='_FXC'//trim(str)
           if (m.eq.1) str='_NLF'//trim(str)
           if (.not.aresdf) str='_NAR'//trim(str)
           if (acont) str='_AC'//trim(str)
           if (tetra) str='_TET'//trim(str)
           filnam2='IDF'//trim(str)
           filnam3=trim(filnam2)
           if (nproc.gt.1) filnam2='.'//trim(filnam2)

           do iproc=1,nproc
              write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
              call getrange(iproc,nproc,nwdf,wpari,wparf)
              open(unit1,file=trim(filnam2)//trim(filextp),form='unformatted',&
                   action='read',status='old',access='direct',recl=recl)
              do iw=wpari,wparf
                 read(unit1,rec=iw-wpari+1) mdf1(iw)
              end do
              close(unit1)
           end do ! iproc

           ! write to file
           open(unit1,file=trim(filnam3)//trim(filext),form='unformatted', &
                action='write',status='replace',access='direct',recl=recl)
           do iw=1,nwdf
              write(unit1,rec=iw) mdf1(iw)
           end do
           close(unit1)

           ! remove partial files
           do iproc=1,nproc
              write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
              call filedel(trim(filnam2)//trim(filextp))
           end do

        end do ! oct

     end do ! m

     write(unitout,'(a,i8)') 'Info('//thisnam//'): inverse dielectric &
          &function gathered for q-point:',iq
  end do

  deallocate(mdf1)

  ! restore offset
  vkloff = vkloff_save
  write(filext,'(".OUT")')

end subroutine idfgather
