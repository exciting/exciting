
subroutine dfgather
  use modmain
  use modtddft
  use modpar
  use m_filedel
  use m_getx0
  use m_putx0
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'dfgather'
  integer :: n,iq,iw,iproc,recl
  real(8) :: vkloff_save(3)
  complex(8), allocatable :: chi0(:,:),chi0wg(:,:,:),chi0hd(:)
  logical :: tq0

  ! save k-point offset
  vkloff_save = vkloff

  ! loop over q-points
  do iq = 1, nqpt
     tq0 = tq1gamma.and.(iq.eq.1)
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td
     ! size of local field effects
     n = ngq(iq)
     ! allocate
     allocate(chi0(n,n),chi0wg(n,2,3),chi0hd(3))

     ! file extension for q-point
     write(filext,'("_Q",i5.5,".OUT")') iq
     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
        call getrange(iproc,nproc,nwdf,wpari,wparf)
        do iw=wpari,wparf
           ! exponential factor matrix elements
           call getx0(tq0,iq,iw-wpari+1,trim(fnchi0p),trim(filextp),chi0,&
                chi0wg,chi0hd)
           call putx0(tq0,iq,iw,trim(fnchi0),trim(filext),chi0,chi0wg,chi0hd)
        end do
     end do
     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
!!$        call filedel(trim(fnchi0p)//trim(filextp))
     end do

     deallocate(chi0,chi0wg,chi0hd)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sham response &
          &function gathered for q-point:',iq
  end do

  ! restore offset
  vkloff = vkloff_save
  write(filext,'(".OUT")')

end subroutine dfgather
