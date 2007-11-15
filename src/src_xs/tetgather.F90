
subroutine tetgather
  use modmain
  use modxs
  use modmpi
  use m_gettetcw
  use m_puttetcw
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tetgather'
  character(256) :: filnam,filnam_t
  integer :: n,iq,iw,iproc,ik,iv,ic,nwdfp
  real(8) :: vkloff_save(3)
  real(8), allocatable :: cw(:),cwa(:),cwsurf(:)
  real(8), allocatable :: cwp(:),cwap(:),cwsurfp(:)
  logical :: tq0

  ! save k-point offset
  vkloff_save = vkloff

  allocate(cw(nwdf),cwa(nwdf),cwsurf(nwdf))

  ! loop over q-points
  do iq = 1, nqpt
     tq0 = tq1gamma.and.(iq.eq.1)
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td
     ! file name for output file
     call genfilname(basename='TETW',iq=iq,filnam=filnam)
     do ik=1,nkpt
        do iv=1,nstval
           do ic=1,nstcon

              ! collect weights from processes
              do iproc=0,procs-1

                 ! filename for input file
                 call genfilname(basename='TETW',iq=iq,rank=iproc,&
                      procs=procs,filnam=filnam_t)
                 
                 wpari=firstofset(iproc,nwdf)
                 wparf=lastofset(iproc,nwdf)

                 nwdfp=wparf-wpari+1
                 allocate(cwp(nwdfp),cwap(nwdfp),cwsurfp(nwdfp))

                 call gettetcw(iq,ik,iv,ic,nwdfp,trim(filnam_t),&
                      cwp,cwap,cwsurfp)

                 cw(wpari:wparf)=cwp(:)
                 cwa(wpari:wparf)=cwap(:)
                 cwsurf(wpari:wparf)=cwsurfp(:)
                 deallocate(cwp,cwap,cwsurfp)
              end do ! iproc

              ! write weights
              call puttetcw(iq,ik,iv,ic,trim(filnam),cw,cwa,cwsurf)

           end do ! ic
        end do ! iv
     end do ! ik

     do iproc=0,procs-1
        call genfilname(basename='TETW',iq=iq,rank=rank,procs=procs,&
             filnam=filnam_t)
        call filedel(trim(filnam_t))
     end do
     
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method gathered for q-point:',iq
     
  end do ! iq

  deallocate(cw,cwa,cwsurf)

  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)

end subroutine tetgather
