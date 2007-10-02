
subroutine tetgather
  use modmain
  use modtddft
  use modpar
  use m_gettetcw
  use m_puttetcw
  use m_filedel
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tetgather'
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
     ! file extension for q-point
     write(filext,'("_Q",i5.5,".OUT")') iq
     do ik=1,nkpt
        do iv=1,nstval
           do ic=1,nstcon

              ! collect weights from processes
              do iproc=1,nproc
                 write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
                 call getrange(iproc,nproc,nwdf,wpari,wparf)
                 nwdfp=wparf-wpari+1
                 allocate(cwp(nwdfp),cwap(nwdfp),cwsurfp(nwdfp))

                 call gettetcw(iq,ik,iv,ic,nwdfp,'.TETW'//trim(filextp),&
                      cwp,cwap,cwsurfp)

                 cw(wpari:wparf)=cwp(:)
                 cwa(wpari:wparf)=cwap(:)
                 cwsurf(wpari:wparf)=cwsurfp(:)
                 deallocate(cwp,cwap,cwsurfp)
              end do ! iproc

              ! write weights
              call puttetcw(iq,ik,iv,ic,'TETW'//trim(filext),cw,cwa,cwsurf)

           end do ! ic
        end do ! iv
     end do ! ik

     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
        call filedel('.TETW'//trim(filextp))
     end do
     
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method gathered for q-point:',iq
     
  end do ! iq

  deallocate(cw,cwa,cwsurf)

  ! restore offset
  vkloff = vkloff_save
  write(filext,'(".OUT")')

end subroutine tetgather
