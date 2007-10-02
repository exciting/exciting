
subroutine ematgather
  use modmain
  use modtddft
  use modpar
  use m_filedel
  use m_getemat
  use m_putemat
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'ematgather'
  integer :: iq,ik,ikr,iproc,recl
  real(8) :: vkloff_save(3)

  ! save k-point offset
  vkloff_save = vkloff

  ! allocate matrix elements array
  if (allocated(xiou)) deallocate(xiou)
  if (allocated(xiuo)) deallocate(xiuo)

  ! loop over q-points
  do iq = 1, nqpt
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td
     allocate(xiou(nstval,nstcon,ngq(iq)))
     allocate(xiuo(nstcon,nstval,ngq(iq)))
     ! file extension for q-point
     write(filext,'("_Q",i5.5,".OUT")') iq
     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
        call getrange(iproc,nproc,nkpt,kpari,kparf)
        do ik=kpari,kparf
           ! exponential factor matrix elements
           call getemat(iq,ik,.false.,trim(fnemat)//trim(filextp),xiou,xiuo)
           call putemat(iq,ik,.true.,'EMAT'//trim(filext),xiou,xiuo)
        end do
     end do
     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
        call filedel(trim(fnemat)//trim(filextp))
     end do

     deallocate(xiou,xiuo)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Matrix elements of &
          &exponential factor gathered for q-point:',iq
  end do

  ! restore offset
  vkloff = vkloff_save
  write(filext,'(".OUT")')

end subroutine ematgather
