
subroutine devalsvgather
  use modmain
  use modtddft
  use modpar
  use m_filedel
  use m_getdevalsv
  use m_putdevalsv
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'devalsvgather'
  integer :: iq,ik,ikr,iproc,recl
  real(8) :: vkloff_save(3)

  ! save k-point offset
  vkloff_save = vkloff

  ! allocate matrix elements array
  if (allocated(deou)) deallocate(deou)
  if (allocated(deuo)) deallocate(deuo)

  ! loop over q-points
  do iq = 1, nqpt
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td
     allocate(deou(nstval,nstcon))
     allocate(deuo(nstcon,nstval))
     ! file extension for q-point
     write(filext,'("_Q",i5.5,".OUT")') iq
     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
        call getrange(iproc,nproc,nkpt,kpari,kparf)
        do ik=kpari,kparf
           ! exponential factor matrix elements
           call getdevalsv(iq,ik,.false.,trim(fndevalsv)//trim(filextp), &
                deou,deuo)
           call putdevalsv(iq,ik,.true.,'DEVALSV'//trim(filext),deou,deuo)
        end do
     end do
     do iproc=1,nproc
        write(filextp,'("_Q",i5.5,"_par",i3.3,".OUT")') iq, iproc
        call filedel(trim(fndevalsv)//trim(filextp))
     end do

     deallocate(deou,deuo)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn-Sham eigenvalue &
          &differences gathered for q-point:',iq
  end do

  ! restore offset
  vkloff = vkloff_save
  write(filext,'(".OUT")')

end subroutine devalsvgather
