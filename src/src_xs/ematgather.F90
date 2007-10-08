
subroutine ematgather
  use modmain
  use modtddft
  use modpar
  use m_filedel
  use m_getemat
  use m_putemat
  use m_genfilname
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
     do iproc=1,nproc
        call genfilname(basename='EMAT',iq=iq,nproc=nproc,rank=iproc-1,&
             filnam=fnemat_t)
        call getrange(iproc,nproc,nkpt,kpari,kparf)
        do ik=kpari,kparf
           ! exponential factor matrix elements
           call getemat(iq,ik,.false.,trim(fnemat_t),xiou,xiuo)
           call putemat(iq,ik,.true.,trim(fnemat),xiou,xiuo)
        end do
     end do
     do iproc=1,nproc
        call genfilname(basename='EMAT',iq=iq,nproc=nproc,rank=iproc-1,&
             filnam=fnemat_t)
        call filedel(trim(fnemat_t))
     end do

     deallocate(xiou,xiuo)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Matrix elements of &
          &exponential factor gathered for q-point:',iq
  end do

  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)

end subroutine ematgather
