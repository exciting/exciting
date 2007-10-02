
subroutine df
  use modmain
  use modtddft
  use modpar
  use m_dfq
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'df'
  integer :: iq,un,j

  if (calledtd.eq.1) call init0

  ! initialise universal variables
  call init1

  ! save Gamma-point variables
  call tdsave0

  ! initialize q-point set
  call init2td

  ! w-point interval for process
  call getrange(rank,nproc,nwdf,wpari,wparf)

  ! loop over q-points
  do iq = 1, nqpt
     ! call for q-point
     if (.not.gather) call dfq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sahm response &
          &function finished for q-point:',iq
  end do

  ! synchronize
  call getunit(un)
  if (.not.gather) call barrier(rank=rank,nproc=nproc,un=un,async=0, &
       string='.barrier')

  if ((nproc.gt.1).and.(rank.eq.1)) call dfgather

  if (.not.gather) call barrier(rank=rank,nproc=nproc,un=un,async=1, &
       string='.barrier')

  write(unitout,'(a)') "Info("//trim(thisnam)//"): Kohn-Sham response &
       &function finished"

  if (gather) then
     write(unitout,'(a)') "Info("//trim(thisnam)//"): gather option: &
          &exiting program"
     call tdepilog
     call terminate
  end if

end subroutine df
