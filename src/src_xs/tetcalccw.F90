
subroutine tetcalccw
  use modmain
  use modtddft
  use modpar
  use modtetra
  use m_tetcalccwq
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tetcalccw'
  integer :: iq,un,j
  logical :: tlfe, tet

  if (calledtd.eq.1) call init0

  ! initialise universal variables
  tet=tetra
  tetra=.true.
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
     call tetcalccwq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method finished for q-point:',iq
  end do

  ! synchronize
  call getunit(un)
  call barrier(rank=rank,nproc=nproc,un=un,async=0,string='.barrier')

  if ((nproc.gt.1).and.(rank.eq.1)) call tetgather

  call barrier(rank=rank,nproc=nproc,un=un,async=1,string='.barrier')

  tetra=tet

  write(unitout,'(a)') "Info("//trim(thisnam)//"): weights for tetrahedron &
       &method finished"

end subroutine tetcalccw
