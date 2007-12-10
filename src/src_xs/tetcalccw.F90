
subroutine tetcalccw
  use modmain
  use modxs
  use modmpi
  use modtetra
  use m_tetcalccwq
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tetcalccw'
  integer :: iq,un,j
  logical :: tlfe, tet

  if (calledxs.eq.1) call init0

  ! initialise universal variables
  tet=tetra
  tetra=.true.
  call init1

  ! save Gamma-point variables
  call tdsave0

  ! initialize q-point set
  call init2xs

  ! w-point interval for process
  wpari=firstofset(rank,nwdf)
  wparf=lastofset(rank,nwdf)

  ! loop over q-points
  do iq = 1, nqpt
     ! call for q-point
     call tetcalccwq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method finished for q-point:',iq
  end do

  ! synchronize
  call getunit(un)
  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  if ((procs.gt.1).and.(rank.eq.0)) call tetgather

  call barrier(rank=rank,procs=procs,un=un,async=1,string='.barrier')

  tetra=tet

  write(unitout,'(a)') "Info("//trim(thisnam)//"): weights for tetrahedron &
       &method finished"

end subroutine tetcalccw
