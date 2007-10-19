
subroutine linoptold
  use modmain
  use modtddft
  use modmpi
  use m_getunit
  integer :: un

  if (rank == 0) then
     call linopt
  end if

  call getunit(un)
  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info(linoptold): linear optics (main version) &
       &finished"

end subroutine linoptold
