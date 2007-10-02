
subroutine linoptold
  use modmain
  use modtddft
  use modpar
  use m_getunit
  integer :: un

  if (rank.eq.1) then
     call linopt
  end if

  call getunit(un)
  call barrier(rank=rank,nproc=nproc,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info(linoptold): linear optics (main version) &
       &finished"

end subroutine linoptold
