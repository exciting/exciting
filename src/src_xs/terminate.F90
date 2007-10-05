
subroutine terminate
  use modpar
  implicit none

  ! abort MPI if necessary
#ifdef MPI
  call mpi_abort(mpi_comm_world,ierr)

  if (ierr.eq.0) then
     write(*,'(a)') 'MPI abort clean. STOP'
  else
     write(*,'(a)') 'MPI abort unclean - zombie processes might remain! STOP'
  end if
#endif

#ifndef MPI
  write(*,'(a)') 'Abort clean. STOP'
#endif

  ! stop
  stop

end subroutine terminate


subroutine terminate_inqr(str)
  use m_getunit
  implicit none
  ! arguments
  character(*) :: str
  ! local variables
  logical :: exis
  ! check for terminating...
  inquire(file='TERMINATE',exist=exis)
  if (exis) then
     write(*,'(a)') 'Error: user termination of program in routine: '// &
          trim(str)
     call terminate
  end if

end subroutine terminate_inqr
