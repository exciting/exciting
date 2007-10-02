
subroutine argparse()
  use modmain
  use modpar
  implicit none
  ! local variables
  character(256) :: str
  
#ifndef MPI
  ! read rank
  call getarg(1,str)
  if (trim(str).eq.'') then
     rank=1
  else
     str=adjustl(str)
     read(str,*) rank
  end if

  ! read nproc
  call getarg(2,str)
  if (trim(str).eq.'') then
     nproc=1
  else
     str=adjustl(str)
     read(str,*) nproc
  end if
#endif

end subroutine argparse
