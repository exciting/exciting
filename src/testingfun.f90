subroutine testingfun
  use modmpi
  use modscl

  implicit none

  integer(4) :: npg
  type(procgroup) :: mygroup
  type(blacsinfo) :: subblacs

  npg = mpiglobal%procs/2
  
  call setup_proc_groups(npg,mygroup) 

  write(*,*) "Global rank:",mpiglobal%rank, "Mygroup%mpi%com",&
   & mygroup%mpi%comm, "mygroup%mpi%procs", mygroup%mpi%procs,&
   & "mygroup%mpi%rank", mygroup%mpi%rank

  call setupblacs(mygroup%mpi, "column", subblacs)

  write(*,*) "Global rank:",mpiglobal%rank, "subblacs%context",&
   & subblacs%context, "subblacs%nprocs", subblacs%nprocs,&
   & "subblacs%mypcol", subblacs%mypcol

  call barrier
  call terminate

end subroutine testingfun
