subroutine testingfun
  use modmpi
  use modscl
  use mod_kpoint, only: nkpt, wkpt
  use mod_qpoint, only: nqpt, vql

  implicit none

  integer(4) :: i
  
  nqpt = 8

  write(*,*) "nofset(",mpiglobal%rank,",",nqpt,",",mpiglobal%procs,")=",&
    & nofset(mpiglobal%rank, nqpt, mpiglobal%procs)
  write(*,*) "nofset(",mpiglobal%rank,",",nqpt,")=",&
    & nofset(mpiglobal%rank, nqpt)
  write(*,*) "firstofset(",mpiglobal%rank,",",nqpt,",",mpiglobal%procs,")=",&
    & firstofset(mpiglobal%rank, nqpt, mpiglobal%procs)
  write(*,*) "firstofset(",mpiglobal%rank,",",nqpt,")=",&
    & firstofset(mpiglobal%rank, nqpt)
  write(*,*) "lastofset(",mpiglobal%rank,",",nqpt,",",mpiglobal%procs,")=",&
    & lastofset(mpiglobal%rank, nqpt, mpiglobal%procs)
  write(*,*) "lastofset(",mpiglobal%rank,",",nqpt,")=",&
    & lastofset(mpiglobal%rank, nqpt)

  do i=1, nqpt
    write(*,*) "procofindex(",i,",",nqpt,",",mpiglobal%procs,")=",&
      & procofindex(i, nqpt, mpiglobal%procs)
    write(*,*) "procofindex(",i,",",nqpt,")=",&
      & procofindex(i, nqpt, mpiglobal%procs)
  end do

  write(*,*) "lastproc(",1,",",nqpt,",",mpiglobal%procs,")=",&
    & lastproc(1, nqpt, mpiglobal%procs)

 ! integer(4) :: npg
 ! type(procgroup) :: mygroup
 ! type(blacsinfo) :: subblacs

 ! npg = mpiglobal%procs/2
 ! 
 ! call setup_proc_groups(npg,mygroup) 

 ! write(*,*) "Global rank:",mpiglobal%rank, "Mygroup%mpi%com",&
 !  & mygroup%mpi%comm, "mygroup%mpi%procs", mygroup%mpi%procs,&
 !  & "mygroup%mpi%rank", mygroup%mpi%rank

 ! call setupblacs(mygroup%mpi, "column", subblacs)

 ! write(*,*) "Global rank:",mpiglobal%rank, "subblacs%context",&
 !  & subblacs%context, "subblacs%nprocs", subblacs%nprocs,&
 !  & "subblacs%mypcol", subblacs%mypcol

 ! call barrier
 ! call terminate

end subroutine testingfun
