! Copyright (C) 2006-2008 C. Ambrosch-Draxl. C. Meisenbichler S. Sagmeister
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
! !MODULE: modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   In case of compiled without MPI support it defines the
!   mpi specific variables such that the code behaves exactly as
!   the unmodified scalar version
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (S. Sagmeister)
!   Added allgatherv interface, August 2010 (S. Sagmeister)
!   Added subroutines/functions to documentation scheme. 2016 (Aurich)
!   Adapted partitioning functions to handle cases with more processes than elements. 2016 (Aurich)
!   Added proc groups functionality. 2016 (Aurich)
!
! TODO(Alex) Issue #23. MPI wrappers are a mix of types and global variables, the latter of which isn't
! required as one can always query MPI variables with MPI calls, there is no overload set for compilation
! in serial, and the routines don't follow single responsibility.
! This needs refactoring if we hope to use more complex MPI distribution schemes.

module modmpi
  use trace, only: trace_back
#ifdef MPI
  use mpi
#endif

  implicit none

  ! MPI info type
  ! Contains basic information regarding
  ! a mpi communicator
  type mpiinfo
    integer(4) :: rank
    integer(4) :: procs
    integer :: comm
    integer(4) :: ierr
    integer(4) :: root = 0   !! Rank of root process 
    logical :: is_root       !! Indicates of the process is root 
  end type mpiinfo

  ! Groups of MPI communicators
  ! connected via inter-communicator
  type procgroup
    ! Total number of process groups
    ! this group belongs to
    integer(4) :: ngroups
    ! Group id
    integer(4) :: id
    ! MPI information for current
    ! process group
    type(mpiinfo) :: mpi
    ! Inter-groups communicator
    type(mpiinfo) :: mpiintercom
  end type procgroup

  ! TODO(Alex) Issue #23. This is a bad idea. Defeats the point of encapsulating an instantiation of the
  ! MPI environment. All MPI variables can be querried with MPI subroutines. Benjamin probably did it
  ! because it was the only way he could get the object passed to BSE without a big refactor of
  ! subroutine APIs
  ! mpiinfo for global scope
  type(mpiinfo) :: mpiglobal

  ! Nodes as procgroup
  type(procgroup) :: mpinodes

  !! Legacy code:
  ! Variables (contained in mpiglobal type)
  integer(4) :: rank
  integer(4) :: procs
  integer(4) :: ierr
  ! Variables (contained in mpinodes)
  integer(4) :: firstinnode_comm

  ! Some parts use these
  logical :: splittfile, firstinnode

  contains

    !> Terminate exciting if condition is false
    subroutine terminate_if_false(condition, message)
      use iso_fortran_env, only: error_unit
      !> Error condition, need to be false to terminate the program and print the message
      logical, intent(in) :: condition
      !> Error message that is printed to the terminal if present and condition is false.
      character(*), optional, intent(in) :: message

      character(256) :: error_message
      
      error_message = 'Error'
      if (present(message)) then
        error_message = error_message//': '//trim(adjustl(message))
      end if 

      if(.not. condition) then
        write(error_unit, *)
        write(error_unit, *) trim(error_message)
        write(error_unit, *)
        call trace_back()
        call terminate()
      end if
    end subroutine terminate_if_false

    !+++++++++++++++++++++++++++++++++++++++++!
    ! Initialization and finalization  of MPI !
    !+++++++++++++++++++++++++++++++++++++++++!
    !BOP
    ! !ROUTINE: initmpi
    ! !INTERFACE:
    subroutine initmpi
    ! !DESCRIPTION:
    !   Initializes MPI and sets procs and rank numbers.
    !   Sets splittfile and initializes the first in node
    !   list. Or if -DMPI is not used sets procs=1 and rank=1
    !   and splittfile=.false.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme (Aurich)
    !EOP
    !BOC
      integer(4) :: ierr
#ifdef MPI
      ! Initialize MPI and get number of processes
      ! and current rank
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, procs, ierr)
      call mpi_comm_rank(mpi_comm_world, rank, ierr)

      ! Set global mpiinfo type
      mpiglobal%procs = procs
      mpiglobal%comm = mpi_comm_world
      mpiglobal%rank = rank
      mpiglobal%ierr = ierr
      mpiglobal%is_root = mpiglobal%rank == mpiglobal%root

      ! Make communicators for
      ! intra- and inter-processor (node) communication.
      call setup_node_groups

      ! Each rank writes its own file (use in pars of GS and XS)
      splittfile = .true.

#endif
#ifndef MPI
      procs = 1
      rank = 0

      mpiglobal%procs = procs
      mpiglobal%rank = rank
      mpiglobal%comm = 0
      mpiglobal%is_root = .true.

      splittfile = .false.
      firstinnode = .true.
#endif
    end subroutine initmpi
    !EOC

    !BOP
    ! !ROUTINE: finitmpi
    ! !INTERFACE:
    subroutine finitmpi
    ! !DESCRIPTION:
    !   If -DMPI calls {\tt mpi\_finalize}, after wainting
    !   for the processes of mpiglobal, mpinodes%mpi and mpinodes%mpiintercom
    !   to finish communicating.
    !
    !   Note: For other communicators you need to make shure that they finished
    !         communication.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Modified to use new mpi types. (Aurich)
    !EOP
    !BOC
      integer(4) :: ierr
      logical :: flag

      ! Wait for everyone to reach this point
      call barrier(mpiglobal)
#ifdef MPI
      call barrier(mpinodes%mpi)
      if(mpinodes%mpiintercom%rank >= 0) then
        call barrier(mpinodes%mpiintercom)
      end if
#endif
#ifdef MPI
      call mpi_finalize(ierr)
      if(ierr /= 0) then
        write(*,*) "Error (finitmpi): ierr =",ierr
      end if
#endif
    end subroutine finitmpi
    !EOC


    !> Terminate an MPI environment
    subroutine terminate_mpi_env(mpi_env, message)
      use iso_fortran_env, only: error_unit
      
      !> MPI environment object
      type(mpiinfo), intent(inout) :: mpi_env
      !> Error message
      character(len=*), optional, intent(in) :: message
      !> Error code for Exciting to return to the invoking environment
      integer, parameter :: error_code = 101
      
#ifdef MPI
      if(mpi_env%rank == 0) then
         if(present(message)) write(error_unit, *) trim(adjustl(message))
      end if
      call mpi_abort(mpi_env%comm, error_code, mpi_env%ierr)
#else
      if(present(message)) write(error_unit, *) trim(adjustl(message))
      stop
#endif
    end subroutine terminate_mpi_env

    
    !> Terminate global MPI environment   
    !>
    !> Developers should not use this, in favour of terminate_mpi_env
    subroutine terminate(user_msg)
      !> Optional error message
      character(len=*), intent(in), optional :: user_msg

      !> Error code for exciting to return to the invoking environment
      integer, parameter :: error_code = 101
      !> Error integer 
      integer :: ierr
      !> Local error message
      character(256) :: err_msg

      if (present(user_msg)) then
        err_msg = user_msg
      else 
        err_msg = "exciting has terminated"
      end if

#ifdef MPI
      if(mpiglobal%rank == 0) then
         write(*,'(a)') trim(adjustl(err_msg))
      end if
      
      call mpi_abort(mpi_comm_world, error_code, ierr)
      
      if(ierr == 0) then
         write(*, '(a)') trim(adjustl(err_msg))
      else
         write(*, '(a)') trim(adjustl(err_msg))//' - zombie processes might remain!'
      end if
#else
      write(*, '(a)') err_msg
      stop
#endif
    end subroutine terminate
 

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Partitioning of N elements to P processes !
    ! in continuous blocks. Each element is     !
    ! associated to one and only one process.   !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: nofset
    ! !INTERFACE:
    function nofset(myrank, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: myrank  ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! integer(4), optional :: nprocs ! Number of processes in communicator
    ! OUT:
    ! integer(4) :: nofset  ! Number of elements for that rank
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the number of elements $N_\text{el}(p)$ a given process is responsible for. \\
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow N_\text{el}(0)=4, N_\text{el}(1)=3, N_\text{el}(2)=3$\\
    !   Example:\\
    !   $N_\text{el}=2, N_\text{p}=3 \rightarrow N_\text{el}(0)=1, N_\text{el}(1)=1, N_\text{el}(2)=0$
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
      integer(4) :: nofset
      integer(4), intent(in) :: myrank, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: np

      ! Sanity checks
      if( myrank < 0 ) then
        write(*,*) "nofset (Error): myrank < 0"
        call terminate
      end if
      if( set < 1 ) then
        write(*,*) "nofset (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "nofset (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "nofset (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(np < myrank+1) then
        write(*,*) "nofset (Error): np < myrank+1"
        call terminate
      end if

      ! Compute number of elements on current rank
      nofset = set / np
      if((mod(set, np) > myrank)) nofset = nofset + 1

    end function nofset
    !EOC

    !BOP
    ! !ROUTINE: firstofset
    ! !INTERFACE:
    function firstofset(myrank, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: myrank  ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! integer(4), optional :: nprocs  ! Number of processes in commuincator
    ! OUT:
    ! integer(4) :: firstofset ! Index of the total set for the first index
    !                          ! of the current subset
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the index of the fist element $i_\text{el}(p)$ a given process is responsible for.
    !   If there are more processes than elements a process responsible for no element
    !   gets the assignment $i_\text{el}(p > N_\text{el}-1) = 0$.\\
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow i_\text{el}(0)=1, i_\text{el}(1)=5, i_\text{el}(2)=8$\\
    !   Example:\\
    !   $N_\text{el}=2, N_\text{p}=4 \rightarrow i_\text{el}(0)=1, i_\text{el}(1)=2, i_\text{el}(2)=0, i_\text{el}(3)=0$
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Changed behaviour if there are more processes than elements. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
      integer(4) :: firstofset
      integer(4), intent(in) :: myrank, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: i, np

      ! Sanity checks
      if( myrank < 0 ) then
        write(*,*) "firstofset (Error): myrank < 0"
        call terminate
      end if
      if( set < 1 ) then
        write(*,*) "firstofset (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "firstofset (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "firstofset (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(np < myrank+1) then
        write(*,*) "firstofset (Error): np < myrank+1"
        call terminate
      end if

      ! Compute first element index on current rank
      firstofset = 1
      do i = 0, min( myrank-1, set-1)
        firstofset = firstofset + nofset(i, set, nprocs=np)
      end do
      if(set <= myrank) firstofset = 0

    end function firstofset
    !EOC

    !BOP
    ! !ROUTINE: lastofset
    ! !INTERFACE:
    function lastofset(myrank, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: myrank  ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! integer(4) :: nprocs  ! Number of processes in commuincator
    ! OUT:
    ! integer(4) :: lastofset  ! Index of the total set for the first index
    !                          ! of the current subset
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the index of the last element $j_\text{el}(p)$ a given process is responsible for.
    !   If there are more processes than elements, a process responsible for no element
    !   gets the assignment $j_\text{el}(p > N_\text{el}-1) = -1$.\\
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow j_\text{el}(0)=4, j_\text{el}(1)=7, j_\text{el}(2)=10$\\
    !   Example:\\
    !   $N_\text{el}=2, N_\text{p}=4 \rightarrow j_\text{el}(0)=1, j_\text{el}(1)=2, j_\text{el}(2)=-1, j_\text{el}(3)=-1$
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Changed behaviour if there are more processes than elements. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
      integer(4) :: lastofset
      integer(4), intent(in) :: myrank, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: i, np

      ! Sanity checks
      if( myrank < 0 ) then
        write(*,*) "lastofset (Error): myrank < 0"
        call terminate
      end if
      if( set < 1 ) then
        write(*,*) "lastofset (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "lastofset (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "lastofset (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(np < myrank+1) then
        write(*,*) "lastofset (Error): np < myrank+1"
        call terminate
      end if

      ! Compute last element index on this rank
      lastofset = 0
      do i = 0, min(myrank, set-1)
         lastofset = lastofset + nofset(i, set, nprocs=np)
      end do
      if(set <= myrank) lastofset = -1

    end function lastofset
    !EOC

    !BOP
    ! !ROUTINE: procofindex
    ! !INTERFACE:
    function procofindex(k, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: k    ! Element number k
    ! integer(4) :: set  ! Total number of distributed elements
    ! integer(4), optional :: nprocs ! Number of processes in communicator
    ! OUT:
    ! integer(4) :: procofindex  ! Rank that holds the element
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the index of the process $i_\text{p}(k)$ that is responsible for the element with index $k$.
    !   If $k$ is larger than $N_\text{el}$ or smaller than $1$ the routine returns terminates the execution.
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow i_\text{p}(1)=1, i_\text{p}(4)=1, i_\text{p}(5)=2, \dots $
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Adapted to changes in lastofset. (Aurich)
    !   Changed behaviour if k is smaller of larger than the set. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
      integer(4) :: procofindex
      integer(4), intent(in) :: k, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: iproc, np

      ! Sanity checks
      if( k < 1 ) then
        write(*,*) "procofindex (Error): k < 1"
        call terminate
      end if
      if( k > set ) then
        write(*,*) "procofindex (Error): set < k"
        call terminate
      end if
      if( set < 1 ) then
        write(*,*) "procofindex (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "procofindex (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "procofindex (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if

      ! Compute rank that holds element k
      procofindex = 0
      do iproc = 0, np - 1
        if(k > lastofset(iproc, set, nprocs=np)&
          & .and. lastofset(iproc, set, nprocs=np) > 0) procofindex = procofindex + 1
      end do

    end function procofindex
    !EOC

    !BOP
    ! !ROUTINE: lastproc
    ! !INTERFACE:
    function lastproc(col, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: col ! ``Column" of process grid (i.e. an element index of rank 0)
    ! integer(4) :: set ! Total number of distributed elements.
    ! integer(4), optional :: nprocs  ! Number of processes in communicator
    ! OUT:
    ! integer(4) :: lastproc ! Number of processes active in process column
    !
    ! !DESCRIPTION:
    !   This functions helps with collecting a set of $N_\text{el}$ elements which were
    !   distributed to $N_\text{p}$ {\tt MPI} processes in continuous blocks and is used for example for
    !   the writing of {\tt PMATXS.OUT}, {\tt EMAT.OUT}, {\tt SCCLI.OUT}
    !   and {\tt EXCLI.OUT}.
    !   For further describing the functionality of this routine, let us consider an example
    !   distribution: \\
    !   Let $N_\text{el} = 13$ and $N_\text{p} = 5$ : \\
    !   \begin{tabular}{c|ccc|c}
    !     rank & firstofset & \dots & lastofset & nofset \\
    !     \hline
    !     0 & 1 & 2 & 3 & 3 \\
    !     1 & 4 & 5 & 6 & 3 \\
    !     2 & 7 & 6 & 9 & 3 \\
    !     3 & 10 & 11 & - & 2 \\
    !     4 & 12 & 13 & - & 2 \\
    !   \end{tabular}
    !
    ! For inputs of $\text{col}=\{1,2,3\}, \text{set}=13$ the routine returns $\{4,4,2\}$, i.e. the process index
    ! of the last active process in the respective column. For all other input for col the routine halts execution.
    ! In the pathological case, where we have more processes than elements the
    ! following example depicts the routines behaviour:\\
    !   Let $N_\text{el} = 3$ and $N_\text{p} = 5$ : \\
    !   \begin{tabular}{c|ccc|c}
    !     rank & firstofset & \dots & lastofset & nofset \\
    !     \hline
    !     0 & 1 & 1 & - & 1 \\
    !     1 & 2 & 2 & - & 1 \\
    !     2 & 3 & 3 & - & 1 \\
    !     3 & 0 & -1 & - & 0 \\
    !     4 & 0 & -1 & - & 0 \\
    !   \end{tabular}
    ! For inputs of $\text{col}=1, \text{set}=3$ the routine returns $2$.
    ! For all other input for col execution is halted.
    ! In the other case pathological case, where we have only one processes the
    ! following example depicts the routines behaviour:\\
    !   Let $N_\text{el} = 3$ and $N_\text{p} = 1$ : \\
    !   \begin{tabular}{c|ccc|c}
    !     rank & firstofset & \dots & lastofset & nofset \\
    !     \hline
    !     0 & 1 & 2 & 3 & 3 \\
    !   \end{tabular}
    !
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Added sanity checks. (Aurich)
    !
    !EOP
    !BOC
      implicit none
      integer(4) :: lastproc
      integer(4), intent(in) :: col, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: np

      ! Sanity checks
      if( set < 1 ) then
        write(*,*) "lastproc (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "lastproc (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "lastproc (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(col > nofset(0, set, nprocs=np) .or. col < 1) then
        write(*,*) "lastproc (Error): col > nofset(0,set,np) or col < 1"
        call terminate
      end if

      ! Only in the last (or only) column less then all
      ! processes can be active.
      ! nofset(0,set,np) gives the number of columns.
      if(col /= nofset(0, set, nprocs=np)) then
        lastproc = np - 1
      ! Processes fit evenly (includes only one row
      else if(modulo(set,np) == 0) then
        lastproc = np - 1
      ! Dangling processes
      else
        ! Rest elements not filling last column
        lastproc = modulo(set, np) - 1
      end if

    end function lastproc
    !EOC

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! MPI wrapper for "convenience"             !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: barrier
    ! !INTERFACE:
    subroutine barrier(mpicom, callername)
    ! !DESCRIPTION:
    !   If -DMPI calls {\tt mpi\_barrier}, else nothing.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      implicit none
      type(mpiinfo), intent(in), optional :: mpicom
      character(*), intent(in), optional :: callername

      character(*), parameter :: thisname = "barrier"

      type(mpiinfo) :: mpinf

      if(present(mpicom)) then
        mpinf = mpicom
      else
        mpinf = mpiglobal
      end if

      if(present(callername)) then
        if(.false.) then
          write(*, '("Info(",a,"): Rank ",i3," of mpicom", i16," called barrier from ", a)')&
            & trim(thisname), mpinf%rank, mpinf%comm, trim(callername)
        end if
      end if

      ! do nothing if only one process
#ifndef MPI
      if(mpinf%procs .eq. 1) return
#endif
      ! call the mpi barrier
#ifdef MPI
      call mpi_barrier(mpinf%comm, mpinf%ierr)
      if(mpinf%ierr /= 0) then
        write(*,*) "Error (barrier): ierr =", mpinf%ierr
      end if
#endif
    end subroutine barrier
    !EOC

    !BOP
    ! !ROUTINE: mpi_allgatherv_ifc
    ! !INTERFACE:
    subroutine mpi_allgatherv_ifc(set, rlen, rlenv, ibuf, rlpbuf, rbuf, zbuf,&
      & inplace, comm)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: set              ! Number of elements in distributed set
    ! integer(4), optional :: rlen   ! Number of data elements per element (constant)
    ! integer(4), optional :: rlenv(set) ! Number of data elements per element
    ! logical, optional :: inplace   ! Use mpi_in_place
    ! type(mpiifo), optional :: comm ! MPI communicator type
    ! In/Out:
    ! integer(4), optional :: ibuf(*) ! Buffers to send/recive
    ! real(4), optional :: rlbuf(*)   ! for different data types
    ! real(8), optional :: rbuf(*)    !
    ! complex(8), optional :: zbuf(*) !
    !
    ! !DESCRIPTION:
    !   Wrapper routine for {\tt MPI\_ALLGATHERV} for different
    !   data types which is adapted for the k-point set
    !   distribution scheme. That is this works, if {\it set} number
    !   of elements (e.g. k-points) is distributed over all
    !   processes in the {\tt MPI} communicator {\it comm} using
    !   a continuous distribution as created by the functions
    !   {\tt nofset, firstofset, lastofset}.
    !   The routine can handle a constant number of data elements per
    !   set element by specifying {\tt rlen} or a set element dependent
    !   number of data elements by passing {\tt rlenv(set)}.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Added input parameter for communicator and
    !   a switch for inplace allgather. (Aurich)
    !   Added support for set element dependent number
    !   of data elements. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: set
      integer(4), intent(in), optional :: rlen
      integer(4), intent(in), optional :: rlenv(set)
      logical, intent(in), optional :: inplace
      type(mpiinfo), intent(in), optional :: comm
      integer(4), intent(inout), optional :: ibuf(*)
      real(4), intent(inout), optional :: rlpbuf(*)
      real(8), intent(inout), optional :: rbuf(*)
      complex(8), intent(inout), optional :: zbuf(*)

      ! Arrays for out of place send
      integer(4), allocatable :: bufi(:)
      real(4), allocatable :: bufrlp(:)
      real(8), allocatable :: bufr(:)
      complex(8), allocatable :: bufz(:)

      type(mpiinfo) :: mpicom
      integer(4) :: ierr
      integer(4), allocatable :: buf_n(:), buf_dspls(:)
      integer(4) :: j
      logical :: ti, tr, trlp, tz, tinplace
      integer(4) :: myrank, myprocs, mycomm

      ! Sanity checks
      if( set < 1 ) then
        write(*,*) "Error (mpi_allgatherv_ifc): set < 1"
        call terminate
      end if
      if(present(inplace)) then
        tinplace = inplace
      else
        tinplace = .false.
      end if
      if(present(comm)) then
        mpicom = comm
      else
        mpicom = mpiglobal
      end if
      ti = present(ibuf)
      tr = present(rbuf)
      trlp = present(rlpbuf)
      tz = present(zbuf)
      if(count((/ti, tr, trlp, tz/)).ne.1) then
        write(*,*)
        write(*,'("Error (mpi_allgatherv_ifc): Exactly one array must be defined.")')
        write(*,*)
        call terminate
      end if
      if(present(rlen) .and. present(rlenv)&
        & .or. .not. present(rlen) .and. .not. present(rlenv)) then
        write(*,*)
        write(*,'("Error (mpi_allgatherv_ifc): Specifiy either rlen or rlenv")')
        write(*,*)
        call terminate
      end if
      if(present(rlen)) then
        if(rlen < 0) then
          write(*,'("Error (mpi_allgatherv_ifc): rlen < 0")')
          call terminate
        end if
      end if
      if(present(rlenv)) then
        if(any(rlenv < 0)) then
          write(*,'("Error (mpi_allgatherv_ifc): rlenv < 0")')
          call terminate
        end if
      end if

      myrank = mpicom%rank
      myprocs = mpicom%procs
      mycomm = mpicom%comm

#ifdef MPI
      allocate(buf_n(myprocs), buf_dspls(myprocs))

      ! Number of elements in send buffer (flattened array)
      if(present(rlen)) then
        buf_n =(/(rlen*nofset(j, set, myprocs), j = 0, myprocs-1)/)
      else
        do j = 0, myprocs-1
          buf_n(j+1) = sum(rlenv(firstofset(j,set,myprocs):lastofset(j,set,myprocs)))
        end do
      end if

      ! Displacements within receive buffer (flattened array)
      if(present(rlen)) then
        buf_dspls =(/(rlen*(firstofset(j, set, myprocs)-1), j = 0, myprocs-1)/)
      else
        do j = 0, myprocs-1
          buf_dspls(j+1) = sum(buf_n(1:j))
        end do
      end if

      ! Integers
      if(ti) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufi(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufi(:)= ibuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufi, &
            buf_n(myrank+1), &
            mpi_integer, &
            ibuf, &
            buf_n, &
            buf_dspls, &
            mpi_integer, &
            mycomm, &
            ierr)
          deallocate(bufi)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_integer, &
            ibuf, &
            buf_n, &
            buf_dspls, &
            mpi_integer, &
            mycomm, &
            ierr)
        end if

      end if

      ! Floats
      if(trlp) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufrlp(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufrlp(:) =&
              & rlpbuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufrlp, &
            buf_n(myrank+1), &
            mpi_real4, &
            rlpbuf, &
            buf_n, &
            buf_dspls, &
            mpi_real4, &
            mycomm, &
            ierr)
          deallocate(bufrlp)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_real4, &
            rlpbuf, &
            buf_n, &
            buf_dspls, &
            mpi_real4, &
            mycomm, &
            ierr)
        end if

      end if

      ! Doubles
      if(tr) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufr(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufr(:)= rbuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufr, &
            buf_n(myrank+1), &
            mpi_double_precision, &
            rbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_precision, &
            mycomm, &
            ierr)
          deallocate(bufr)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_double_precision, &
            rbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_precision, &
            mycomm, &
            ierr)
        end if

      end if

      ! Complex doubles
      if(tz) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufz(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufz(:)= zbuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufz, &
            buf_n(myrank+1), &
            mpi_double_complex, &
            zbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_complex, &
            mycomm, &
            ierr)
          deallocate(bufz)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_double_complex, &
            zbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_complex, &
            mycomm, &
            ierr)
        end if

      end if

      deallocate(buf_n, buf_dspls)

#endif
    end subroutine
    !EOC

    !++++++++++++++++++++++++++++++++++++++++++++!
    ! Setup routine for distributing N elements  !
    ! to P processes where each element may have !
    ! more than one process associated to it.    !
    !++++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: setup_proc_groups
    ! !INTERFACE:
    subroutine setup_proc_groups(ngroups, mygroup)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! input(4) :: ngroups ! number of groups to be formed
    !                     ! form the available MPI threads
    ! OUT:
    ! type(procgroup) :: mygroup ! Type containing MPI information
    !
    ! !DESCRIPTION:
    !   Given a total of $N_\text{p}$ {\tt MPI} processes
    !   this routine generates $N_\text{g}$ process groups
    !   of size $N_\text{p}/N_\text{g}$ (no accutal {\tt MPI}
    !   groups, but rather full {\tt MPI} communicators).
    !   $N_\text{p} \geq N_\text{g}$ is required.
    !   If $\text{mod}(N_\text{p}, N_\text{g}) \neq 0$ then
    !   one more group is created with the rest of the processes.
    !   Also a commuicator between the first processes in each such
    !   group is created.
    !
    ! !REVISION HISTORY:
    !   Based of setup_ProcGroups form mpi_mortadella branch. (Aurich)
    !
    !EOP
    !BOC

      implicit none

      integer(4), intent(in) :: ngroups
      type(procgroup), intent(out) :: mygroup

      integer(4) :: color, key, i
      integer(4) :: global_group, interprocs_group
      integer(4) :: ngprocs, dangling_procs
      integer(4), allocatable :: proclist(:)

#ifdef MPI

      ! Sanity check
      if(ngroups > mpiglobal%procs) then
        write(*,*) "setup_proc_groups (ERROR): ngroups > procs"
        call terminate
      end if

      ! Integer part of Np/Ng
      ! Ex: Np = 10, Ng =3 --> Ngp = 3
      ngprocs = mpiglobal%procs/ngroups
      ! Rest part of Np/Ng
      ! Ex: Np = 10, Ng =3 --> Ndp = 1
      dangling_procs = mod(mpiglobal%procs, ngroups)
      ! Resulting number of process groups
      if(dangling_procs > 0) then
        mygroup%ngroups = ngroups+1
      else
        mygroup%ngroups = ngroups
      end if

      ! Group number
      ! Ex: Np = 10, Ng = 3 --> 3+1 groups
      ! world rank p = 0,1,2,3,4,5,6,7,8,9
      !       Gid(p) = 0,0,0,1,1,1,2,2,2,3
      mygroup%id = mpiglobal%rank/ngprocs

      ! Dividing the processes
      ! Ex: Np = 10, Ng = 3 --> 3+1 groups
      !     Group com/id  size  ranks  global ranks
      !             1/0    3    0,1,2   0,1,2
      !             2/1    3    0,1,2   3,4,5
      !             3/2    3    0,1,2   6,7,8
      !             4/3    1    0       9
      color = mygroup%id
      key = mpiglobal%rank
      ! Split global comm intra-group comms
      call mpi_comm_split(mpiglobal%comm, color, key,&
        & mygroup%mpi%comm, mygroup%mpi%ierr)

write (*, '("setup_proc_groups@rank",i3,"mycolor=", i3," mygroup%mpi%comm=",i16)')&
  & mpiglobal%rank, color, mygroup%mpi%comm

      !   Error checking
      if(mygroup%mpi%ierr .ne. 0) then
        write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mygroup%mpi%ierr, "com split failed."
        call terminate
      end if
      ! Get number of procs in group comm
      call mpi_comm_size(mygroup%mpi%comm, mygroup%mpi%procs, mygroup%mpi%ierr)
      !   Error checking
      if(mygroup%mpi%ierr .ne. 0) then
        write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mygroup%mpi%ierr, "comm size failed."
        call terminate
      end if
      ! Get rank
      call mpi_comm_rank(mygroup%mpi%comm, mygroup%mpi%rank,  mygroup%mpi%ierr)
      !   Error checking
      if(mygroup%mpi%ierr .ne. 0) then
        write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mygroup%mpi%ierr, "comm rank failed."
        call terminate
      end if

      ! Inter-groups communicator
      ! Ex: Np = 10, Ng = 3 --> 3+1 groups
      !     intercom  size  ranks    global ranks
      !         5      4    0,1,2,3  0,3,6,9

      ! Global ranks of first group ranks
      allocate(proclist(ngroups))
      proclist = (/(i*ngprocs, i=0,ngroups-1)/)

      ! Make MPI group from global MPI communicator, i.e. a group of all processes
      call mpi_comm_group(mpiglobal%comm, global_group, mpiglobal%ierr)
      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "comm group failed."
        call terminate
      end if
      ! From that global group make a new sub group only containing the ranks in proclist
      call mpi_group_incl(global_group, ngroups, proclist, interprocs_group,&
        & mpiglobal%ierr)
      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "group incl failed."
        call terminate
      end if

      ! Inter-groups communicator
      mygroup%mpiintercom%comm = mpi_comm_null
      mygroup%mpiintercom%rank = -1
      mygroup%mpiintercom%procs = ngroups

      ! Take that sub group and create a new communicator for only that sub group
      call mpi_comm_create(mpiglobal%comm, interprocs_group, mygroup%mpiintercom%comm,&
        & mpiglobal%ierr)
      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "comm create failed."
        call terminate
      end if
      if(mygroup%mpiintercom%comm /= mpi_comm_null) then
        ! Collect information for intercommunication communicator
        call mpi_comm_size(mygroup%mpiintercom%comm, mygroup%mpiintercom%procs,&
          & mygroup%mpiintercom%ierr)
        !   Error checking
        if(mpiglobal%ierr .ne. 0) then
          write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
            & mpiglobal%rank, mpiglobal%ierr, "comm create failed."
          call terminate
        end if
        call mpi_comm_rank(mygroup%mpiintercom%comm, mygroup%mpiintercom%rank,&
          & mygroup%mpiintercom%ierr)
        !   Error checking
        if(mpiglobal%ierr .ne. 0) then
          write (*, '("setup_proc_groups@rank",i3," (ERROR ",i3,"):",a)')&
            & mpiglobal%rank, mpiglobal%ierr, "comm create failed."
          call terminate
        end if
      end if

      ! Cleanup
      deallocate(proclist)
      ! Free group handles, the communicators are all we need form here on
      call mpi_group_free(global_group, mpiglobal%ierr)
      call mpi_group_free(interprocs_group, mpiglobal%ierr)

#endif

    end subroutine setup_proc_groups
    !EOC

    !++++++++++++++++++++++++++++++++++++++!
    ! Primitive node information gathering !
    ! Sets up intra- and inter-node        !
    ! communicators.                       !
    !++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: setup_node_groups
    ! !INTERFACE:
    subroutine setup_node_groups
    ! !INPUT/OUTPUT PARAMETERS:
    ! Module OUT:
    ! type(procgroup) :: mpinodes ! Partitioning according to processor name
    !
    ! !DESCRIPTION:
    !   WARNING: This can only ever work if you pin the {\tt MPI} threads
    !   to the processors during program execution.\\
    !   Given a total of $N_\text{p}$ {\tt MPI} processes
    !   this routine generates $N_\text{nodes}$ process groups.
    !   The size of each process group depends on the respective node size
    !   and node utilization.
    !   Also a communicator between the first processes in each node is created.
    !   A node in the context of this routine is a processor
    !   (return value of {\tt mpi\_get\_processor\_name}).
    !
    ! !REVISION HISTORY:
    !   Based of setup_ProcGroups form mpi_mortadella branch
    !   and get_isfirstinnode form master. (Aurich)
    !
    !EOP
    !BOC

      implicit none

      ! Vars for splitting the MPI com
      integer(4) :: color, key, i
      integer(4) :: global_group, interprocs_group

      ! Vars for determining "node" layout
      integer(4) :: strsize
      integer(4) :: procnamelen
      integer(4) :: pos1, pos2, n
#ifdef MPI
      integer(4) :: recvstatus(mpi_status_size)
#endif
      character(200) :: procname, myprocname
      character(len=200*mpiglobal%procs) :: neighbors, neighborssend
      logical :: lbuffer
      integer(4) :: mynodesize, mynode, nnodes
      integer(4), allocatable :: mynoderanks(:)
      integer(4), allocatable :: nodechefs(:)
      integer(4), parameter :: tag = 25

#ifdef MPI

      !++++++++++++++++++++++++++++++++++++++++++!
      ! Stupidly simplistic hardware inspection. !
      ! One should do this with hwloc.           !
      !++++++++++++++++++++++++++++++++++++++++++!

      ! Assume maximal processor name length to be 200
      strsize = 200*mpiglobal%procs

      myprocname = ''
      neighbors = ''
      neighborssend = ''

      ! Our definition of "node" is a processor.
      call mpi_get_processor_name(myprocname, procnamelen, mpiglobal%ierr)
      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "get processor name failed."
        call terminate
      end if

      ! The following code is illustrated by the example
      ! with 4 threads and 2 physical processors A and B.
      !   rank 0 has an empty neighbors string (ns) "" and a procname (pn) of "A",
      !   it sends ",A" to rank 1
      !   rank 1 now has ns=",A" and pn="A", it sends ",A,A" to rank 2
      !   rank 2 now has ns=",A,A" and pn="B", it sends ",A,A,B" to rank 3
      !   rank 3 now has ns=",A,A,B" and pn="B", it has no one to send to.
      !   Each rank checks if pn is in ns, if not it is the nodechef.
      if(mpiglobal%rank > 0) then
        call mpi_recv(neighbors, strsize, mpi_character,&
          & mpiglobal%rank-1, tag, mpiglobal%comm, recvstatus, mpiglobal%ierr)
      end if
      if(mpiglobal%rank < mpiglobal%procs-1) then
        write(neighborssend,*) trim(adjustl(neighbors))//","//trim(adjustl(myprocname))
        call mpi_send(neighborssend, strsize, mpi_character,&
          & mpiglobal%rank+1, tag, mpiglobal%comm, mpiglobal%ierr)
      endif
      if(index(trim(adjustl(neighbors)), trim(adjustl(myprocname)) ) > 0) then
        firstinnode = .false.
      else
        firstinnode = .true.
      endif

      ! Send the full list of processors to all
      if(mpiglobal%rank == mpiglobal%procs-1) then
        write(neighborssend,*) trim(adjustl(neighbors))//","//trim(adjustl(myprocname))
        ! Cut leading comma
        neighborssend = trim(adjustl(neighborssend))
        neighborssend = neighborssend(2:)
      end if
      call mpi_bcast(neighborssend, strsize, mpi_character,&
        mpiglobal%procs-1, mpiglobal%comm, mpiglobal%ierr)

      ! Figure out how many MPI threads are on current processor
      mynodesize = countsubstring(trim(adjustl(neighborssend)),&
        & trim(adjustl(myprocname)))

      ! Make list of ranks on current processor
      allocate(mynoderanks(mynodesize))
      procname=''
      pos1=1
      n=0
      do i = 0, mpiglobal%procs-1
        pos2 = index(neighborssend(pos1:),",")
        if(pos2 == 0) then
          procname = neighborssend(pos1:)
          if( trim(adjustl(myprocname)) == trim(adjustl(procname)) ) then
            n = n+1
            mynoderanks(n) = i
          end if
          exit
        end if
        procname = neighborssend(pos1:pos1+pos2-2)
        if( trim(adjustl(myprocname)) == trim(adjustl(procname)) ) then
          n = n+1
          mynoderanks(n) = i
        end if
        pos1=pos1+pos2
      end do

      ! The following collects all the rank numbers that are nodechefs
      ! in the first 1:nnodes elements of nodechefs and creates
      ! a mpi group and communicator for inter-node communication.
      allocate(nodechefs(mpiglobal%procs))
      nnodes = 0
      do i = 0, mpiglobal%procs-1
        lbuffer = firstinnode
        call mpi_bcast(lbuffer, 1, mpi_logical, i, mpiglobal%comm, mpiglobal%ierr)
        if(lbuffer)then
           nnodes = nnodes+1
           nodechefs(nnodes)= i
        endif
      end do

      ! Figure out on which node the current process is
      mynode = -1
      do i = 1, nnodes
        if(mynoderanks(1) == nodechefs(i)) then
          exit
        end if
      end do
      mynode = i-1

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! We somehow determined how many nodes we are using,     !
      ! how may and which MPI threads are running on each node !
      ! and determined a master thread for each node.          !
      ! Now we make MPI communicators for each node and one    !
      ! communicator that facilitates inter-node communication.!
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! nnodes, mynode, mynoderanks, nodechefs

      mpinodes%ngroups = nnodes
      mpinodes%id = mynode

      color = mpinodes%id
      key = mpiglobal%rank

      ! Split global comm intra-group comms
      call mpi_comm_split(mpiglobal%comm, color, key,&
        & mpinodes%mpi%comm, mpinodes%mpi%ierr)

      !   Error checking
      if(mpinodes%mpi%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpinodes%mpi%ierr, "com split failed."
        call terminate
      end if

      ! Get number of procs in group comm
      call mpi_comm_size(mpinodes%mpi%comm, mpinodes%mpi%procs, mpinodes%mpi%ierr)

      !   Error checking
      if(mpinodes%mpi%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpinodes%mpi%ierr, "comm size failed."
        call terminate
      end if

      ! Get rank
      call mpi_comm_rank(mpinodes%mpi%comm, mpinodes%mpi%rank,  mpinodes%mpi%ierr)
      call barrier(mpinodes%mpi)

      !   Error checking
      if(mpinodes%mpi%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpinodes%mpi%ierr, "comm rank failed."
        call terminate
      end if

      ! Inter-groups communicator
      mpinodes%mpiintercom%comm = mpi_comm_null
      mpinodes%mpiintercom%rank = -1
      mpinodes%mpiintercom%procs = nnodes

      ! Make MPI group from global MPI communicator, i.e. a group of all processes
      call mpi_comm_group(mpiglobal%comm, global_group, mpiglobal%ierr)

      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "comm group failed."
        call terminate
      end if

      ! From that global group make a new sub group only containing
      ! the ranks in nodechefs
      call mpi_group_incl(global_group, nnodes, nodechefs(1:nnodes),&
        & interprocs_group, mpiglobal%ierr)

      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "group incl failed."
        call terminate
      end if

      ! Take that sub group and create a new commuicator for only that sub group
      call mpi_comm_create(mpiglobal%comm, interprocs_group,&
        & mpinodes%mpiintercom%comm, mpiglobal%ierr)


      !   Error checking
      if(mpiglobal%ierr .ne. 0) then
        write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
          & mpiglobal%rank, mpiglobal%ierr, "comm create failed."
        call terminate
      end if

      ! Collect information for intercommunication communicator
      !   Note: Not only the ranks 0 in the nodes have a valid (non null)
      !         intercommunicator
      if(mpinodes%mpiintercom%comm /= mpi_comm_null) then

        call mpi_comm_size(mpinodes%mpiintercom%comm, mpinodes%mpiintercom%procs,&
          & mpinodes%mpiintercom%ierr)

        !   Error checking
        if(mpiglobal%ierr .ne. 0) then
          write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
            & mpiglobal%rank, mpiglobal%ierr, "comm create failed."
          call terminate
        end if

        call mpi_comm_rank(mpinodes%mpiintercom%comm, mpinodes%mpiintercom%rank,&
          & mpinodes%mpiintercom%ierr)

        !   Error checking
        if(mpiglobal%ierr .ne. 0) then
          write (*, '("setup_node_groups@rank",i3," (ERROR ",i3,"):",a)')&
            & mpiglobal%rank, mpiglobal%ierr, "comm create failed."
          call terminate
        end if

      end if

      ! Free group handles, the communicators are all we need form here on
      call mpi_group_free(global_group, mpiglobal%ierr)
      call mpi_group_free(interprocs_group, mpiglobal%ierr)

      ! Set to global variables
      firstinnode_comm = mpinodes%mpiintercom%comm

#endif
      contains

        function countsubstring(s1, s2) result(c)
          character(*), intent(in) :: s1, s2
          integer :: c, p, posn

          c = 0
          if(len(s2) == 0) return
          p = 1
          do
            posn = index(s1(p:), s2)
            if(posn == 0) return
            c = c + 1
            p = p + posn + len(s2)
          end do
        end function

    end subroutine setup_node_groups
    !EOC

    !++++++++++++++++++++++++++++++!
    ! Older code should be removed !
    !++++++++++++++++++++++++++++++!

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Primitive node information gathering (Old version)   !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!    !BOP
!    ! !ROUTINE: get_isfirstinnode
!    ! !INTERFACE:
!    subroutine get_isfirstinnode(strsize)
!    ! !INPUT/OUTPUT PARAMETERS:
!    ! IN:
!    ! integer(4) :: strsize ! Size of neighbour string
!    !
!    ! !DESCRIPTION:
!    ! The routine uses the processor name returned by
!    ! {\tt mpi\_get\_processor\_name} to determine the first rank
!    ! on each node (nodechef). For those ranks a mpi group
!    ! and communicator is created. Also sets the firstinnode
!    ! flag in the modmpi module.
!    !
!    ! !REVISION HISTORY:
!    ! Added to documentation scheme. (Aurich)
!    !EOP
!    !BOC
!
!      integer(4), intent(in)::strsize
!#ifdef MPI
!      integer(4) :: ierr
!      integer(4) :: recvstatus(mpi_status_size), i
!      character(len = strsize) :: neighbors, neighborssend
!      character(200) :: procname
!      logical :: lbuffer
!      integer(4) :: procnamelen
!      integer(4) :: nnodes
!      integer(4) :: world_group, firstinnode_group
!      integer(4), allocatable :: nodechefs(:)
!      integer(4), parameter :: tag = 25
!
!      allocate(nodechefs(procs))
!
!      procname = ''
!      neighbors = ''
!      neighborssend = ''
!
!      call mpi_get_processor_name(procname, procnamelen, ierr)
!
!      ! The following code is illustrated by the example
!      ! with 4 threads and 2 physical processors A and B.
!      !   rank 0 has an empty neighbors string (ns) "" and a procname (pn) of "A",
!      !   it sends ",A" to rank 1
!      !   rank 1 now has ns=",A" and pn="A", it sends ",A,A" to rank 2
!      !   rank 2 now has ns=",A,A" and pn="B", it sends ",A,A,B" to rank 3
!      !   rank 3 now has ns=",A,A,B" and pn="B", it has no one to send to.
!      !   Each rank checks if pn is in ns, if not it is the nodechef.
!      if(rank .gt. 0) then
!        call mpi_recv(neighbors, strsize, mpi_character,&
!          & rank-1, tag, mpi_comm_world, recvstatus, ierr)
!      end if
!      if(rank .lt. procs-1) then
!        write(neighborssend,*) adjustl(trim(neighbors))//","//adjustl(trim(procname))
!        call mpi_send(neighborssend, strsize, mpi_character,&
!          & rank+1, tag, mpi_comm_world, ierr)
!      endif
!      if(index(neighbors, adjustl(trim(procname))).gt.0) then
!        firstinnode = .false.
!      else
!        firstinnode = .true.
!      endif
!
!      ! The following collects all the rank numbes that are nodechefs
!      ! in the first 1:nnodes elements of nodechefs and chreates
!      ! a mpi group and communicator for internode communication.
!      nnodes = 0
!      do i = 0, procs-1
!        lbuffer = firstinnode
!        call mpi_bcast(lbuffer, 1, mpi_logical, i, mpi_comm_world, ierr)
!        if(lbuffer)then
!           nnodes = nnodes+1
!           nodechefs(nnodes)= i
!        endif
!      end do
!      call mpi_comm_group(mpi_comm_world, world_group, ierr)
!      call mpi_group_incl(world_group, nnodes, nodechefs, firstinnode_group, ierr)
!      call mpi_comm_create(mpi_comm_world, firstinnode_group, firstinnode_comm, ierr)
!
!      deallocate(nodechefs)
!      call mpi_group_free(world_group)
!      call mpi_group_free(firstinnode_group)
!
!#endif
!    end subroutine
!    !EOC

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Older k-point partitioning routines.      !
    !+++++++++++++++++++++++++++++++++++++++++++!

    ! Service functions to partition k points
    ! still kept around but should be replaced by generig functions
    ! firstofset nofset lastofset ....

    function nofk(process)
      use mod_kpoint, only: nkpt
      integer(4) :: nofk
      integer(4), intent(in) :: process
      nofk = nkpt / procs
      if((mod(nkpt, procs) .gt. process)) nofk = nofk + 1
    end function nofk

    function firstk(process)
       integer(4) :: firstk
       integer(4), intent(in) :: process
       integer(4)::i
       firstk = 1
       do i = 0, process - 1
          firstk = firstk + nofk(i)
       end do
    end function firstk

    function lastk(process)
       integer(4) :: lastk, i
       integer(4), intent(in) :: process
       lastk = 0
       do i = 0, process
          lastk = lastk + nofk(i)
       end do
    end function lastk

    function procofk(k)
       integer(4) :: procofk
       integer(4), intent(in) :: k
       integer(4) :: iproc
       procofk = 0
       do iproc = 0, procs - 1
          if(k .gt. lastk(iproc)) procofk = procofk + 1
       end do
    end function procofk

end module modmpi
