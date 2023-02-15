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

! TODO(Alex) Issue #23.
! * MPI wrappers are a mix of types and global variables, the latter of which isn't
!   required as one can always query MPI variables with MPI calls.
! * There is no overload set for compilation in serial.
! * Routines don't follow single responsibility.
! * Many routines use if statements to determine type: This should be done at compile-time with overloads
!   and a generic interface.
! * Preprocessor variables are not obfuscated.

! This needs refactoring if we hope to use more complex MPI distribution schemes.
! The approach should be to gradually migrate the routines out of this module and into routines/
! with serial overloads in serial/

!> Main exciting MPI module, in the process of being depreciated
module modmpi
  use trace, only: trace_back
#ifdef MPI
  use mpi
#endif
#ifdef SIRIUS
  ! Have to call the external lib rather than the API
  ! to avoid a circular dependency (as this module is far too large)
  use sirius, only: sirius_initialize, sirius_finalize, sirius_print_timers
#endif
  use exciting_mpi, only: mpiinfo
  implicit none

  ! TODO(Alex) Issue #23. Type for refactoring
  !> Groups of MPI communicators connected via inter-communicator
  type procgroup
    !> Total number of process groups this group belongs to
    integer(4) :: ngroups
    !> Group id
    integer(4) :: id
    !> MPI information for current process group
    type(mpiinfo) :: mpi
    !> Inter-groups communicator
    type(mpiinfo) :: mpiintercom
  end type procgroup


  ! TODO(Alex) Issue #23. Legacy data for refactoring
  !> mpiinfo with global scope. Instantiating a global instance of the MPI env is pointless
  type(mpiinfo) :: mpiglobal
  !> MPI communicator which is used to distribute k-points.
  type(mpiinfo) :: mpi_env_k
  !> MPI communicator which is used to parallelized bands and diagonalization in SIRIUS.
  !> In pure Exciting this communicator is set to the trivial MPI_COMM_SELF
  type(mpiinfo) :: mpi_env_band

  !> Nodes as procgroup
  type(procgroup) :: mpinodes
  !> Variables (contained equivalently contained in mpiglobal)
  integer(4) :: rank
  integer(4) :: procs
  integer(4) :: ierr

  !> Variables (contained in mpinodes)
  integer(4) :: firstinnode_comm
  !> Some parts use these
  logical :: splittfile, firstinnode


contains

   !> Initialise a global instance of the MPI environment, mpiglobal.
   !>
   !> Also initialises global legacy variables: procs, rank. These should be scrapped.
  subroutine initmpi()
    integer :: ierr
    !> MPI communicator
    integer :: comm
    comm = 0

#ifdef MPI
#ifdef SIRIUS
    call sirius_initialize(.true.)
#else
    call mpi_init(ierr)
#endif
    comm = mpi_comm_world
#endif

    call mpiglobal%init(comm)

   ! Set legacy globals, because who knows where they're used
    procs = mpiglobal%procs
    rank = mpiglobal%rank
    ierr = mpiglobal%ierr

    ! Make communicators for intra-and inter-processor (node) communication.
    call setup_node_groups()
    ! Each rank writes its own file (use in parts of GS and XS)
    splittfile = .true.

    ! TODO(Alex) Issue 26. How should this behave for processes = 1?
#ifndef MPI
    splittfile = .false.
    firstinnode = .true.
#endif
 end subroutine initmpi


  !> Finalize the MPI environment
  !> Routine waits for all communicators: mpiglobal, mpinodes%mpi and mpinodes%mpiintercom,
  !> before finalizing.
  !>
  !> Note: For other communicators you need to make sure they have finished communication.
  subroutine finitmpi()
    integer(4) :: ierr = 0
    logical :: flag

    ! Wait for everyone to reach this point
    call barrier(mpiglobal)
#ifdef MPI
    !call barrier(mpinodes%mpi)
    !if (mpinodes%mpiintercom%rank >= 0) then
    !   call barrier(mpinodes%mpiintercom)
    !end if
#ifdef SIRIUS
    call sirius_finalize(call_mpi_fin=.true.)
    if ( mpiglobal%rank == 0 ) then
      call sirius_print_timers( .false. )
    endif
#else
    call mpi_finalize(ierr)
    if (ierr /= 0) then
       write (*, *) "Error (finitmpi): ierr =", ierr
    end if
#endif
#endif
 end subroutine finitmpi

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

  !> Terminate exciting if condition is false
  subroutine terminate_if_false(condition, message)
    use iso_fortran_env, only: error_unit
    !> Error condition, need to be false to terminate the program and print the message
    logical, intent(in) :: condition
    !> Error message that is printed to the terminal if present and condition is false.
    character(*), optional, intent(in) :: message

    character(512) :: error_message

    error_message = 'Error'
    if (present(message)) then
      error_message = trim(error_message)//': '//trim(adjustl(message))
    end if

    if(.not. condition) then
      write(error_unit, *)
      write(error_unit, '(a)') trim(error_message)
      write(error_unit, *)
      call trace_back()
      call terminate()
    end if
  end subroutine terminate_if_false

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

    !> Resolve MPI error code to error message.
    !> Appends `errstring` with `errmsg` and the resolved MPI error message
    !> corresponding to error code `ierr`.
    subroutine mpi_error_to_string( errmsg, ierr, errstring)
      use precision, only: str_1024
      !> custom error message to add
      character(*), intent(in) :: errmsg
      !> MPI error code
      integer, intent(in) :: ierr
      !> error string to be appended
      character(:), allocatable, intent(inout) :: errstring

      integer :: errlen, i
      character(len=str_1024) :: mpi_error

      if (ierr == 0) return

#ifdef MPI
      call MPI_error_string( ierr, mpi_error, errlen, i)
      errstring = trim(errstring)//trim(errmsg)//" (MPI Error: "//trim(mpi_error)//")"//NEW_LINE('a')
#else
      errstring = trim(errstring)//trim(errmsg)//" (No MPI used.)"//new_line('a')
#endif
    end subroutine mpi_error_to_string


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

    !> Assuming we have a bunch of work that can be split into \(N_p\) parts that can
    !> be executed independently of another and in parallel. Each part amounts to a
    !> computational load of \(L_p\) units of work, where each unit of work represents 
    !> a piece of work that can be executed in parallel and all \(L_p\) units are assumed 
    !> to take approximately the same amount of time to execute (e.g. k-points).
    !> This function distributes the given parts among the number of available processes 
    !> and uses further an internal parallelization over units of work inside of each part.
    !> The resulting execution schedule is returned as a 2D integer matrix, where each
    !> row corresponds to one process and each column to a unit of work. The value of elements
    !> in the matrix give the index of a part and a value of 0 means that nothing will
    !> be executed at this point in the schedule.
    !>
    !> **Example:**
    !> We have `np=7` parts with loads `load=[30, 20, 20, 40, 32, 32, 16]` and
    !> want to distribute them among `num_procs=11` processes. Further, we specify that
    !> at least `min_procs_per_part=3` and at most `max_procs_per_part=7` processes should work
    !> on one part.
    !> The resulting schedule will look like (0s not displayed)
    !>
    !>```
    !> 0| 2 2 2 3 3 3 5 5 5 5 5 1 1 1 1 1
    !> 1| 2 2 2 3 3 3 5 5 5 5 5 1 1 1 1 1
    !> 2| 2 2 2 3 3 3 5 5 5 5 5 1 1 1 1 1
    !> 3| 2 2 2 3 3 3 5 5 5 5 5 1 1 1 1 1
    !> 4| 2 2 2 3 3 3 5 5 5 5 5 1 1 1 1 1
    !> 5| 2 2 2 3 3 3 5 5 5 5 5 1 1 1 1 1
    !> 6| 2 2 2 3 3 3 5 5 5 5 5   4 4 4 4 4 4 4 4
    !> 7| 7 7 7 7 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4
    !> 8| 7 7 7 7 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4
    !> 9| 7 7 7 7 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4
    !>10| 7 7 7 7 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4
    !>```
    !>
    !> This means that part 7 will be executed by processes 7 - 10 and its work load of 16 units
    !> will be further separated into 4 blocks of 4 among the participating processes.
    !> The total load of 32 units of part 5 will be separated into 7 blocks of 5 which will be 
    !> executed by processes 0 - 6.
    !>
    !> Another way to view the schedule is to say that, e.g., process 6 will work on parts 2, 3, 5 and 4
    !> in this order.
    function gen_schedule( num_procs, np, load, min_procs_per_part, max_procs_per_part ) result( schedule )
      use sorting
      !> number of available processes
      integer, intent(in) :: num_procs
      !> number of parts to distribute
      integer, intent(in) :: np
      !> computational load per part
      integer, intent(in) :: load(:)
      !> minimum number of processes per part
      integer, intent(in) :: min_procs_per_part
      !> maximum number of processes per part
      integer, intent(in) :: max_procs_per_part
      !> schedule
      integer, allocatable :: schedule(:,:)

      integer :: i, ip, icol, iblock, tot_load

      integer, allocatable :: procs(:,:), excess(:,:), sort(:,:), order(:), blocks(:,:), res(:,:)
      logical, allocatable :: done(:)

      tot_load = sum( load )

      allocate( procs(tot_load, np), excess(tot_load, np) )
      
      do ip = 1, np
        do i = 1, tot_load
          procs(i, ip) = ceiling( dble( load(ip) ) / i )
          if( procs(i, ip) > max_procs_per_part .or. procs(i ,ip) < min_procs_per_part ) &
            procs(i, ip) = 0
          excess(i, ip) = procs(i, ip) * i - load(ip)
        end do
      end do

      allocate( res(num_procs, tot_load), source=0 )
      allocate( done(np), source=.false. )
      allocate( sort(3, np), order(np) )

      icol = 0
      do while( .not. all( done ) )
        icol = icol + 1
        ! find blocks
        blocks = find_blocks( res(:, 1:icol) )
        iblock = 0
        do while( iblock < size( blocks, dim=2 ) )
          iblock = iblock + 1
          sort(1, :) = -procs(blocks(4, iblock), :)
          sort(2, :) = excess(blocks(4, iblock), :)
          sort(3, :) = [(ip, ip=1, np)]
          order = sort_index_2d( 3, np, sort, 3 )
          do i = 1, np
            ip = order(i)
            if( done(ip) .or. procs(blocks(4, iblock), ip) == 0 .or. &
                procs(blocks(4, iblock), ip) > blocks(3, iblock) ) cycle
            ! add part to block
            res(blocks(1, iblock):blocks(1, iblock)+procs(blocks(4, iblock), ip)-1, &
                blocks(2, iblock):blocks(2, iblock)+blocks(4, iblock)-1) = ip
            done(ip) = .true.
            ! update blocks
            blocks = find_blocks( res(:, 1:icol) )
            iblock = 0
            exit
          end do
        end do
      end do

      if( allocated( schedule ) ) deallocate( schedule )
      allocate( schedule(num_procs, icol) )
      schedule = res(:, 1:icol)

      deallocate( procs, excess, res, done, sort, order )
      if( allocated( blocks ) ) deallocate( blocks )

      contains
        
        function find_blocks( schedule ) result( blocks )
          integer, intent(in) :: schedule(:,:)
          integer, allocatable :: blocks(:,:)

          integer, parameter :: num_max_blocks = 1000

          integer :: nrow, ncol, nblock, irow, icol, iblock, width, height

          integer, allocatable :: tmp(:,:), sort(:,:), order(:)

          nrow = size( schedule, dim=1 )
          ncol = size( schedule, dim=2 )
          allocate( tmp(4, num_max_blocks), source=0 )

          nblock = 0
          do icol = 1, ncol
            irow = 1
            do while( irow <= nrow ) 
              if( schedule(irow, icol) /= 0 ) then
                irow = irow + 1
                cycle
              end if
              ! find block height
              height = 0
              do while( schedule(irow+height, icol) == 0 )
                height = height + 1
                if( irow + height > nrow ) exit
              end do
              ! find block width
              width = 0
              do while( all( schedule(irow:irow+height-1, icol+width) == 0 ) )
                width = width + 1
                if( icol + width > ncol ) exit
              end do
              ! check if block is part of another block
              do iblock = 1, nblock
                if( tmp(1, iblock) <= irow .and. tmp(2, iblock) <= icol .and. &
                    tmp(3, iblock) >= height .and. tmp(4, iblock) >= width ) exit
              end do
              ! add block
              if( iblock > nblock ) then
                nblock = nblock + 1
                tmp(:, nblock) = [irow, icol, height, width]
              end if
              ! increment row
              irow = irow + height
            end do
          end do

          if( allocated( blocks ) ) deallocate( blocks )
          if( nblock > 0 ) then
            allocate( order(nblock), sort(2, nblock) )
            sort(1, :) = tmp(3, 1:nblock) * tmp(4, 1:nblock)
            sort(2, :) = tmp(3, 1:nblock)
            order = sort_index_2d( 2, nblock, -sort, 2 )
            allocate( blocks, source=tmp(:, order) )
            deallocate( order, sort )
          else
            allocate( blocks(4, 0) )
          end if

          deallocate( tmp )
        end function find_blocks
    end function gen_schedule


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
    ! Old k-point partitioning routines.      !
    !+++++++++++++++++++++++++++++++++++++++++++!
    ! These have be superceded but still require replacing in:
    !
    ! For procofk:
    ! * getevalfv.f90, putevalfv.f90
    ! * getevecfv.f90, getevecsv.f90, putevecfv.f90, putevalsv.f90
    ! * getoccsv.f90, putoccsv.f90
    !
    ! For firstk and lastk:
    ! * src_hybrids/putvxnl.f90
    ! * Several src_xs/src_rttddft/ routines (sigh)

   !> Number of k-points on MPI process.
   !> Only used by the depreciated routines, firstk and lastk.
   function nofk(process, nkpt)
     integer(4) :: nofk
     integer(4), intent(in) :: process
     integer, intent(in) :: nkpt
     nofk = nkpt / procs
     if((mod(nkpt, procs) .gt. process)) nofk = nofk + 1
   end function nofk

   function firstk(process, nkpt)
      integer(4) :: firstk
      integer(4), intent(in) :: process
      integer, intent(in) :: nkpt
      integer(4)::i
      firstk = 1
      do i = 0, process - 1
         firstk = firstk + nofk(i, nkpt)
      end do
   end function firstk

   function lastk(process, nkpt)
      integer(4) :: lastk, i
      integer(4), intent(in) :: process
      integer, intent(in) :: nkpt
      lastk = 0
      do i = 0, process
         lastk = lastk + nofk(i, nkpt)
      end do
   end function lastk

   function procofk(k, nkpt)
      integer(4) :: procofk
      integer(4), intent(in) :: k
      integer, intent(in) :: nkpt
      integer(4) :: iproc
      procofk = 0
      do iproc = 0, procs - 1
         if(k .gt. lastk(iproc, nkpt)) procofk = procofk + 1
      end do
   end function procofk

!> Find a 2D grid distribution, given n_processes and n_groups.
!>
!> The number of processes must be equally divisible by the number of groups,
!> else a 2D grid cannot be constructed.
!>
subroutine find_2d_grid(n_processes, n_groups, rows_per_group, cols_per_group, group_label)
  use precision, only: dp

  !> Number of processes
  integer, intent(in) :: n_processes
  !> Number of groups
  integer, intent(in) :: n_groups
  !> Number of rows per group
  integer, intent(out) :: rows_per_group
  !> Number of columns per group
  integer, intent(out) :: cols_per_group
  !> Group label for error message
  character(len=*), intent(in), optional :: group_label

  !> Small, positive tolerance
  real(dp), parameter :: tol = 1.e-8_dp
  !> Processes per group
  integer :: processes_per_group
  !> Error message
  character(len=100) :: error_message

  !TODO(Alex/Peter) Issue 79. SIRUS Integration.
  ! Is this the correct behaviour, or should the routine throw an error?
  if (n_groups < 1) then
    ! One process per group
    cols_per_group = 1
    rows_per_group = 1
    return
  endif

  if (mod(n_processes, n_groups) /= 0) then
    error_message =  "Error(find_2d_grid): Number of processes not divisible by number of groups."
    if (present(group_label)) then
      error_message = error_message//" "//trim(adjustl(group_label))
    endif
    call terminate_mpi_env(mpiglobal, error_message//" groups")
  endif

  processes_per_group = n_processes / n_groups
  cols_per_group = sqrt(real(processes_per_group, dp) + tol)

  do while (mod(processes_per_group, cols_per_group) /= 0)
    cols_per_group = cols_per_group - 1
  end do
  rows_per_group = processes_per_group / cols_per_group

end subroutine find_2d_grid


!> Distribute n_elements between the total number of processes assigned to the communicator of mpi_env.
!>
!> The routine returns the index of the first and last element of a continguous sub-vector, found by
!> dividing up a contiguous vector of n_elements amongst N mpi processes. Indexing begins at 1 for
!> rank = 0, and ends at n_elements for rank = (n_processes - 1).
!>
!> If n_elements is divisible by the number of processes (mpi_env%procs), each process gets the same
!> number of elements. If n_elements is not evenly divisible by the number of processes, remaining elements
!> are distributed amongst the lowest ranks. The total number of elements for a given rank is always
!> defined as (last - first + 1).
!>
!> If n_elements is less than mpi_env%procsm or n_elements == 0, the surplus processes will be assigned
!> (first = 0 ; last = -1) such that the body of a do loop using these limits will not be evaluated, and
!> the total number of elements (last - first + 1) = 0.
!>
!> The routine will terminate in firstofset if:
!>   rank < 0            ! MPI rank is negative. Should be impossible with this API.
!>   n_procs < 1         ! 0 or negative MPI ranks. Should be impossible with this API.
!>   n_procs < rank + 1  ! Rank exceeds total number of MPI processes assigned to this comm.
!>                       ! Should be impossible  with this API.
!
  subroutine distribute_loop(mpi_env, n_elements, first, last)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Total number of elements to distribute.
    integer, intent(in) :: n_elements
    !> First index of the current subset.
    integer, intent(out) :: first
    !> Last index of the current subset.
    integer, intent(out) :: last

    if (n_elements == 0) then
      first = 0
      last = -1
      return
    endif

    first = firstofset(mpi_env%rank, n_elements, mpi_env%procs)
    last = lastofset(mpi_env%rank, n_elements, mpi_env%procs)
  end subroutine distribute_loop




end module modmpi
