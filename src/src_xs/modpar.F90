
module modpar
  implicit none

#ifdef MPI
#include "mpif.h"
#endif

  !---------------------------!
  !     general variables     !
  !---------------------------!
  ! current parallel process index
  integer :: myrank,rank
  ! total number of processs
  integer :: nproc
  ! maximum allowed number of processs 
  integer, parameter :: maxproc=100
  ! error
  integer :: ierr

  ! cumulative counter for barriers
  integer :: barcntcum
  data barcntcum /0/
  ! counter for barriers
  integer :: barcnt
  data barcnt /0/
  ! cumulative wall time for barriers
  real(8) :: bartimcum
  data bartimcum /0.d0/
  ! wall time for barriers
  real(8) :: bartim
  data bartim /0.d0/

  ! file extentsion for parallelization
  character(256) :: filextp
  ! characteristic string for parallelization
  character(256) :: spar

  ! idle times for barrier
  integer :: baridl,baridl2

  !--------------------------!
  !     resume variables     !
  !--------------------------!
  ! true if calculation is resumed from checkpoint
  logical :: tresume
  ! task to resume
  integer :: resumetask
  ! maximum number of checkpoints
  integer, parameter :: maxchkpts=2
  ! loop checkpoints to resume
  integer :: resumechkpts(maxchkpts,3)
  ! filename for resuming
  character(256) :: fnresume

  !---------------------!
  !     q-point set     !
  !---------------------!
  ! current initial q-point index
  integer :: qpari
  ! current final q-point index
  integer :: qparf

  !---------------------!
  !     k-point set     !
  !---------------------!
  ! current initial k-point index
  integer :: kpari
  ! current final k-point index
  integer :: kparf

  !---------------------!
  !     w-point set     !
  !---------------------!
  ! current initial w-point index
  integer :: wpari
  ! current final w-point index
  integer :: wparf

contains

  subroutine initpar
    implicit none
    ! default values
    nproc=1
    rank=1
    ! interface to parallel environment
    ! use MPI environment here
#ifdef MPI
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
    rank=myrank+1
#endif
  end subroutine initpar

  subroutine finitpar
    implicit none
    ! interface to parallel environment
    ! use MPI environment here
#ifdef MPI
    call MPI_Finalize(ierr)
    write(*,*) 'MPI finalized with rank,ierr:',rank,ierr
#endif
  end subroutine finitpar

  subroutine resupd(un,tsk,val,str)
    implicit none
    ! arguments
    integer, intent(in) :: un,tsk,val(maxchkpts,3)
    character(*), intent(in) :: str
    ! local variables
    character(*), parameter :: thisnam='resupd'
    character(32) :: fmt
    logical :: existent

    open(un,file=trim(fnresume),action='write',status='replace')
    write(fmt,*) 1+3*maxchkpts
    fmt='('//trim(adjustl(fmt))//'i9,a)'
    write(un,fmt=trim(fmt)) tsk, val, str
    close(un)
  end subroutine resupd
  
  subroutine resread(un,tsk,val,tres)
    implicit none
    ! arguments
    integer, intent(in) :: un
    integer, intent(out) :: tsk,val(maxchkpts,3)
    logical :: tres
    ! local variables
    character(*), parameter :: thisnam='resread'

    inquire(file=trim(fnresume),exist=tres)
    if (.not.tres) then
       ! no checkpoint present
       return
    end if
    open(un,file=trim(fnresume),action='read')
    read(un,*) tsk, val
    close(un)
  end subroutine resread


  subroutine getrange(rank,nproc,npt,i,f)
    implicit none
    ! arguments
    integer, intent(in) :: rank,nproc,npt
    integer, intent(out) :: i,f
    ! local variables
    character(*), parameter :: thisnam = 'getrange'
    integer, allocatable :: t(:), it(:), ft(:)
    integer :: a,b,s

    ! only one process
    if (nproc.eq.1) then
       i=1
       f=npt
       return
    end if
    allocate(t(0:nproc),it(nproc),ft(nproc))
    t(:) = 0
    ! number of points for each process
    do a=1,npt
       b=mod(a,nproc)
       if (b.eq.0) b=nproc
       t(b) = t(b)+1
    end do
    ! first and last point for each process
    do b=1,nproc
       s = sum(t(:b-1))
       it(b) = s + 1
       ft(b) = s + t(b) 
    end do
    i = it(rank)
    f = ft(rank)
    deallocate(it,ft)
  end subroutine getrange


  subroutine getridx(nproc,npt,i,irec)
    implicit none
    ! arguments
    integer, intent(in) :: nproc,npt,i
    integer, intent(out) :: irec
    ! local variables
    character(*), parameter :: thisnam = 'getridx'
    integer :: irnk,pi,pf

    ! only one process
    if (nproc.eq.1) then
       irec=i
       return
    end if
    ! get range for process
    do irnk=1,nproc
      call getrange(irnk,nproc,npt,pi,pf)
      if ((pi.le.i).and.(i.le.pf)) then
        ! record position for current k-point
        irec = i-pi+1
        exit
      end if
    end do
  end subroutine getridx


  subroutine barrier(rank,nproc,un,async,string)
    implicit none
    ! arguments
    integer, intent(in) :: rank,nproc,un,async
    character(*) :: string
    ! local variables
    character(*), parameter :: thisnam = 'barrier'
    integer :: irank, cc
    integer :: cntr,ti,tf
    character(256) :: ranks
    logical :: existent,wait

    ! do nothing if only one process
    if (nproc.eq.1) return
    
    call system_clock(COUNT_RATE=cntr)
    call system_clock(COUNT=ti)

    ! increase counter
    barcntcum=barcntcum+1
    barcnt=barcnt+1
    
    ! call the MPI barrier
#ifdef MPI
    call MPI_barrier(mpi_comm_world,ierr)
#endif
    
#ifndef MPI
    if (barcntcum.gt.1000) then
       write(*,'(a)') "Info("//trim(thisnam)//"): maximum barrier &
            &count exceeded"
       stop
    end if

    write(*,*) 'begin of barrier for counter/process:', barcntcum,rank

    if (rank.eq.1) then
       ! wait until other processes have reached barrier
       wait=.true.
       do while (wait)
          cc=0
          do irank=2,nproc
             write(ranks,'(i3.3,".",i3.3)') irank,barcntcum
             inquire(file=trim(string)//trim(ranks),exist=existent)
             if (existent) cc = cc + 1
          end do
          if (cc.eq.(nproc-1)) wait=.false.
          call sleepifc(1)
       end do
       ! mark this process as finished
       write(ranks,'(i3.3,".",i3.3)') rank,barcntcum
       open(un,file=trim(string)//trim(ranks),action='write')
       close(un)
       ! wait unitl barrier disappears for other processes
       wait=.true.
       do while (wait)
          cc=0
          do irank=2,nproc
             write(ranks,'(i3.3,".",i3.3)') irank,barcntcum
             inquire(file=trim(string)//trim(ranks),exist=existent)
             if (.not.existent) cc = cc + 1
          end do
          if (cc.eq.(nproc-1)) wait=.false.
          call sleepifc(1)
       end do
       ! remove tag
       write(ranks,'(i3.3,".",i3.3)') rank,barcntcum
       open(un,file=trim(string)//trim(ranks),action='write')
       close(un,status='delete') 
    else
       ! mark reach of barrier
       write(ranks,'(i3.3,".",i3.3)') rank,barcntcum
       open(un,file=trim(string)//trim(ranks),action='write')
       close(un)
       ! wait for process one to finish barrier
       wait=.true.
       write(ranks,'(i3.3,".",i3.3)') 1,barcntcum
       do while (wait)
          inquire(file=trim(string)//trim(ranks),exist=existent)
          if (existent) wait=.false.
          ! sleep for one second
          call sleepifc(1)
       end do
       ! remove tag
       write(ranks,'(i3.3,".",i3.3)') rank,barcntcum
       open(un,file=trim(string)//trim(ranks),action='write')
       close(un,status='delete') 
       call sleepifc(baridl2)
    end if
#endif

    ! for I/O release
    call sleepifc(async*rank)
    ! take overhead into account
    call sleepifc(baridl)
    call system_clock(COUNT=tf)
    bartimcum=bartimcum+dble(tf-ti)/dble(cntr)
    bartim=bartim+dble(tf-ti)/dble(cntr)

#ifndef MPI
    write(*,'(a,2i9)') 'end of barrier for counter/process  :', barcntcum,rank
#endif

  end subroutine barrier

end module modpar
