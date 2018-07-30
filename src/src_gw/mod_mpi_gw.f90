!----------------------------------------------------
! MPI variables and interface functions
! Can be compiled without mpi libraries installed
!
! Author: Hong Jiang (FHI-gap code)
!----------------------------------------------------
module mod_mpi_gw

    use modmpi
    implicit none
    
    integer, parameter :: nmax_procs=1000

    integer :: myrank     ! the rank in MPI_COMM_WORLD
    integer :: myrank_row ! the rank of the current process in the row commicator 
    integer :: myrank_col ! the rank of the current process in the column commicator
    
    integer :: mycomm
    integer :: mycomm_row ! the communicator corresponding to rows  
    integer :: mycomm_col ! the communicator corresponding to columns 

    integer :: nproc_tot  ! total number of processes
    integer :: nproc_row  ! the number of process of each row (the number of colunms) 
    integer :: nproc_col  ! the number of process of each colunm (the number of rows) 
 
    integer :: iqstart, iqend
    integer :: iqcnt(0:nmax_procs), iqdsp(0:nmax_procs)
    integer :: iomstart, iomend
    integer :: iomcnt(0:nmax_procs), iomdsp(0:nmax_procs)

#ifdef MPI
    interface mpi_sum_array
      module procedure mpi_sum4c
      module procedure mpi_sum3c
      module procedure mpi_sum2c
      module procedure mpi_sum1c
    end interface 
#endif

contains

  !-----------------------------------------------------------------------------
    subroutine init_mpi_gw()
#ifdef MPI
        integer :: fid
        integer :: namelen
        character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
        
        mycomm = MPI_COMM_WORLD

        ! Global MPI environment is already initialized in exciting (modmpi.f90)
        call MPI_Comm_size(mycomm,nproc_tot,ierr)
        call MPI_Comm_rank(mycomm,myrank,ierr)
        
        call MPI_get_processor_name(processor_name,namelen,ierr)
        !write(*,*)
        !write(*,*) "Process ", myrank, " of ", nproc_tot, &
        !&          " running on ", trim(processor_name)
        !write(*,*)
#else
        nproc_tot = 1
        nproc_col = 1
        nproc_row = 1
        myrank = 0
        myrank_col = 0
        myrank_row = 0
#endif
    end subroutine

  !-----------------------------------------------------------------------------    
    subroutine mpi_set_range(np,ip,num,i0,istart,iend,icnts,idisp)
    ! This is a general subroutine that separate i0..i0+num-1 equally to np processes 
    ! in case mod(num,nranks)==0, the residual is put to last mod(num,np) processes 
    ! This subroutine can be used to set k-point or frequency point range  
        integer, intent(in)  :: np, ip, num, i0
        integer, intent(out) :: istart, iend 
        integer, intent(out), optional :: icnts(0:), idisp(0:)
        integer :: icount(nmax_procs)
        integer :: i
    
        ! Since the calculation for Gamma point, iq=1, at most takes twice 
        ! cpu time than other q-points, we take the folllowing strategy:
        ! if nkpt is equal to an integral times nproc_tot, then divide nkpt equally
        ! otherwise, the residual k-points are equally assigned to non-root processes
        ! starting from the last one
        if (np==1) then 
          istart = i0
          iend = i0+num-1 
          icount(1) = num
        else
          icount(:) = num/np
          do i = 1, mod(num,np)
            icount(np-i+1) = icount(np-i+1)+1
          end do
          istart = i0
          do i = 1, ip
            istart = istart+icount(i)
          end do
          iend = istart+icount(ip+1)-1
        end if 
        if (present(icnts)) then 
          icnts(0:np-1) = icount(1:np)
          if (present(idisp)) then 
            idisp(0) = 0
            do i = 1, np-1
              idisp(i) = idisp(i-1)+icnts(i-1)
            end do
          end if 
        end if 
        return
    end subroutine 

#ifdef MPI
  !-----------------------------------------------------------------------------
    subroutine set_mpi_group(nkp)
    ! this subroutine is to set up the MPI communicators 
    ! all processes are divided into groups
        implicit none
        integer, intent(in) :: nkp 
        integer :: myrow, mycol
        character(len=10) :: str
        
        if (myrank==0) then 
        
          call getenv("MPI_NPROC_COL",str)
          ! write(*,*) "Get MPI_NPROC_COL"
          ! write(*,*) "=>", str
          if (str=='') then 
            nproc_col = 1
          else 
            call str2int(str,nproc_col,ierr) 
            if (ierr.ne.0) then 
              write(*,*) "WARNING: invalid value for MPI_NPROC_IN_GRP"
              nproc_col = 1
            end if 
          end if 

          ! if nproc_col is not defined then try to set it in terms of nproc_tot
          ! and nkp automatically
          if ((nproc_col==1) .and. (nproc_tot>nkp)) then 
            if (mod(nproc_tot,nkp)==0) then
              nproc_col = nproc_tot/nkp
            else
              nproc_col = nproc_tot/nkp+1
            end if
          end if 
          ! write(*,*) "Broadcast MPI_NPROC_COL"
        end if 

        call MPI_Bcast(nproc_col,1,MPI_INTEGER,0,mycomm,ierr)
        call MPI_Barrier(mycomm,ierr)
    
        if (mod(nproc_tot,nproc_col)==0) then 
          nproc_row = nproc_tot/nproc_col
        else 
          nproc_row = nproc_tot/nproc_col+1
        end if 

        ! if (myrank==0) then 
        !   write(*,*) " nproc_tot =", nproc_tot
        !   write(*,*) " nproc_row =", nproc_row
        !   write(*,*) " nproc_col =", nproc_col
        ! end if

        myrow = mod(myrank,nproc_col)
        mycol = myrank/nproc_col
        
        call MPI_Comm_split(mycomm,myrow,myrank,mycomm_row,ierr)
        call MPI_Comm_size(mycomm_row,nproc_row,ierr)
        call MPI_Comm_rank(mycomm_row,myrank_row,ierr)
        
        call MPI_Comm_split(mycomm,mycol,myrank,mycomm_col,ierr)
        call MPI_Comm_size(mycomm_col,nproc_col,ierr)
        call MPI_Comm_rank(mycomm_col,myrank_col,ierr)
        
        return
    end subroutine 

  !-----------------------------------------------------------------------------
    subroutine mpi_free_group
      call MPI_Comm_free(mycomm_row,ierr)
      call MPI_Comm_free(mycomm_col,ierr)
    end subroutine 

  !-----------------------------------------------------------------------------
  ! Routines to reduce the complex arrays of different dimentionalities
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
    subroutine mpi_sum4c(iflag,a,n1,n2,n3,n4,comm)
    ! This subroutine  is used to sum complex arrays of the form, 
    ! a(1:n1,1:n2,1:n3,1:n4) that is calculated at different processes
    ! iflag=0 --- reduce to the root process
        implicit none
        integer, intent(in) :: iflag
        integer, intent(in) :: n1, n2, n3, n4
        complex(8), intent(inout) :: a(n1,n2,n3,n4)    
        integer, intent(inout) :: comm
        integer(4) :: rank, ierr

        if (comm<0) comm = MPI_COMM_WORLD
        call MPI_Comm_rank(comm,rank,ierr)
        call MPI_Barrier(comm,ierr)
        if (iflag==0) then
          if (rank==0) then 
            call MPI_Reduce(MPI_IN_PLACE,a,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          else 
            call MPI_Reduce(a,0,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          end if 
        else
          call MPI_AllReduce(MPI_IN_PLACE,a,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX, &
          &  MPI_SUM,comm,ierr)
        end if
        return
    end subroutine

  !-----------------------------------------------------------------------------
    subroutine mpi_sum3c(iflag,a,n1,n2,n3,comm)
    ! This subroutine  is used to sum complex arrays of the form, 
    ! a(1:n1,1:n2,1:n3) that is calculated at different processes
    ! iflag=0 --- reduce to the root process 
        implicit none
        integer, intent(in) :: iflag
        integer, intent(in) :: n1, n2, n3
        complex(8), intent(inout) :: a(n1,n2,n3)    
        integer, intent(inout) :: comm
        integer(4) :: rank, ierr

        if (comm<0) comm = MPI_COMM_WORLD
        call MPI_COMM_RANK(comm,rank,ierr)
        call MPI_BARRIER(comm,ierr)
        if (iflag==0) then
          if (rank==0) then 
            call MPI_REDUCE(MPI_IN_PLACE,a,n1*n2*n3,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          else 
            call MPI_REDUCE(a,0,n1*n2*n3,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          end if
        else
          call MPI_ALLREDUCE(MPI_IN_PLACE,a,n1*n2*n3,MPI_DOUBLE_COMPLEX, &
          &  MPI_SUM,comm,ierr)
        endif
        return
    end subroutine
    
  !-----------------------------------------------------------------------------
    subroutine mpi_sum2c(iflag,a,n1,n2,comm)
    ! This subroutine  is used to sum complex arrays of the form, 
    ! a(1:n1,1:n2) that is calculated at different processes
    ! iflag=0 --- reduce to the root process
        implicit none
        integer, intent(in) :: iflag
        integer, intent(in) :: n1, n2
        complex(8), intent(inout) :: a(n1,n2)    
        integer, intent(inout) :: comm
        integer(4) :: rank, ierr  

        if (comm<0) comm = MPI_COMM_WORLD
        call MPI_COMM_RANK(comm,rank,ierr)
        call MPI_BARRIER(comm,ierr)
        if (iflag==0) then 
          if (rank==0) then 
            call MPI_REDUCE(MPI_IN_PLACE,a,n1*n2,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          else 
            call MPI_REDUCE(a,0,n1*n2,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          end if 
        else
          call MPI_ALLREDUCE(MPI_IN_PLACE,a,n1*n2,MPI_DOUBLE_COMPLEX, &
          &  MPI_SUM,comm,ierr)
        end if 
        return
    end subroutine

  !-----------------------------------------------------------------------------
    subroutine mpi_sum1c(iflag,a,n1,comm)
    ! This subroutine  is used to sum complex arrays of the form, 
    ! a(1:n1) that is calculated at different processes
    ! iflag=0 --- reduce to the root process
        implicit none
        integer, intent(in) :: iflag
        integer, intent(in) :: n1
        complex(8), intent(inout) :: a(n1)    
        integer, intent(inout) :: comm
        integer(4) :: rank, ierr

        if (comm<0) comm = MPI_COMM_WORLD
        call MPI_COMM_RANK(comm,rank,ierr)
        call MPI_BARRIER(comm,ierr)
        if (iflag==0) then
          if (rank==0) then 
            call MPI_REDUCE(MPI_IN_PLACE,a,n1,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          else
            call MPI_REDUCE(a,0,n1,MPI_DOUBLE_COMPLEX, &
            &  MPI_SUM,0,comm,ierr)
          end if 
        else
          call MPI_ALLREDUCE(MPI_IN_PLACE,a,n1,MPI_DOUBLE_COMPLEX, &
          &  MPI_SUM,comm,ierr)
        end if
        return
    end subroutine

#endif

end module
