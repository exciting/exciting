!> Module to add status reports to any loop.
module m_status_report
  use precision, only: dp, i32
  use modmpi, only: mpiglobal
  use mod_omp_utils, only: omp_thread_num

  implicit none

  private
  public :: status_report_t

  !> Type that contains all information needed for performing the status reports.
  type status_report_t
    !> index to keep track of when to print a progress message
    integer(i32) :: next_report
    !> full loop index
    integer(i32) :: index
    !> number of status reports to print
    integer(i32) :: nreports
    !> total number of elements in the loop
    integer(i32) :: niter
    !> where to print the status messages
    integer(i32) :: out_unit
    !> name of the caller, used as identifier in the status message
    character(:), allocatable :: calling_loop_name
    !> when the loop started
    real(dp) :: start_time
    !> stores at which full loop indices a message should be printed
    integer(i32), allocatable :: report_indices(:)
  contains
    procedure, public :: init
    procedure, public :: update
    procedure, public :: delete
    procedure, private :: print_progress
  end type

  contains
    !> Initialize the status report type
    subroutine init(self, nreports, niter, calling_loop_name, out_unit, start_time)
      class(status_report_t), intent(inout) :: self
      !> number of reports requested, might be lowered to the total number of elements in the loop
      integer(i32), intent(in) :: nreports
      !> total number of elements in the loop
      integer(i32), intent(in) :: niter
    !> where to print the status messages
      integer(i32), intent(in) :: out_unit
    !> name of the caller, used as identifier in the status message
      character(*), intent(in) :: calling_loop_name
    !> when the loop started
      real(dp), intent(in) :: start_time

      integer(i32) :: i

      ! use less steps if not enough are available
      self%niter = niter
      self%nreports = min(nreports, self%niter)
      self%out_unit = out_unit
      self%calling_loop_name = calling_loop_name
      self%start_time = start_time

      self%index = 0

      self%next_report = 1
      allocate(self%report_indices(self%nreports))
      do i=1,self%nreports
        self%report_indices(i) = ceiling( dble(i) * dble(self%niter) / dble(self%nreports) )
      end do

      ! give some hint in the out_unit to structure the output
      if (mpiglobal%rank == 0 .and. self%nreports > 0 .and. omp_thread_num() == 0) then
        write(self%out_unit,'("Starting ", A, " loop ...")') self%calling_loop_name
        flush(self%out_unit)
      end if

    end subroutine

    !> Update, should be called during each iteration. Decides on whether to print a message.
    subroutine update(self)
      class(status_report_t), intent(inout) :: self

      self%index = self%index + 1
      if (self%next_report <= self%nreports) then
        if (self%index == self%report_indices(self%next_report)) then
          call self%print_progress()
          self%next_report = self%next_report + 1
        end if
      end if

    end subroutine

    !> Call after the loop to free memory.
    subroutine delete(self)
      class(status_report_t), intent(inout) :: self

      if (allocated(self%calling_loop_name)) deallocate(self%calling_loop_name)
      if (allocated(self%report_indices)) deallocate(self%report_indices)

      ! just an empty line to give some structure in the out_unit
      if (mpiglobal%rank == 0 .and. self%nreports > 0 .and. omp_thread_num() == 0) then
        write(self%out_unit,'("")')
        flush(self%out_unit)
      end if

    end subroutine

    !> Print a progress message with current percentage and timing.
    subroutine print_progress(self)
      class(status_report_t), intent(inout) :: self

      real(dp) :: percentage
      real(dp) :: current_time
      integer :: h, m, s
      character(len=8) :: elapsed_str

      percentage = 100.d0 * dble(self%index) / dble(self%niter)
      call timesec(current_time)

      s = int(current_time - self%start_time)
      h = s / 3600
      m = mod(s,3600) / 60
      s = mod(s,60)
      write(elapsed_str,'(I2.2,":",I2.2,":",I2.2)') h, m, s

      if (mpiglobal%rank == 0 .and. omp_thread_num() == 0) then
        write(self%out_unit,'("Progress info(",A," loop): ",I0,"/",I0," (",F5.1,"%), elapsed: ",A)') &
                trim(self%calling_loop_name), self%index, self%niter, percentage, trim(elapsed_str)
        flush(self%out_unit)
      end if

    end subroutine print_progress

end module m_status_report
