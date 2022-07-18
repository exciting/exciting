!> Module generating and setting seeds for random number generation
module seed_generation
  use asserts, only: assert
  use exciting_mpi, only: mpiinfo, xmpi_bcast
  use modmpi, only: terminate

  implicit none

  private
  public :: set_seed

  !> Set the seed for random number generation.
  interface set_seed
    procedure :: select_set_seed, set_seed_int
  end interface

  contains

  !> Wrapper for seed generation.
  !> For details see
  !>
  !> - [[set_seed_to_system_time]]
  !>
  !> - [[set_seed_to_system_time_MPI]]
  !>
  !> This routine is not thread save!
  subroutine select_set_seed(source)
    !> Source for seed generation. This can either be
    !>
    !> - "fixed" for setting a fixed number (see [[default_seed]])
    !>
    !> - "clock" for "system_time" for setting the system time as seed (not thread safe)
    character(*), intent(in) :: source

    integer, parameter :: default_seed = 0 ! set fixed seed
    integer :: system_time                 ! set system time as seed

    select case(source)
      ! Use a fixed seed
      case('fixed')
        call set_seed(default_seed)

      ! Use system time as seed
      case("clock")
        call system_clock(count = system_time)
        call set_seed(system_time)

      case default
        call terminate('source is none of "fixed" or "clock".')
    end select
  end subroutine select_set_seed

  !> Set an integer seed for random number generation.
  subroutine set_seed_int(seed_in)
    !> Seed to set.
    integer, intent(in) :: seed_in

    integer :: len_seed
    integer, allocatable :: seed(:)

    call random_seed(size = len_seed)
    allocate(seed(len_seed))
    seed = seed_in
    call random_seed(put=seed)
  end subroutine set_seed_int

end module seed_generation