module bse_transitions 


  implicit none 

  private 
  public :: transition_type

  !> Index of the first unoccupied band index in the band limits lookup table per k-point.
  integer, parameter, public :: uf=1
  !> Index of the last unoccupied band index in the band limits lookup table per k-point.
  integer, parameter, public :: ul=2
  !> Index of the first occupied band index in the band limits lookup table per k-point.
  integer, parameter, public :: of=3
  !> Index of the last noccupied band index in the band limits lookup table per k-point.
  integer, parameter, public :: ol=4
  !> Index of the first transition index in the band limits lookup table per k-point.
  integer, parameter, public :: tf=5
  !> Index of the last transition index in the band limits lookup table per k-point.
  integer, parameter, public :: tl=6

  !> Integer for uninitialized type
  integer, parameter, private :: uninitialized = -1

  type transition_type
    integer :: n_tot, n_allowed, n_k, n_o, n_u
    integer, allocatable :: n_o_per_k(:), n_u_per_k(:), n_per_k(:)
    integer, allocatable :: o_first(:), o_last(:), u_first(:), u_last(:), first(:), last(:)
    logical, allocatable :: allowed(:)

    contains 

    procedure :: initialize
    procedure :: finalize 
  end type transition_type

  contains 

  subroutine initialize(this, band_idx, allowed)
    class(transition_type) :: this 
    integer, intent(in) :: band_idx(:, :)
    integer, intent(in) :: allowed(:)

    this%allowed = allowed

    this%o_first = band_idx(of, :)
    this%o_last = band_idx(ol, :)

    this%u_first = band_idx(uf, :)
    this%u_last = band_idx(ul, :)

    this%first = band_idx(tf, :)
    this%last = band_idx(tl, :)

    this%n_o_per_k = band_idx(ol, :) - band_idx(of, :) + 1
    this%n_u_per_k = band_idx(ul, :) - band_idx(uf, :) + 1
    this%n_per_k = this%n_o_per_k * this%n_u_per_k

    this%n_k = size(band_idx, 2)
    this%n_tot = sum(this%n_per_k)
    this%n_allowed = count(this%allowed)
    this%n_o = sum(this%n_o_per_k)
    this%n_u = sum(this%n_u_per_k)
  end subroutine 

  subroutine finalize(this)
    class(transition_type) :: this

    if(allocated(this%o_first)) deallocate(this%o_first)
    if(allocated(this%o_last)) deallocate(this%o_last)

    if(allocated(this%u_first)) deallocate(this%u_first)
    if(allocated(this%u_last)) deallocate(this%u_last)

    if(allocated(this%first)) deallocate(this%first)
    if(allocated(this%last)) deallocate(this%last)

    if(allocated(this%n_o_per_k)) deallocate(this%n_o_per_k)
    if(allocated(this%n_u_per_k)) deallocate(this%n_u_per_k)
    if(allocated(this%n_per_k)) deallocate(this%n_per_k)

    this%n_k = uninitialized
    this%n_tot = uninitialized
    this%n_allowed = uninitialized
    this%n_o = uninitialized
    this%n_u = uninitialized
  end subroutine finalize

end module bse_transitions