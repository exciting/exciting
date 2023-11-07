!> Module provides type for Brillouin zone paths.
module bz_path
  use precision, only: dp
  use asserts, only: assert

  implicit none
  private

  ! numbers smaller that this are treated as zero
  real(dp), parameter :: eps_zero = 1e-64_dp

  !> Vertex type.
  type vertex
    !> vertex in lattice coordinates
    real(dp) :: coord_lat(3)
    !> vertex in Cartensian coordinates
    real(dp) :: coord_cart(3)
    !> 1d distance of vertex along path
    real(dp) :: distance
    !> vertex label
    character(:), allocatable :: label
    !> break path after this vertex
    logical :: break = .false.
  end type vertex

  !> Point type.
  type point
    !> point in lattice coordinates
    real(dp) :: coord_lat(3)
    !> point in Cartesian coordinates
    real(dp) :: coord_cart(3)
    !> 1d distance of point along path
    real(dp) :: distance
  end type point

  !> Brillouin zone path type.
  type bz_path_type
    !> number of points along path
    integer, public :: num_points = 0
    !> list of points
    type(point), allocatable, public :: points(:)
    !> number of vertices along path
    integer, public :: num_vertices = 0
    !> list of vertices
    type(vertex), allocatable, public :: vertices(:)
  end type bz_path_type
  !> constructor
  interface bz_path_type
    module procedure :: new_bz_path_from_vertices, new_bz_path_from_path_type
  end interface

  public :: bz_path_type

contains

  !> Construct path from a given set of vertices.
  function new_bz_path_from_vertices( vertices, num_points, gamma_offset ) result( this )
    use grid_utils, only: linspace
    !> list of vertices that define the path
    type(vertex), intent(in) :: vertices(:)
    !> number of points along path
    integer, intent(in) :: num_points
    !> minimum distance from Gamma point do avoid Gamma sampling (default: `0.0_dp`)
    real(dp), optional, intent(in) :: gamma_offset
    !> BZ path
    type(bz_path_type) :: this

    integer :: iv, ip, nv, np, seg_num_points, vi(3)
    real(dp) :: offset, path_length, seg_length, length, vsl(3), vsc(3), vl(3)
    logical :: break_before
    real(dp), allocatable :: segment(:)
    type(vertex), allocatable :: verts(:)

    offset = 0.0_dp
    if( present( gamma_offset ) ) offset = gamma_offset

    call assert( size( vertices ) > 0, &
      'No vertices given.' )
    call assert( num_points >= size( vertices ), &
      'Number of points must not be smaller than number of vertices.' )
    call assert( (size( vertices ) == 1 .and. num_points == 1) .or. (size( vertices ) > 1), &
      'Number of points must equal 1 of only 1 vertex is given.' )

    this%num_points = num_points
    this%num_vertices = size( vertices )
    allocate( this%points(this%num_points) )
    allocate( this%vertices, source=vertices )

    ! set vertex distances and get total path length
    path_length = 0.0_dp
    this%vertices(1)%distance = 0.0_dp
    do iv = 2, this%num_vertices
      if( .not. this%vertices(iv-1)%break ) &
        path_length = path_length + norm2( this%vertices(iv)%coord_cart - this%vertices(iv-1)%coord_cart )
      this%vertices(iv)%distance = path_length
    end do

    ! add additional vertices for Gamma offset
    if( offset == 0.0_dp ) then
      nv = this%num_vertices
      allocate( verts, source=this%vertices )
    else
      break_before = .true.
      allocate( verts(2*this%num_vertices) )
      nv = 0
      do iv = 1, this%num_vertices
        vl = this%vertices(iv)%coord_lat
        call r3frac( eps_zero, vl, vi )
        if( norm2( vl ) < eps_zero .and. .not. break_before .and. .not. this%vertices(iv)%break .and. iv < this%num_vertices ) then
          nv = nv + 1
          verts(nv) = this%vertices(iv)
          verts(nv)%break = .true.
          nv = nv + 1
          verts(nv) = this%vertices(iv)
          break_before = verts(nv)%break
        else
          nv = nv + 1
          verts(nv) = this%vertices(iv)
          break_before = verts(nv)%break
        end if
      end do
    end if
    verts(nv)%break = .true.
    call assert( this%num_points >= nv, &
      'Number of points must not be smaller than number of vertices. &
      (Additional vertices have been added due to Gamma point offset.)' )
    call assert( .not. (all( [(verts(iv)%break, iv=1, nv)] ) .and. this%num_points /= nv), &
      'All vertices are isolated. Number of points must equal number of vertices.' )

    ! set points
    iv = 1; np = 0; length = 0.0_dp; break_before = .true.
    do while( iv <= nv )
      ! isolated vertex
      if( verts(iv)%break .and. break_before ) then
        np = np + 1
        this%points(np)%coord_lat = verts(iv)%coord_lat
        this%points(np)%coord_cart = verts(iv)%coord_cart
        this%points(np)%distance = verts(iv)%distance
      ! segment
      else
        vsl = verts(iv+1)%coord_lat - verts(iv)%coord_lat
        vsc = verts(iv+1)%coord_cart - verts(iv)%coord_cart
        seg_length = norm2( vsc )
        if( verts(iv+1)%break ) then
          vl = verts(iv+1)%coord_lat
          call r3frac( eps_zero, vl, vi )
          seg_num_points = nint( (this%num_points - np - 1) * (seg_length / (path_length - length)) )
          segment = [(dble(ip-1)/dble(seg_num_points), ip=1, seg_num_points+1)]
          if( norm2( vl ) < eps_zero ) segment(seg_num_points+1) = 1.0_dp - offset / norm2( vsl )
        else
          seg_num_points = nint( (this%num_points - np) * (seg_length / (path_length - length)) )
          segment = [(dble(ip-1)/dble(seg_num_points), ip=1, seg_num_points)]
        end if
        vl = verts(iv)%coord_lat
        call r3frac( eps_zero, vl, vi )
        if( norm2( vl ) < eps_zero ) segment(1) = offset / norm2( vsl )
        do ip = 1, size( segment )
          this%points(np+ip)%coord_lat = verts(iv)%coord_lat + segment(ip) * vsl
          this%points(np+ip)%coord_cart = verts(iv)%coord_cart + segment(ip) * vsc
          this%points(np+ip)%distance = verts(iv)%distance + segment(ip) * seg_length
        end do
        this%points(np+1)%distance = verts(iv)%distance
        np = np + size( segment )
        if( verts(iv+1)%break ) then
          iv = iv + 1
          this%points(np)%distance = verts(iv)%distance
        end if
        break_before = verts(iv)%break
        length = length + seg_length
      end if
      iv = iv + 1
    end do
    
    deallocate( verts )
  end function new_bz_path_from_vertices

  !> Construct path from exciting [[path_type(type)]] object.
  function new_bz_path_from_path_type( bvec, path, gamma_offset ) result( this )
    use modinput, only: path_type
    real(dp), intent(in) :: bvec(3, 3)
    !> `path_type` object
    type(path_type), intent(in) :: path
    !> minimum distance from Gamma point do avoid Gamma sampling (default: `0.0_dp`)
    real(dp), optional, intent(in) :: gamma_offset
    !> BZ path
    type(bz_path_type) :: this

    integer :: iv
    type(vertex), allocatable :: verts(:)

    allocate( verts(size( path%pointarray )) )

    do iv = 1, size( verts )
      verts(iv)%coord_lat = path%pointarray(iv)%point%coord
      verts(iv)%coord_cart = matmul( bvec, verts(iv)%coord_lat )
      verts(iv)%label = trim( adjustl( path%pointarray(iv)%point%label ) )
      verts(iv)%break = path%pointarray(iv)%point%breakafter
    end do

    if( present( gamma_offset ) ) then
      this = new_bz_path_from_vertices( verts, path%steps, gamma_offset=gamma_offset )
    else
      this = new_bz_path_from_vertices( verts, path%steps )
    end if
  end function new_bz_path_from_path_type

end module bz_path
