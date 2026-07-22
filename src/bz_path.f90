!> Module provides type for Brillouin zone paths.
module bz_path
  use precision, only: dp
#include "asserts.fpp"
  implicit none
  private

  ! numbers smaller that this are treated as zero
  real(dp), parameter :: eps_zero = epsilon( 1.0_dp )

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
    
    contains
      !> return JSON string
      procedure, public :: to_json => vertex_to_json
  end type vertex

  !> Point type.
  type point
    !> point in lattice coordinates
    real(dp) :: coord_lat(3)
    !> point in Cartesian coordinates
    real(dp) :: coord_cart(3)
    !> 1d distance of point along path
    real(dp) :: distance
    
    contains
      !> return JSON string
      procedure, public :: to_json => point_to_json
  end type point

  !> Segment type.
  type segment
    !> bounding vertices
    type(vertex), public :: vertices(2)
    !> indices of endpoints
    integer, public :: end_point_idx(2)
    !> number of points in segment (including endpoints)
    integer, public :: num_points = 0
    !> segment length
    real(dp), public :: length
    
    contains
      !> return JSON string
      procedure, public :: to_json => segment_to_json
  end type segment

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
    !> number of segments path is made of
    integer, public :: num_segments
    !> list of segments
    type(segment), allocatable, public :: segments(:)
    !> total length of path
    real(dp), public :: length
    
    contains
      !> return JSON string
      procedure, public :: to_json => path_to_json
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

    integer :: iv, ip, is, nv, np, num_seg_points, vi(3)
    real(dp) :: offset, length, vsl(3), vsc(3), vl(3)
    logical :: break_before
    real(dp), allocatable :: segment(:)
    type(vertex), allocatable :: verts(:)

    offset = 0.0_dp
    if( present( gamma_offset ) ) offset = gamma_offset

    CALL_ASSERT( size( vertices ) > 0,  'No vertices given.' )
    CALL_ASSERT( num_points >= size( vertices ),  'Number of points must not be smaller than number of vertices.' )
    CALL_ASSERT( (size( vertices ) == 1 .and. num_points == 1) .or. (size( vertices ) > 1),  'Number of points must equal 1 of only 1 vertex is given.' )

    this%num_points = num_points
    this%num_vertices = size( vertices )
    this%num_segments = this%num_vertices - 1 - count( [(vertices(iv)%break, iv=1, this%num_vertices-1)] )
    allocate( this%points(this%num_points) )
    allocate( this%vertices, source=vertices )
    allocate( this%segments(this%num_segments) )

    ! set vertex distances and get total path length
    this%length = 0.0_dp
    this%vertices(1)%distance = 0.0_dp
    do iv = 2, this%num_vertices
      if( .not. this%vertices(iv-1)%break ) &
        this%length = this%length + norm2( this%vertices(iv)%coord_cart - this%vertices(iv-1)%coord_cart )
      this%vertices(iv)%distance = this%length
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
    CALL_ASSERT( this%num_points >= nv,  'Number of points must not be smaller than number of vertices.  (Additional vertices have been added due to Gamma point offset.)' )
    CALL_ASSERT( .not. (all( [(verts(iv)%break, iv=1, nv)] ) .and. this%num_points /= nv),  'All vertices are isolated. Number of points must equal number of vertices.' )

    ! set points and segments
    iv = 1; is = 0; np = 0; length = 0.0_dp; break_before = .true.
    do while( iv <= nv )
      ! isolated vertex
      if( verts(iv)%break .and. break_before ) then
        np = np + 1
        this%points(np)%coord_lat = verts(iv)%coord_lat
        this%points(np)%coord_cart = verts(iv)%coord_cart
        this%points(np)%distance = verts(iv)%distance
        is = is + 1
        this%segments(is)%vertices = verts(iv)
        this%segments(is)%end_point_idx = np
        this%segments(is)%num_points = 1
        this%segments(is)%length = 0.0_dp
      ! segment
      else
        vsl = verts(iv+1)%coord_lat - verts(iv)%coord_lat
        vsc = verts(iv+1)%coord_cart - verts(iv)%coord_cart
        is = is + 1
        this%segments(is)%vertices = verts(iv:iv+1)
        this%segments(is)%end_point_idx(1) = np + 1
        this%segments(is)%length = norm2( vsc )
        if( verts(iv+1)%break ) then
          vl = verts(iv+1)%coord_lat
          call r3frac( eps_zero, vl, vi )
          num_seg_points = nint( (this%num_points - np - 1) * (this%segments(is)%length / (this%length - length)) + this%num_points * eps_zero )
          segment = [(dble(ip-1)/dble(num_seg_points), ip=1, num_seg_points+1)]
          if( norm2( vl ) < eps_zero ) segment(num_seg_points+1) = 1.0_dp - offset / norm2( vsl )
          this%segments(is)%num_points = size( segment )
        else
          num_seg_points = nint( (this%num_points - np) * (this%segments(is)%length / (this%length - length)) + this%num_points * eps_zero )
          segment = [(dble(ip-1)/dble(num_seg_points), ip=1, num_seg_points)]
          this%segments(is)%num_points = size( segment ) + 1
        end if
        this%segments(is)%end_point_idx(2) = this%segments(is)%end_point_idx(1) + this%segments(is)%num_points - 1
        vl = verts(iv)%coord_lat
        call r3frac( eps_zero, vl, vi )
        if( norm2( vl ) < eps_zero ) segment(1) = offset / norm2( vsl )
        do ip = 1, size( segment )
          this%points(np+ip)%coord_lat = verts(iv)%coord_lat + segment(ip) * vsl
          this%points(np+ip)%coord_cart = verts(iv)%coord_cart + segment(ip) * vsc
          this%points(np+ip)%distance = verts(iv)%distance + segment(ip) * this%segments(is)%length
        end do
        this%points(np+1)%distance = verts(iv)%distance
        np = np + size( segment )
        if( verts(iv+1)%break ) then
          iv = iv + 1
          this%points(np)%distance = verts(iv)%distance
        end if
        break_before = verts(iv)%break
        length = length + this%segments(is)%length
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

  !> Return JSON string representing a `vertex` object.
  function vertex_to_json( this ) result( json )
    use xjson, only: to_json
    !> `vertex` object
    class(vertex), intent(in) :: this
    !> JSON string
    character(:), allocatable :: json 

    character(256) :: string

    json = "{"
    write( string, '(a,": ",a)' ) '"coord_lat"', to_json( this%coord_lat )
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"coord_cart"', to_json( this%coord_cart )
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"distance"', to_json( this%distance )
    json = json // trim( string ) // ', '
    if( this%label(1:1) == '\' ) then
      write( string, '(a,": ",a)' ) '"label"', '"\'//trim( this%label )//'"'
    else
      write( string, '(a,": ",a)' ) '"label"', '"'//trim( this%label )//'"'
    end if
    json = json // trim( string ) // ', '
    if( this%break ) then
      write( string, '(a,": ",a)' ) '"break"', 'true'
    else
      write( string, '(a,": ",a)' ) '"break"', 'false'
    end if
    json = json // trim( string ) // '} '
  end function

  !> Return JSON string representing a `point` object.
  function point_to_json( this ) result( json )
    use xjson, only: to_json
    !> `point` object
    class(point), intent(in) :: this
    !> JSON string
    character(:), allocatable :: json 

    character(256) :: string

    json = "{"
    write( string, '(a,": ",a)' ) '"coord_lat"', to_json( this%coord_lat )
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"coord_cart"', to_json( this%coord_cart )
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"distance"', to_json( this%distance )
    json = json // trim( string ) // '} '
  end function

  !> Return JSON string representing a `segment` object.
  function segment_to_json( this ) result( json )
    use xjson, only: to_json
    !> `segment` object
    class(segment), intent(in) :: this
    !> JSON string
    character(:), allocatable :: json 

    character(1024) :: string

    json = "{"
    write( string, '(a,": [",a,", ",a,"]")' ) '"vertices"', this%vertices(1)%to_json(), this%vertices(2)%to_json()
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"end_point_idx"', to_json( this%end_point_idx )
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"num_points"', to_json( this%num_points )
    json = json // trim( string ) // ', '
    write( string, '(a,": ",a)' ) '"length"', to_json( this%length )
    json = json // trim( string ) // '} '
  end function

  !> Return JSON string representing a `bz_path_type` object.
  function path_to_json( this ) result( json )
    use xjson, only: to_json
    !> `bz_path_type` object
    class(bz_path_type), intent(in) :: this
    !> JSON string
    character(:), allocatable :: json 

    integer :: iv, ip, is

    character(1024) :: string

    json = "{"
    ! points
    write( string, '(a,": ",a)' ) '"num_points"', to_json( this%num_points )
    json = json // trim( string ) // ', '
    write( string, '(a,": [")' ) '"points"'
    json = json // trim( string )
    do ip = 1, this%num_points
      json = json // this%points(ip)%to_json()
      if( ip /= this%num_points ) json = json // ', '
    end do
    json = json // '], '
    ! vertices
    write( string, '(a,": ",a)' ) '"num_vertices"', to_json( this%num_vertices )
    json = json // trim( string ) // ', '
    write( string, '(a,": [")' ) '"vertices"'
    json = json // trim( string )
    do iv = 1, this%num_vertices
      json = json // this%vertices(iv)%to_json()
      if( iv /= this%num_vertices ) json = json // ', '
    end do
    json = json // '], '
    ! segments
    write( string, '(a,": ",a)' ) '"num_segments"', to_json( this%num_segments )
    json = json // trim( string ) // ', '
    write( string, '(a,": [")' ) '"segments"'
    json = json // trim( string )
    do is = 1, this%num_segments
      json = json // this%segments(is)%to_json()
      if( is /= this%num_segments ) json = json // ', '
    end do
    json = json // '], '
    ! length
    write( string, '(a,": ",a)' ) '"length"', to_json( this%length )
    json = json // trim( string ) // '} '
  end function

end module bz_path
