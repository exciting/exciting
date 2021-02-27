


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine leblaik(n, v, w)
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(8), intent(out) :: v(3, n), w(n)
  ! local variables
  integer :: m
  real(8), allocatable :: x(:), y(:), z(:)
  if (n.lt.1) then
     write(*, *)
     write(*, '("Error(leblaik): n < 1")')
     write(*, *)
     stop
  end if
  m=n
  allocate(x(n), y(n), z(n)) 
  select case(n)
  case(6)
    call ld0006(x, y, z, w, m)
  case(14)
    call ld0014(x, y, z, w, m)
  case(26)
    call ld0026(x, y, z, w, m)
  case(38)
    call ld0038(x, y, z, w, m)
  case(50)
    call ld0050(x, y, z, w, m)
  case(74)
    call ld0074(x, y, z, w, m)
  case(86)
    call ld0086(x, y, z, w, m)
  case(110)
    call ld0110(x, y, z, w, m)
  case(146)
    call ld0146(x, y, z, w, m)
  case(170)
    call ld0170(x, y, z, w, m)
  case(194)
    call ld0194(x, y, z, w, m)
  case(230)
    call ld0230(x, y, z, w, m)
  case(266)
    call ld0266(x, y, z, w, m)
  case(302)
    call ld0302(x, y, z, w, m)
  case(350)
    call ld0350(x, y, z, w, m)
  case(434)
    call ld0434(x, y, z, w, m)
  case(590)
    call ld0590(x, y, z, w, m)
  case(770)
    call ld0770(x, y, z, w, m)
  case(974)
    call ld0974(x, y, z, w, m)
  case(1202)
    call ld1202(x, y, z, w, m)
  case(1454)
    call ld1454(x, y, z, w, m)
  case(1730)
    call ld1730(x, y, z, w, m)
  case(2030)
    call ld2030(x, y, z, w, m)
  case(2354)
    call ld2354(x, y, z, w, m)
  case(2702)
    call ld2702(x, y, z, w, m)
  case(3074)
    call ld3074(x, y, z, w, m)
  case(3470)
    call ld3470(x, y, z, w, m)
  case(3890)
    call ld3890(x, y, z, w, m)
  case(4334)
    call ld4334(x, y, z, w, m)
  case(4802)
    call ld4802(x, y, z, w, m)
  case(5294)
    call ld5294(x, y, z, w, m)
  case(5810)
    call ld5810(x, y, z, w, m)
  case default
    write(*, *)
    write(*, '("Error(leblaik): invalid number of covering points: ", i6)') m
    write(*, '(" possible numbers are: 6, 14, 26, 38, 50, 74, 86, 110, 146, ")')
    write(*, '(" 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, ")')
    write(*, '(" 1454, 1730, 2030, 2354, 2702, 3074, 3740, 3890, 4334, 4802, ")')
    write(*, '(" 5294, 5810")')
    write(*, *)
    stop
  end select
  v(1, :)=x(:)
  v(2, :)=y(:)
  v(3, :)=z(:)
  deallocate(x, y, z)
end subroutine leblaik
