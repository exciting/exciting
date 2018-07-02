
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3ws
! !INTERFACE:
subroutine r3ws( eps, b, v, iv)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in, real)
!   b   : basis vectors (in, real(3,3))
!   v   : input lattice vector (inout, real(3))
!   iv  : lattice vector that maps v in WS-cell (out,integer(3))
! !DESCRIPTION:
!   Finds the lattice vector {\tt iv} that maps the real 3-vector {\tt v}
!   into the Wigner-Seitz cell. On exit, {\tt v} contains a vector within
!   the Wigner-Seitz cell.
!
! !REVISION HISTORY:
!   Created May 2018 (SeTi)
!EOP
!BOC
implicit none
! arguments
real(8), intent( in) :: eps
real(8), intent( in) :: b(3,3)
real(8), intent( inout) :: v(3)
integer, intent( out) :: iv(3)
! local variables
integer :: i, j, k, iv0(3)
real(8) :: d, vl0(3), vl(3), vc(3)

call r3frac( eps, v, iv)
vl0 = v
call r3mv( b, vl0, vc)
d = norm2( vc)

do i = -3, 3
  do j = -3, 3
    do k = -3, 3
      vl = vl0 + dble( (/i, j, k/))
      call r3mv( b, vl, vc)
      if( norm2( vc) .lt. d) then
        d = norm2( vc)
        iv0 = (/i, j, k/)
        v = vl
      end if
    end do
  end do
end do

iv = -iv + iv0

return
end subroutine
!EOC
