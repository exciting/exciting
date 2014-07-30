! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
! -----------------------------------------------------------------------------
subroutine INTERPOLATE_DIEL (diel, e_diel, rlas, diel_intpl)
! -----------------------------------------------------------------------------
!
! reads 7 consecutive values for the dielectric function, interpolates via
! a 4th order polynomial and returns the interpolated value for rlas
!
use raman_params, only: zone
implicit none
! arguments
complex (8), intent(in) :: diel(7)       ! 7 values of the dielectric function
real (8), intent(in) :: e_diel(7)        ! corresponding energy values (assumed ascending)
real (8), intent(in) :: rlas             ! the laser energy to interpolate for
complex (8), intent(out) :: diel_intpl   ! the interpolated value of the dielectric function
! local variables
complex (8) :: A(7, 5), B(7)
complex (8) :: work(8*7)
integer :: i, j, info
real(8) :: ener, de
!
!
!  check validity
if (rlas .lt. e_diel(1) .or. rlas .gt. e_diel(7)) then
   write(*, '("Error(Interpolate_diel): Laser energy not within data range for interpolation")')
   stop
endif
!
!
do i = 1, 7
   A(i,1) = zone
   do j = 2, 5
      A(i,j) = cmplx (e_diel(i)**(j-1), 0.d0, 8)
   enddo
   B(i) = diel(i)
enddo
!
! linear least squares fit
call ZGELS ('N', 7, 5, 1, A, 7, B, 7, work, 8*7, info)
!
!
de = (e_diel(7) - e_diel(1)) / 200.d0
do i = 1, 200
   ener = e_diel(1) + de*dble(i-1)
   diel_intpl = B(1) + B(2)*ener + B(3)*ener**2 + B(4)*ener**3 + B(5)*ener**4
!  if (i .le. 7) then
!     write(111, *) ener, dble(diel_intpl), aimag(diel_intpl), e_diel(i), dble(diel(i)), aimag(diel(i))
!  else
!     write(111, *) ener, dble(diel_intpl), aimag(diel_intpl)
!  endif
enddo
!write(111, *)
!
diel_intpl = B(1) + B(2)*rlas + B(3)*rlas**2 + B(4)*rlas**3 + B(5)*rlas**4
!
return
end subroutine INTERPOLATE_DIEL
