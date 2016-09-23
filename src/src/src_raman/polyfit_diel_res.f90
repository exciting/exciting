! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
! -----------------------------------------------------------------------------
subroutine POLYFIT_DIEL_RES (imode, iw)
! -----------------------------------------------------------------------------
!
use modinput
use raman_coeff
use raman_params, only: zzero
implicit none
! arguments
integer, intent(in) :: imode
integer, intent(in) :: iw
! local variables
complex (8) :: A(numpot,input%properties%raman%degree+1) 
complex (8) :: B(numpot)
complex (8) :: work(8*numpot)
complex (8), allocatable :: dieldat(:, :, :)
character(88) :: funcout
integer :: i, j, k, m, n, oct1, oct2, deg, info
character(3) :: str_mode
!
allocate( dieldat(3, 3, numpot) )
!
if (input%properties%raman%nstep .le. 2) then
   deg = 1
else
   deg = input%properties%raman%degree
endif
!deg = 1
!
!
!  check validity
if (numpot .lt. deg+1) then
   write(*, '("Error(Polyfit_diel_silent): Too few data points given!")')
   stop
endif
!
!
deq = zzero; d2eq = zzero; d3eq = zzero; d4eq = zzero; d5eq = zzero; d6eq = zzero
!
Do oct1 = 1, 3
   Do oct2 = 1, 3
      if (oct1 .ne. oct2 .and. .not.offdiag) cycle
      do i = 1, numpot
         dieldat (oct1, oct2, i) = df (imode, i, oct1, oct2, iw)
      enddo
!
      if (input%properties%raman%nstep .le. 2) then
         ! 2 is the equilibirum geometry result
         B(1) =  dieldat(oct1, oct2, 2)
         B(2) = (dieldat(oct1, oct2, 3) - dieldat(oct1, oct2, 2)) / potinx(imode, 3)
      else
         do i = 1,numpot
            A(i,1) = zone
            do j = 2, deg + 1
               A(i,j) = cmplx(potinx(imode, i)**(j-1), 0.d0, 8)
            enddo
            B(i) = dieldat(oct1, oct2, i)
         enddo
!
!  perform linear least squares fit through LAPACK routine ZGELS
         m = numpot
         n = deg + 1
         call ZGELS ('N', m, n, 1, A, m, B, m, work, 8*numpot, info)
      endif
!
!  save solution vector to polynomial coefficients
!  d*eq are assumed to be derivatives; b_k = f^(k)(0) / k!
!  take first order only, leave higher order 0
      if (deg .ge. 1) deq  (oct1, oct2) =       B(2)
!     if (deg .ge. 2) d2eq (oct1, oct2) =  2.d0*B(3)
!     if (deg .ge. 3) d3eq (oct1, oct2) =  6.d0*B(4)
!     if (deg .ge. 4) d4eq (oct1, oct2) = 24.d0*B(5)
!     if (deg .ge. 5) d5eq (oct1, oct2) = 12.d1*B(6)
!     if (deg .eq. 6) d6eq (oct1, oct2) = 72.d1*B(7)
!
! end loops over optical components
   enddo
enddo
!
deallocate( dieldat )
!
return
end subroutine POLYFIT_DIEL_RES
