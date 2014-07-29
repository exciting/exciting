! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
! -----------------------------------------------------------------------------
subroutine POLYFIT (imode)
! -----------------------------------------------------------------------------
!
use modinput
use raman_coeff
implicit none
real(8) :: A(numpot,input%properties%raman%degree+1),B(numpot),work(8*numpot)
real(8) :: minpot
integer :: i,j,info,m,n, imode
character(84) :: funcout
!
! 
write(66,4)
!
if (input%properties%raman%useforces) then
 write(66, '(/," Fit from forces: ")')
 write(funcout,'(a8,7(a3,i1,a2,i1,a4,i1))') 'y = -a_1',(' - ',i,'a_',i,'*x**',i-1,i=2,input%properties%raman%degree)
 write(66,5) funcout
 write(66,6) input%properties%raman%degree,numpot
 !  check validity
 if (numpot .lt. input%properties%raman%degree) then
    write(*, '("Error(Polyfit): Too few data points given!")')
    stop
 endif
else
 write(66, '(/," Fit from total energies: ")')
 write(funcout,'(a7,7(a5,i1,a4,i1))') 'y = a_0',(' + a_',i,'*x**',i,i=1,input%properties%raman%degree)
 write(66,5) funcout
 write(66,6) input%properties%raman%degree+1,numpot
 !  check validity
 if (numpot .lt. input%properties%raman%degree+1) then
    write(*, '("Error(Polyfit): Too few data points given!")')
    stop
 endif
endif
!
if (.not. input%properties%raman%useforces) then
   ! shift potential minimum to y(x_min) = 0
   minpot = 1.d20
   do i = 1,numpot
      if (potiny(imode, i) .lt. minpot) minpot = potiny(imode, i)
   enddo
   do i = 1,numpot
      potiny(imode, i) = potiny(imode, i) - minpot
      do j = 1,input%properties%raman%degree+1
         A(i,j) = potinx(imode, i)**(j-1)
      enddo
      B(i) = potiny(imode, i)
      write(66,7) potinx(imode, i), potiny(imode, i)
   enddo
   n = input%properties%raman%degree + 1
else
   do i = 1, numpot
      do j = 1,input%properties%raman%degree
         A(i,j) = -dble(j)*potinx(imode, i)**(j-1)
      enddo
      B(i) = potiny(imode, i)
      write(66,7) potinx(imode, i), potiny(imode, i)
   enddo
   n = input%properties%raman%degree
endif
!
!  perform linear least squares fit through LAPACK routine DGELS
m = numpot
call DGELS ('N', m, n, 1, A, m, B, m, work, 8*numpot, info)
if (info .ne. 0) then
   write(*, '("Warning(Polyfit): Least squares polynomial fit encountered a problem")')
   write(*, '("                  INFO :",i3)') info 
endif
!
!
a0 = 0.d0; a1 = 0.d0; a2 = 0.d0; a3 = 0.d0; a4 = 0.d0; a5 = 0.d0; a6 = 0.d0
!
!  save solution vector to polynomial coefficients
if (.not. input%properties%raman%useforces) then
   ! for total energies
   a0 = B(1)
   if (input%properties%raman%degree .ge. 1) a1 = B(2)
   if (input%properties%raman%degree .ge. 2) a2 = B(3)
   if (input%properties%raman%degree .ge. 3) a3 = B(4)
   if (input%properties%raman%degree .ge. 4) a4 = B(5)
   if (input%properties%raman%degree .ge. 5) a5 = B(6)
   if (input%properties%raman%degree .eq. 6) a6 = B(7)
   write(66,8) (B(i),i=1,input%properties%raman%degree+1)
else
   ! for forces
   a0 = 0.d0
   if (input%properties%raman%degree .ge. 1) a1 = B(1)
   if (input%properties%raman%degree .ge. 2) a2 = B(2)
   if (input%properties%raman%degree .ge. 3) a3 = B(3)
   if (input%properties%raman%degree .ge. 4) a4 = B(4)
   if (input%properties%raman%degree .ge. 5) a5 = B(4)
   if (input%properties%raman%degree .eq. 6) a6 = B(6)
   write(66,8) (B(i),i=1,input%properties%raman%degree)
endif
!
!
4 format (/,49('*'),'   DATA FITTING   ',49('*'),' ',/, &
&        54(' '),'potential')
5 format (/,'  Fit function: polynom ',a84)
6 format (/,i8,' coefficients are fitted to ',i4,' data points:',/)
7 format (f12.4,f20.8)
8 format (/,'  Coefficients from polynomial fit a_i: ',/7f15.6)
return
end subroutine polyfit
!
!
