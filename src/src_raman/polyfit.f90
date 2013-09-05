! -----------------------------------------------------------------------------
subroutine POLYFIT (imode)
! -----------------------------------------------------------------------------
!
use modinput
use raman_coeff
implicit none
real(8) :: A(numpot,input%properties%raman%degree+1),B(numpot),work(8*numpot),sing(input%properties%raman%degree+1)
real(8) :: minpot
integer :: i,j,info,m,n, imode
character(84) :: funcout
!
! 
!
write(66,4)
write(funcout,'(a7,7(a5,i1,a4,i1))') 'y = a_0',(' + a_',i,'*x**',i,i=1,input%properties%raman%degree)
write(66,5) funcout
write(66,6) input%properties%raman%degree+1,numpot
!  check validity
if (numpot .lt. input%properties%raman%degree+1) then
   write(*, '("Error(Polyfit): Too few data points given!")')
   stop
endif
!
minpot = 1.d20
do i = 1,numpot
   if (potiny(imode, i) .lt. minpot) minpot = potiny(imode, i)
enddo
!
do i = 1,numpot
   potiny(imode, i) = potiny(imode, i) - minpot           ! shift potential minimum to y(x_min) = 0
   do j = 1,input%properties%raman%degree+1
      A(i,j) = potinx(imode, i)**(j-1)
   enddo
   B(i) = potiny(imode, i)
   write(66,7) potinx(imode, i), potiny(imode, i)
enddo
!
if (input%properties%raman%nstep .le. 2) then
!  do simplicistic determination of quadratic coefficients
   b(1) = 0.d0
   b(2) = 0.d0
   b(3) = potiny(imode, 1) / potinx(imode, 1)**2
else
!  perform linear least squares fit through LAPACK routine DGELS
   m = numpot
   n = input%properties%raman%degree + 1
   call DGELS ('N', m, n, 1, A, m, B, m, work, 8*numpot, info)
   if (info .ne. 0) then
      write(*, '("Warning(Polyfit): Least squares polynomial fit encountered a problem")')
      write(*, '("                  INFO :",i3)') info 
   endif
endif
!
!
a0 = 0.d0; a1 = 0.d0; a2 = 0.d0; a3 = 0.d0; a4 = 0.d0; a5 = 0.d0; a6 = 0.d0
!
!  save solution vector to polynomial coefficients
a0 = B(1)
if (input%properties%raman%degree .ge. 1) a1 = B(2)
if (input%properties%raman%degree .ge. 2) a2 = B(3)
if (input%properties%raman%degree .ge. 3) a3 = B(4)
if (input%properties%raman%degree .ge. 4) a4 = B(5)
if (input%properties%raman%degree .ge. 5) a5 = B(6)
if (input%properties%raman%degree .eq. 6) a6 = B(7)
!
write(66,8) (B(i),i=1,input%properties%raman%degree+1)
!
4 format (/,49('*'),'   DATA FITTING   ',49('*'),' ',/, &
&        54(' '),'potential')
5 format (//,'  Fit function: polynom ',a84)
6 format (/,i8,' coefficients are fitted to ',i4,' data points:',/)
7 format (f12.4,f20.8)
8 format (/,'  Coefficients from polynomial fit a_i: ',/7f15.6)
return
end subroutine polyfit
!
!
