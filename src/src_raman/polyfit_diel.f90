! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
! -----------------------------------------------------------------------------
subroutine POLYFIT_DIEL (imode, rlas)
! -----------------------------------------------------------------------------
!
use modinput
use raman_coeff
use raman_params, only: zzero
implicit none
! arguments
integer, intent(in) :: imode
real (8), intent(in) :: rlas
! local variables
complex (8) :: A(numpot,input%properties%raman%degree+1) 
complex (8) :: B(numpot)
complex (8) :: work(8*numpot)
complex (8), allocatable :: dieldat(:, :, :)
complex (8) :: diel_intpl, eps_calc
real (8) :: t1, e_diel(input%xs%energywindow%points), xx
character(88) :: funcout
integer :: i, j, k, m, n, iw, oct1, oct2, deg, info
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
if (input%properties%raman%molecule) then
   write(66,'(/,49("*"),"   DATA FITTING   ",49("*")," ",/,51(" "),"polarizability")')
   write(funcout,'(a11,7(a5,i1,a4,i1))') 'alpha = b_0',(' + b_',i,'*u**',i,i=1,deg)
else
   write(66,'(/,49("*"),"   DATA FITTING   ",49("*")," ",/,48(" "),"dielectric function")')
   write(funcout,'(a9,7(a5,i1,a4,i1))') 'eps = b_0',(' + b_',i,'*u**',i,i=1,deg)
endif
write(66,'(//,"  Fit function: polynom ",a84)') funcout
write(66,'(/,i8," coefficients are fitted to ",i4," datapoints:",/)') deg+1,numpot
!
!  check validity
if (numpot .lt. deg+1) then
   write(*, '("Error(Polyfit_diel): Too few data points given!")')
   stop
endif
!
! determine data point of df
t1 = (input%xs%energywindow%intv(2) - input%xs%energywindow%intv(1)) / dble( input%xs%energywindow%points )
iw = nint((rlas - input%xs%energywindow%intv(1)) / t1) + 1
do i = 1, input%xs%energywindow%points
   e_diel(i) = input%xs%energywindow%intv(1) + t1*dble(i-1)
enddo
!
deq = zzero; d2eq = zzero; d3eq = zzero; d4eq = zzero; d5eq = zzero; d6eq = zzero
!
Do oct1 = 1, 3
   Do oct2 = 1, 3
      if (oct1 .ne. oct2 .and. .not.offdiag) cycle
      write(66,'(/,"  Fitting to data for optical component ",a2,/)') comp(oct1, oct2)
!   interpolate for rlas
      write(66,'("  Values are interpolated for laser energy: ",f12.4," eV",/)') rlas*fhaev
      if (input%properties%raman%molecule) then
         write(66,'("    u [Bohr]   Re alpha [Bohr^3]   Im alpha [Bohr^3]")')
      else
         write(66,'("    u [Bohr]          Re epsilon          Im epsilon")')
      endif
      do i = 1, numpot
         if ((iw - 3) .le. 0) iw = 4 
         if ((iw + 3) .gt. input%xs%energywindow%points) iw = input%xs%energywindow%points - 3
         call interpolate_diel ( df (imode, i, oct1, oct2, (iw-3):(iw+3)),  &
           &                     e_diel ((iw-3):(iw+3)), &
           &                     rlas, &
           &                     dieldat (oct1, oct2, i) )
         write(66,'(f12.4,2f20.8)') potinx(imode, i), dble( dieldat(oct1, oct2, i) ), aimag( dieldat(oct1, oct2, i) )
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
! write results for each component to output file
      write(66,'(/,"  Coefficients from polymomial fit b_i (real part): ",/5f20.6)') &
          &   ( dble(B(i)),i=1,deg+1)
      write(66,'(/,"  Coefficients from polymomial fit b_i (imaginary part): ",/5f20.6,/)') &
          &   (aimag(B(i)),i=1,deg+1)
!
! write fit functions for plotting
      write(str_mode,'(i3.3)') imode
      open(unit=80,file='RAMAN_EPSILON_OC'//comp(oct1,oct2)//'_MOD'//str_mode//'.OUT',status='unknown',form='formatted')
      write(80,'(/,"# Fitted function and raw data (Re and Im) for optical component ",a2)') comp(oct1, oct2)
      do i = 1,input%properties%raman%ninter+1
         xx = input%properties%raman%xmin + (i - 1)*(input%properties%raman%xmax - &
            & input%properties%raman%xmin)/input%properties%raman%ninter
         eps_calc = zzero
         do k = 1, deg + 1
            eps_calc = eps_calc + B(k)*xx**(k-1)
         enddo
         if (i .le. numpot) then
            write(80,'(2(f12.6,2f26.8))') xx, dble(eps_calc), aimag(eps_calc), potinx(imode, i), &
             &                  dble(dieldat(oct1, oct2, i)), aimag(dieldat(oct1, oct2, i))
         else
            write(80,'(f12.6,2f26.8)')    xx, dble(eps_calc), aimag(eps_calc)
         endif
      enddo
      close (80)
! end loops over optical components
   enddo
enddo
!
deallocate( dieldat )
!
return
end subroutine POLYFIT_DIEL
