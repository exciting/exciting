! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! adapted from spectrum.f90, July 2014 by Stefan Kontur
!
! ------------------------------------------------------------------------------
subroutine SPECTRUM_M(rlas, filext)
! ------------------------------------------------------------------------------
use raman_ew
use raman_coeff
use modinput
use mod_lattice, only: omega
implicit none
! arguments
real (8), intent(in) :: rlas
character(*), intent(in) :: filext
! local variables
integer :: i,j,ind,ii,iw
integer :: oct1, oct2
real(8) :: arg,beta,dde,gam2
real(8) :: sc, total
real(8) :: spec_par, spec_perp, spec_abs
real(8) :: w, dw, zfact, zsum, zzsum, w_scat, w_fact
real(8) :: anisoinvar2, activity, depol
complex(8) :: isoinvar
character(80) :: fname
real(8), allocatable :: rcont_par(:), rcont_perp(:)
real(8), allocatable :: ex(:)
complex(8) :: trans1
!
allocate( ex(input%properties%raman%nstate) )
allocate( rcont_par(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont_perp(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!
rcont_par = 0.d0; rcont_perp = 0.d0
!
fname = 'RAMAN_SPEC_OC'//filext
open (unit=73, file=fname, status='unknown',form='formatted')
!
! scaling factor for output of scattering contributions
sc = 1.d32*fau3cm3**2
!
write(66,'(//,51("*"),"   SPECTRUM   ",51("*"))')
write(66,'(/," Temperature [ K ]: ",17x,f8.2,/)') input%properties%raman%temp
if (input%properties%raman%temp .lt. 1.0e-9) then
   input%properties%raman%temp = 1.d-9
   write(66,'(  " (adjusted to     : ",17x,e9.2,")"/)') input%properties%raman%temp
endif
write(66,'(/," Partition sums for eigenstates:",/)')  ! header for part sums in output
write(66,'("  state ","         eigenvalue ","             factor ", &
     &           "      partition sum ","         population  "/)')
!
gam2 = input%properties%raman%broad**2/4.d0
zsum = 0.d0
beta = 1.d0/fkwn/input%properties%raman%temp
do i = 1,input%properties%raman%nstate
   arg = -beta*eigen(i)*fhawn
   ex(i) = dexp(arg)
   zsum = zsum + ex(i)                        ! zsum = Sum( exp(-E/kT) )
enddo
if (ex(1) .eq. 0.d0) then                     ! handle numerical problems with low temperatures
   ex(1) = 1.d0
   zsum = 1.d0
endif
write(66,'(i7,3e20.6,f20.5)') (i,eigen(i),ex(i),zsum,ex(i)/zsum,i=1,input%properties%raman%nstate)
dw = (input%properties%raman%energywindow%intv(2) - &
   &  input%properties%raman%energywindow%intv(1))*fhawn / &
   &  dble(input%properties%raman%energywindow%points)
write(73,'("# Laser energy:  ",f9.6," Ha = ",f8.1," cm^-1 = ",f5.2," eV = ",  &
  &        f6.1, " nm",/,"# Temperature: ",f8.1," K",/,"# Damping: ",f8.1," cm^-1   ",/, &
  &       "# Number of eigenvalues: ",i4,/)') rlas,rlas*fhawn,rlas*fhaev, &
  &              frnmha/rlas,input%properties%raman%temp, &
  &              input%properties%raman%broad,input%properties%raman%nstate
!
isoinvar = 1.d0/3.d0 * (deq(1,1) + deq(2,2) + deq(3,3))
anisoinvar2 = 0.5 * (abs(deq(1,1)-deq(2,2))**2 + abs(deq(2,2)-deq(3,3))**2 + abs(deq(3,3)-deq(1,1))**2) + &
  & 3.d0 * (abs(deq(1,2))**2 + abs(deq(1,3))**2 + abs(deq(2,3))**2)
activity = (45.d0*abs(isoinvar)**2 + 7.d0*anisoinvar2) / fgew
depol = 3.d0 * anisoinvar2 / (45.d0*abs(isoinvar)**2 + 4.d0*anisoinvar2)
write(66, '(//," Raman activity [ a.u. ]       : ",f14.6,/, &
       &       " Raman activity [ A^4 amu^-1 ] : ",f14.6,/, &
       &       " depolarization                : ",f14.6,/)') activity, activity*faua**4*famuau, depol
!
ind = 0
! drop V, divide by epsilon_0 = 1/4pi (au)
zzsum = zsum/(2.d0*pi)**4
write(66,'(/," Differential cross section (d sigma) / (d Omega), in [ 10^-32 cm^2 sr^-1 ] " &
  & ,//,2x,"Transition",4x,"energy [ cm^-1 ]",15x,"polarization",/,            &
  & 38x,"parallel",8x,"perpendicular",13x,"absolute",/)')
do i = 1,input%properties%raman%nstate
   do j = i-1,1,-1
      dde = (eigen(i) - eigen(j))*fhawn
      w_scat = rlas*fhawn - dde              ! Stokes scattering
      w_fact = rlas*fhawn*w_scat**3
      zfact = ex(j)/zzsum                    ! prefactors * exp(-Ej(kT) / Sum( exp(E/kT) )
      ind = ind + 1
!
      rcont_par(ind)  = (45.d0*abs(isoinvar)**2 + 4.d0*anisoinvar2) / 45.d0 &
     &                  * (transme1(ind)/factorial(1))**2                   &
     &                  * zfact
      rcont_perp(ind) = 3.d0*anisoinvar2 / 45.d0                            &
     &                  * (transme1(ind)/factorial(1))**2                   &
     &                  * zfact
      write(66,'(3x,I3," ->",I3,5x,f8.2,3(7x,f14.8))') j,i,dde,             &
     &        rcont_par(ind)*sc*w_fact, rcont_perp(ind)*sc*w_fact,          &
     &       (rcont_par(ind)+rcont_perp(ind))*sc*w_fact
   enddo
enddo
!
! write out spectra
total = 0.d0
do iw = 1, input%properties%raman%energywindow%points+1
   ind = 0
   spec_par = 0.0d0
   spec_perp = 0.d0
   w = input%properties%raman%energywindow%intv(1)*fhawn + dble(iw - 1)*dw
   do i = 1,input%properties%raman%nstate
      do j = i-1,1,-1
         ind = ind + 1
         dde = (eigen(i) - eigen(j))*fhawn
         w_scat = rlas*fhawn - dde
         w_fact = rlas*fhawn*w_scat**3
         spec_par = spec_par + rcont_par(ind)*w_fact/((dde-w)**2 + gam2)
         spec_perp = spec_perp + rcont_perp(ind)*w_fact/((dde-w)**2 + gam2)
      enddo
   enddo
!   spec_par = spec_par*(rlas*fhawn - w)**3/(2.d0*pi)
!   spec_perp = spec_perp*(rlas*fhawn - w)**3/(2.d0*pi)
   spec_par = spec_par*sc*input%properties%raman%broad/(2.d0*pi)
   spec_perp = spec_perp*sc*input%properties%raman%broad/(2.d0*pi)
   spec_abs = spec_par + spec_perp
   total = total + spec_abs*dw
   write(73,'(f12.2,3g15.8)') w, spec_par, spec_perp, spec_abs
enddo
write(*,*) 'total ',total
!
close(73)
!
deallocate( ex )
deallocate( rcont_par, rcont_perp )
!
return
!
end
! 
