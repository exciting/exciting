! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014, Stefan Kontur
!
! ------------------------------------------------------------------------------
subroutine SPECTRUM(rlas, oct1, oct2, comp_ch, filext)
! ------------------------------------------------------------------------------
use raman_ew
use raman_coeff
use modinput
use mod_lattice, only: omega
implicit none
! arguments
real (8), intent(in) :: rlas
integer, intent(in) :: oct1, oct2
character(2), intent(in) :: comp_ch
character(*), intent(in) :: filext
! local variables
integer :: i,j,ind,iw
real(8) :: arg,beta,dde,factl,gam2,dnc
real(8) :: sc,sqn1,sqn2,sqn3
real(8) :: sqn4,sqn5,sqn6
real(8) :: spec,w,dw,zfact,zsum,zzsum, w_scat, w_fact
character(80) :: fname
real(8), allocatable :: rcont1(:),rcont2(:),rcont3(:),rcont4(:),rcont5(:),rcont6(:)
real(8), allocatable :: ex(:)
complex(8) :: trans1,trans2,trans3,trans4,trans5,trans6
!
allocate( ex(input%properties%raman%nstate) )
allocate( rcont1(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont2(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont3(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont4(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont5(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont6(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!
rcont1 = 0.d0; rcont2 = 0.d0; rcont3 = 0.d0; rcont4 = 0.d0; rcont5 = 0.d0; rcont6 = 0.d0
!
fname = 'RAMAN_SPEC_OC'//comp_ch//filext
open (unit=73, file=fname, status='unknown',form='formatted')
!
! scaling factor for output of scattering contributions
   sc = 1.d07
dnc = dble(ncell)
write(66,'(//,46("*"),"   SPECTRUM FOR OC ",a2,"   ",46("*"))') comp(oct1, oct2)
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
sqn1 = dsqrt(dnc)                             ! powers of sqrt( Nc )
sqn2 = dnc
sqn3 = sqn1*sqn2
sqn4 = sqn2*sqn2
sqn5 = sqn2*sqn3
sqn6 = sqn3*sqn3
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
! energy window for Raman spectrum given in Ha, spectrum computed in wave numbers
dw = (input%properties%raman%energywindow%intv(2) - &
   &  input%properties%raman%energywindow%intv(1))*fhawn / &
   &  dble(input%properties%raman%energywindow%points)
write(73,'("# Laser energy:  ",f9.6," Ha = ",f8.1," cm^-1 = ",f5.2," eV = ",  &
  &        f6.1, " nm",/,"# Temperature: ",f8.1," K",/,"# Damping: ",f8.1," cm^-1   ",/, &
  &       "# Number of eigenvalues: ",i4,/)') rlas,rlas*fhawn,rlas*fhaev, &
  &              frnmha/rlas,input%properties%raman%temp, &
  &              input%properties%raman%broad,input%properties%raman%nstate
!
ind = 0
   zzsum = zsum/(omega*fau3cm3*pi*pi)
   write(66,'(/," Scattering contributions in [ 10^-5 sr^-1 m^-1 ]",// &
     & ,4x," Transition energy [ cm^-1 ]",/)')
do i = 1,input%properties%raman%nstate
   do j = i-1,1,-1
      dde = (eigen(i) - eigen(j))*fhawn
      w_scat = rlas*fhawn - dde              ! Stokes scattering
      w_fact = rlas*fhawn*w_scat**3
      zfact = ex(j)/zzsum                    ! prefactors * exp(-Ej(kT) / Sum( exp(E/kT) )
      ind = ind + 1
!
      trans1 = deq (oct1,oct2)*transme1(ind)/factorial(1)/sqn1 + &
     &         d2eq(oct1,oct2)*transme2(ind)/factorial(2)/sqn2 + &
     &         d3eq(oct1,oct2)*transme3(ind)/factorial(3)/sqn3 + &
     &         d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4 + &
     &         d5eq(oct1,oct2)*transme5(ind)/factorial(5)/sqn5 + &
     &         d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
!
      rcont1(ind) = dnc*zfact*              abs(trans1)**2
!
!
      write(66,'(3x,I3," ->",I3,5x,f8.2,7x,f14.8)') j,i,dde,             &
     &           rcont1(ind)*sc*w_fact
   enddo
enddo
!
! write out spectra
!factl = 2.d0*pi*(rlas*fhawn)**3
!factl = 2.d0*pi*(rlas*fhawn)**4
factl = 4.d0*pi*pi
do iw = 1, input%properties%raman%energywindow%points+1
   ind = 0
   spec = 0.0d0
   w = input%properties%raman%energywindow%intv(1)*fhawn + dble(iw - 1)*dw
   do i = 1,input%properties%raman%nstate
      do j = i-1,1,-1
         ind = ind + 1
         dde = (eigen(i) - eigen(j))*fhawn
         w_scat = rlas*fhawn - dde
         w_fact = rlas*fhawn*w_scat**3
         spec = spec + &                        ! line broadening by Lorentzian
     &          rcont1(ind)*w_fact/((dde-w)**2 + gam2) + &
     &          rcont2(ind)*w_fact/((dde-w)**2 + gam2) + &
     &          rcont3(ind)*w_fact/((dde-w)**2 + gam2) + &
     &          rcont4(ind)*w_fact/((dde-w)**2 + gam2) + &
     &          rcont5(ind)*w_fact/((dde-w)**2 + gam2) + &
     &          rcont6(ind)*w_fact/((dde-w)**2 + gam2)
      enddo
   enddo
   spec = spec*input%properties%raman%broad*sc/(2.d0*pi)
   write(73,'(2g15.8)') w,spec
enddo
!
close(73)
!
deallocate( ex )
deallocate( rcont1,rcont2,rcont3,rcont4,rcont5,rcont6 )
!
return
!
end
! 
