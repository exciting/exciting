! ------------------------------------------------------------------------------
subroutine SPECTRUM_M(temper, rlas, filext)
! ------------------------------------------------------------------------------
use raman_ew
use raman_input
use raman_coeff
use modinput
use mod_lattice, only: omega
implicit none
! arguments
real (8), intent(in) :: temper, rlas
character(*), intent(in) :: filext
! local variables
integer :: i,j,ind,ii,iw
integer :: oct1, oct2
real(8) :: arg,beta,dde,factl,gam1,gam2,gam3,gam4,dnc
real(8) :: sc,sqn1,sqn2,sqn3
real(8) :: sqn4,sqn5,sqn6
real(8) :: spec,w,dw,zfact,zsum,zzsum, temp, w_scat, w_fact
real(8) :: isoinvar, anisoinvar2, activity, depol
character(80) :: fname
real(8), allocatable :: rcont_par(:), rcont_perp(:)
real(8), allocatable :: ex(:)
complex(8) :: trans1!,trans2,trans3,trans4,trans5,trans6
!
allocate( ex(input%properties%raman%nstate) )
allocate( rcont_par(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
allocate( rcont_perp(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!allocate( rcont2(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!allocate( rcont3(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!allocate( rcont4(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!allocate( rcont5(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!allocate( rcont6(input%properties%raman%nstate*(input%properties%raman%nstate-1)/2) )
!
rcont_par = 0.d0; rcont_perp = 0.d0
!
fname = 'RAMAN_SPEC_OC'//filext
open (unit=73, file=fname, status='unknown',form='formatted')
!
! scaling factor for output of scattering contributions
sc = 1.d32*fau3cm3**2
!
!dnc = dble(ncell)
dnc = 1.d0
write(66,'(//,51("*"),"   SPECTRUM   ",51("*"))')
temp = temper
write(66,'(/," Temperature [ K ]: ",17x,f8.2,/)') temp
if (temp .lt. 1.0e-9) then
   temp = 1.d-9
   write(66,'(  " (adjusted to     : ",17x,e9.2,")"/)') temp
endif
write(66,'(/," Partition sums for eigenstates:",/)')  ! header for part sums in output
write(66,'("  state ","         eigenvalue ","   Boltzmann factor ", &
     &           "      partition sum ","         population  "/)')
!
gam1 = gamma1**2/4.d0                         ! broadening parameters from input
gam2 = gamma1**2/4.d0
gam3 = gamma1**2/4.d0
gam4 = gamma1**2/4.d0
zsum = 0.d0
beta = 1.d0/fkwn/temp
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
dw = (input%properties%raman%energywindow%intv(2) - &
   &  input%properties%raman%energywindow%intv(1)) / &
   &  dble(input%properties%raman%energywindow%points)
write(73,'("# Laser energy:  ",f9.6," Ha = ",f8.1," cm^-1 = ",f5.2," eV = ",  &
  &        f6.1, " nm",/,"# Temperature: ",f8.1," K",/,"# Damping: ",f8.1," cm^-1   ",/, &
  &       "# Number of eigenvalues: ",i4,/)') rlas,rlas*fhawn,rlas*fhaev,frnmha/rlas,temp, &
  &              gamma1,input%properties%raman%nstate
!
isoinvar = 1.d0/3.d0 * (deq(1,1) + deq(2,2) + deq(3,3))
anisoinvar2 = 0.5 * ((deq(1,1)-deq(2,2))**2 + (deq(2,2)-deq(3,3))**2 + (deq(3,3)-deq(1,1))**2) + &
  & 3.d0 * (deq(1,2)**2 + deq(1,3)**2 + deq(2,3)**2)
activity = (45.d0*abs(isoinvar)**2 + 7.d0*abs(anisoinvar2)) / fgew
depol = 3.d0 * abs(anisoinvar2) / (45.d0*abs(isoinvar)**2 + 4.d0*abs(anisoinvar2))
write(66, '(/," Raman activity [ A^4 amu^-1 ] : ",f14.6,/, &
       &      " depolarization                : ",f14.6,/)') activity*faua**4*famuau, depol
!
ind = 0
! drop V, divide by epsilon_0 = 1/4pi (au)
zzsum = zsum/(2.d0*pi)**4
write(66,'(/," Scattering contributions in [ 10^-32 cm^2 sr^-1 ] " &
  & ,//,4x," Transition energy [ cm^-1 ]",            &
  & 6x,"parallel",8x,"perpendicular",13x,"absolute",/)')
do i = 1,input%properties%raman%nstate
   do j = i-1,1,-1
      dde = (eigen(i) - eigen(j))*fhawn
      w_scat = rlas*fhawn - dde              ! Stokes scattering
      w_fact = rlas*fhawn*w_scat**3
      zfact = ex(j)/zzsum                    ! prefactors * exp(-Ej(kT) / Sum( exp(E/kT) )
      ind = ind + 1
!
!     trans_iso = isoinvar*transme1(ind)/factorial(1)/sqn1  !+ &
!    &               d2eq(oct1,oct2)*transme2(ind)/factorial(2)/sqn2 + &
!    &               d3eq(oct1,oct2)*transme3(ind)/factorial(3)/sqn3 + &
!    &               d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4 + &
!    &               d5eq(oct1,oct2)*transme5(ind)/factorial(5)/sqn5 + &
!    &               d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
!     if (input%properties%raman%intphonon) then
!        trans2 = d2eq(oct1,oct2)*transme2(ind)/factorial(2)/sqn2 + &
!    &            d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4 + &
!    &            d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
!        trans3 = d3eq(oct1,oct2)*transme3(ind)/factorial(3)/sqn3 + &
!    &            d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
!        trans4 = d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4
!        trans5 = d5eq(oct1,oct2)*transme5(ind)/factorial(5)/sqn5 
!        trans6 = d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
!     endif
!
      rcont_par(ind)  = (45.d0*abs(isoinvar)**2 + 4.d0*abs(anisoinvar2)) / 45.d0 &
     &                  * (transme1(ind)/factorial(1)/sqn1)**2                   &
     &                  * dnc*zfact
      rcont_perp(ind) = 3.d0*abs(anisoinvar2) / 45.d0                            &
     &                  * (transme1(ind)/factorial(1)/sqn1)**2                   &
     &                  * dnc*zfact
!     rcont1(ind) = dnc*zfact*              abs(trans1)**2
!     if (input%properties%raman%intphonon) then
!        rcont2(ind) = dnc*zfact*(dnc    - 1.)*abs(trans2)**2
!        rcont3(ind) = dnc*zfact*(dnc**2 - 1.)*abs(trans3)**2
!        rcont4(ind) = dnc*zfact*(dnc**3 - 1.)*abs(trans4)**2
!        rcont5(ind) = dnc*zfact*(dnc**4 - 1.)*abs(trans5)**2
!        rcont6(ind) = dnc*zfact*(dnc**5 - 1.)*abs(trans6)**2
!     endif
!
!
!     if (input%properties%raman%intphonon) then
!        write(66,'(1x,I3," ->",I3,3x,f8.2,5x,6(f14.8,"  |"))') j,i,dde,             &
!    &           rcont1(ind)*sc*w_fact,rcont2(ind)*sc*w_fact, &
!    &           rcont3(ind)*sc*w_fact,rcont4(ind)*sc*w_fact, &
!    &           rcont5(ind)*sc*w_fact,rcont6(ind)*sc*w_fact
!     else
         write(66,'(3x,I3," ->",I3,5x,f8.2,3(7x,f14.8))') j,i,dde,             &
     &           rcont_par(ind)*sc*w_fact, rcont_perp(ind)*sc*w_fact,          &
     &          (rcont_par(ind)+rcont_perp(ind))*sc*w_fact
!     endif
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
   w = input%properties%raman%energywindow%intv(1) + dble(iw - 1)*dw
   do i = 1,input%properties%raman%nstate
      do j = i-1,1,-1
         ind = ind + 1
         dde = (eigen(i) - eigen(j))*fhawn
         spec = spec + &                        ! line broadening by Lorentzian
     &          (rcont_par(ind)+rcont_perp(ind))*gamma1/((dde-w)**2 + gam1)!+ &
!    &          rcont2(ind)*gamma2/((dde-w)**2 + gam2) + &
!    &          rcont3(ind)*gamma3/((dde-w)**2 + gam3) + &
!    &          rcont4(ind)*gamma4/((dde-w)**2 + gam4) + &
!    &          rcont5(ind)*gamma4/((dde-w)**2 + gam4) + &   ! 5th and 6th use gamma4
!    &          rcont6(ind)*gamma4/((dde-w)**2 + gam4)
      enddo
   enddo
   spec = spec*(rlas*fhawn - w)**3/factl/(speed_of_light*100.d0)
   write(73,'(2g15.8)') w,spec
enddo
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
