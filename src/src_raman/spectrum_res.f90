! ------------------------------------------------------------------------------
subroutine SPECTRUM_RES(temper, rlas, oct1, oct2, Sab, ws)
! ------------------------------------------------------------------------------
use raman_ew
use raman_input
use raman_coeff
use modinput
use mod_lattice, only: omega
implicit none
! arguments
real (8), intent(in) :: temper, rlas
integer, intent(in) :: oct1, oct2
real(8), intent(out) :: Sab
real(8), intent(out) :: ws
! local variables
integer :: i,j,ind
real(8) :: arg,beta,dde,factl,gam1,gam2,gam3,gam4,dnc
real(8) :: sc,sqn1,sqn2,sqn3
real(8) :: sqn4,sqn5,sqn6
real(8) :: spec,zfact,zsum,zzsum, temp, w_scat, w_fact
real(8) :: rcont1,rcont2,rcont3,rcont4,rcont5,rcont6
real(8) :: rcont_par, rcont_perp, isoinvar, anisoinvar2
real(8), allocatable :: ex(:)
complex(8) :: trans1,trans2,trans3,trans4,trans5,trans6
!
!
allocate( ex(input%properties%raman%nstate) )
rcont1 = 0.d0; rcont2 = 0.d0; rcont3 = 0.d0; rcont4 = 0.d0; rcont5 = 0.d0; rcont6 = 0.d0
!
!
dnc = dble(ncell)
temp = temper
if (temp .lt. 1.0e-9) then
   temp = 1.d-9
endif
!
gam1 = gamma1**2/4.d0                         ! broadening parameters from input
gam2 = gamma2**2/4.d0
gam3 = gamma3**2/4.d0
gam4 = gamma4**2/4.d0
zsum = 0.d0
beta = 1.d0/fkwn/temp
sqn1 = dsqrt(dnc)                             ! powers of sqrt( Nc )
sqn2 = dnc
sqn3 = sqn1*sqn2
sqn4 = sqn2*sqn2
sqn5 = sqn2*sqn3
sqn6 = sqn3*sqn3
do i = 1, input%properties%raman%nstate
   arg = -beta*eigen(i)*fhawn
   ex(i) = dexp(arg)
   zsum = zsum + ex(i)                        ! zsum = Sum( exp(-E/kT) )
enddo
if (ex(1) .eq. 0.d0) then                     ! handle numerical problems with low temperatures
   ex(1) = 1.d0
   zsum = 1.d0
endif
!
if (input%properties%raman%molecule) then
   sc = 1.d32*fau3cm3**2
   isoinvar = 1.d0/3.d0 * (deq(1,1) + deq(2,2) + deq(3,3))
   anisoinvar2 = 0.5 * ((deq(1,1)-deq(2,2))**2 + (deq(2,2)-deq(3,3))**2 + (deq(3,3)-deq(1,1))**2) + &
     & 3.d0 * (deq(1,2)**2 + deq(1,3)**2 + deq(2,3)**2)
!
! drop V, divide by epsilon_0 = 1/4pi (au)
   zzsum = zsum/(2.d0*pi)**4
   dde = (eigen(2) - eigen(1))*fhawn
   w_scat = rlas*fhawn - dde              ! Stokes scattering
   w_fact = rlas*fhawn*w_scat**3
   zfact = ex(1)/zzsum                    ! prefactors * exp(-Ej(kT) / Sum( exp(E/kT) )
!
   rcont_par  = (45.d0*abs(isoinvar)**2 + 4.d0*abs(anisoinvar2)) / 45.d0 &
     &                  * (transme1(1)/factorial(1))**2                &
     &                  * zfact
   rcont_perp = 3.d0*abs(anisoinvar2) / 45.d0                            &
     &                  * (transme1(1)/factorial(1))**2                &
     &                  * zfact
   Sab = rcont_par + rcont_perp
   Sab = Sab*sc*w_fact
!
else
   sc = 1.d07
   zzsum = zsum/(omega*fau3cm3*pi*pi)
   dde = (eigen(2) - eigen(1))*fhawn
   w_scat = rlas*fhawn - dde              ! Stokes scattering
   w_fact = rlas*fhawn*w_scat**3
   zfact = ex(j)/zzsum                    ! prefactors * exp(-Ej(kT) / Sum( exp(E/kT) )
!
   trans1 = deq (oct1,oct2)*transme1(1)/factorial(1)/sqn1 + &
        &   d2eq(oct1,oct2)*transme2(1)/factorial(2)/sqn2 + &
        &   d3eq(oct1,oct2)*transme3(1)/factorial(3)/sqn3 + &
        &   d4eq(oct1,oct2)*transme4(1)/factorial(4)/sqn4 + &
        &   d5eq(oct1,oct2)*transme5(1)/factorial(5)/sqn5 + &
        &   d6eq(oct1,oct2)*transme6(1)/factorial(6)/sqn6
   if (input%properties%raman%intphonon) then
      trans2 = d2eq(oct1,oct2)*transme2(1)/factorial(2)/sqn2 + &
        &      d4eq(oct1,oct2)*transme4(1)/factorial(4)/sqn4 + &
        &      d6eq(oct1,oct2)*transme6(1)/factorial(6)/sqn6
      trans3 = d3eq(oct1,oct2)*transme3(1)/factorial(3)/sqn3 + &
        &      d6eq(oct1,oct2)*transme6(1)/factorial(6)/sqn6
      trans4 = d4eq(oct1,oct2)*transme4(1)/factorial(4)/sqn4
      trans5 = d5eq(oct1,oct2)*transme5(1)/factorial(5)/sqn5 
      trans6 = d6eq(oct1,oct2)*transme6(1)/factorial(6)/sqn6
   endif
!
   rcont1    = dnc*zfact*              abs(trans1)**2
   if (input%properties%raman%intphonon) then
      rcont2 = dnc*zfact*(dnc    - 1.)*abs(trans2)**2
      rcont3 = dnc*zfact*(dnc**2 - 1.)*abs(trans3)**2
      rcont4 = dnc*zfact*(dnc**3 - 1.)*abs(trans4)**2
      rcont5 = dnc*zfact*(dnc**4 - 1.)*abs(trans5)**2
      rcont6 = dnc*zfact*(dnc**5 - 1.)*abs(trans6)**2
   endif
!
!
   Sab = rcont1 + rcont2 + rcont3 + rcont4 + rcont5 + rcont6
   Sab = Sab*sc*w_fact
endif
!return scattering energy in Ha
ws = w_scat*fwnha
!
!
deallocate(ex)
!
return
!
end
! 
