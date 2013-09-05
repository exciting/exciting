! ------------------------------------------------------------------------------
subroutine SPECTRUM(temper, rlas, oct1, oct2, comp_ch, filext)
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
character(2), intent(in) :: comp_ch
character(*), intent(in) :: filext
! local variables
integer :: i,j,ind,ii,iw,niw
real(8) :: arg,beta,dde,factl,gam1,gam2,gam3,gam4,dnc
real(8) :: sc,sqn1,sqn2,sqn3
real(8) :: sqn4,sqn5,sqn6
real(8) :: spec,w,dw,zfact,zsum,zzsum, temp
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
fname = 'RAMAN_SPEC_OC'//comp_ch//filext
open (unit=73, file=fname, status='unknown',form='formatted')
!
! scaling factor for output of scattering contributions
if (input%properties%raman%molecule) then
   ! includes the conversion factor for the polarizability Bohr^3 -> cm^3
   sc = 1.d32*fau3cm3**2
else
   sc = 1.d07
endif
dnc = dble(ncell)
write(66,'(//,46("*"),"   SPECTRUM FOR OC ",a2,"   ",46("*"))') comp(oct1, oct2)
temp = temper
write(66,'(/," Temperature [ K ]: ",17x,f8.2,/)') temp
if (temp .lt. 1.0e-9) then
   temp = 1.d-9
   write(66,'(  " (adjusted to     : ",17x,e8.2,")"/)') temp
endif
write(66,'(/," Partition sums for eigenstates:",/)')  ! header for part sums in output
write(66,'("  state ","         eigenvalue ","   Boltzmann factor ", &
     &           "      partition sum ","         population  "/)')
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
!!     sqnn = dsqrt(dble(nc - 1)/dble(nc))
!sqnn = dsqrt(dble(ncell - 1))/sqn1
!sqnn2 = sqnn/ncell
!sqnn22 = sqnn2/ncell
!!     sqnn3 = dsqrt(dble(nc*nc - 1)/dble(nc*nc))
!sqnn3 = dsqrt(sqn4 - 1.d0)/sqn2
!!     sqnn33 = dsqrt(dble(nc*nc - 1)/dble(nc**5))
!sqnn33 = dsqrt(sqn4 - 1.d0)/sqn5
!!     sqnn4 = dsqrt(dble(nc**3 - nc)/dble(nc**3))
!csnpw(1) = cmplx(sqrt(dnc),0.)                 ! powers of sqrt(N)
!csnpw(2) = csnpw(1)*csnpw(1)
!csnpw(3) = csnpw(2)*csnpw(1)
!csnpw(4) = csnpw(2)*csnpw(2)
!csnpw(5) = csnpw(3)*csnpw(2)
!csnpw(6) = csnpw(3)*csnpw(3)
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
write(73,'("# Laser energy:  ",f8.6," Ha = ",f8.1," cm^-1 = ",f5.2," eV = ",  &
  &        f6.1, " nm",/,"# Temperature: ",f8.1," K",/,"# Damping: ",4(f8.1," cm^-1   |"),/, &
  &       "# Number of eigenvalues: ",i4,/)') rlas,rlas*fhawn,rlas*fhaev,frnmha/rlas,temp, &
  &              gamma1,gamma2,gamma3,gamma4,input%properties%raman%nstate
!
ind = 0
if (input%properties%raman%molecule) then
! drop V, divide by epsilon_0 = 1/4pi (au)
   zzsum = 4.d0*zsum/(pi*(rlas*fhawn)**4)
   write(66,'(/," Scattering contributions in [ 10^-32 cm^2 sr^-1 ] from ", &
     & "process of order",//,4x," Transition energy [ cm^-1 ]",            &
     & 6x,"1",14x,"2",14x,"3",14x,"4",14x,"5",14x,"6",/)')
else
   zzsum = zsum/(omega*fau3cm3*pi*pi*(rlas*fhawn)**4)
   write(66,'(/," Scattering contributions in [ 10^-7 sr^-1 cm^-1 ] from ", &
     & "process of order",//,4x," Transition energy [ cm^-1 ]",            &
     & 6x,"1",14x,"2",14x,"3",14x,"4",14x,"5",14x,"6",/)')
endif
do i = 1,input%properties%raman%nstate
   do j = i-1,1,-1
      zfact = ex(j)/zzsum                    ! exp(-Ej(kT) / Sum( exp(E/kT) )
      ind = ind + 1
!
      trans1 = deq (oct1,oct2)*transme1(ind)/factorial(1)/sqn1 + &
     &         d2eq(oct1,oct2)*transme2(ind)/factorial(2)/sqn2 + &
     &         d3eq(oct1,oct2)*transme3(ind)/factorial(3)/sqn3 + &
     &         d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4 + &
     &         d5eq(oct1,oct2)*transme5(ind)/factorial(5)/sqn5 + &
     &         d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
      if (input%properties%raman%intphonon) then
         trans2 = d2eq(oct1,oct2)*transme2(ind)/factorial(2)/sqn2 + &
     &            d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4 + &
     &            d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
         trans3 = d3eq(oct1,oct2)*transme3(ind)/factorial(3)/sqn3 + &
     &            d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
         trans4 = d4eq(oct1,oct2)*transme4(ind)/factorial(4)/sqn4
         trans5 = d5eq(oct1,oct2)*transme5(ind)/factorial(5)/sqn5 
         trans6 = d6eq(oct1,oct2)*transme6(ind)/factorial(6)/sqn6
      endif
!
      rcont1(ind) = dnc*zfact*              abs(trans1)**2
      if (input%properties%raman%intphonon) then
         rcont2(ind) = dnc*zfact*(dnc    - 1.)*abs(trans2)**2
         rcont3(ind) = dnc*zfact*(dnc**2 - 1.)*abs(trans3)**2
         rcont4(ind) = dnc*zfact*(dnc**3 - 1.)*abs(trans4)**2
         rcont5(ind) = dnc*zfact*(dnc**4 - 1.)*abs(trans5)**2
         rcont6(ind) = dnc*zfact*(dnc**5 - 1.)*abs(trans6)**2
      endif
!
!
      dde = (eigen(i) - eigen(j))*fhawn
      if (input%properties%raman%intphonon) then
         write(66,'(1x,I3," ->",I3,3x,f8.2,5x,6(f14.8,"  |"))') j,i,dde,             &
     &           rcont1(ind)*sc,rcont2(ind)*sc, &
     &           rcont3(ind)*sc,rcont4(ind)*sc, &
     &           rcont5(ind)*sc,rcont6(ind)*sc
      else
         write(66,'(3x,I3," ->",I3,5x,f8.2,7x,f14.8)') j,i,dde,             &
     &           rcont1(ind)*sc
      endif
   enddo
enddo
!
! write out spectra
niw = 0
!factl = 2.d0*pi*(rlas*fhawn)**3
factl = 2.d0*pi*(rlas*fhawn)**4
do iw = 1, input%properties%raman%energywindow%points+1
   ind = 0
   spec = 0.0d0
   w = input%properties%raman%energywindow%intv(1) + dble(iw - 1)*dw
   do i = 1,input%properties%raman%nstate
      do j = i-1,1,-1
         ind = ind + 1
         dde = (eigen(i) - eigen(j))*fhawn
         spec = spec + &                        ! line broadening by Lorentzian
     &          rcont1(ind)*gamma1/((dde-w)**2 + gam1) + &
     &          rcont2(ind)*gamma2/((dde-w)**2 + gam2) + &
     &          rcont3(ind)*gamma3/((dde-w)**2 + gam3) + &
     &          rcont4(ind)*gamma4/((dde-w)**2 + gam4) + &
     &          rcont5(ind)*gamma4/((dde-w)**2 + gam4) + &   ! 5th and 6th use gamma4
     &          rcont6(ind)*gamma4/((dde-w)**2 + gam4)
      enddo
   enddo
   spec = spec*(rlas*fhawn - w)**3/factl
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
