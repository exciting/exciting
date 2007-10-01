subroutine constants
use modeos
implicit none
pi=2.d0*asin(1.d0)
! CODATA 98 units
planck_si=1.054571596d-34
eps0_si=8.854187817d-12
kb_si=1.3806503d-23
aumass_si=9.10938188d-31
aucharge_si=1.602176462d-19
aulength_si=0.5291772083d-10
! Avogadro's number
avogadro=6.02214199d23
! Hartree in SI
auenergy_si=(aumass_si/(planck_si**2))*(aucharge_si**2/(4.d0*pi*eps0_si))**2
! atomic time unit in SI
autime_si=planck_si/auenergy_si
! Boltzmann constant in atomic units
kb_au=kb_si/auenergy_si
! atomic mass unit in electron masses
mu_au=1822.888480
! atomic pressure unit in GPa
aupress_gpa=1.d-9*aumass_si/(aulength_si*autime_si**2)
return
end subroutine
