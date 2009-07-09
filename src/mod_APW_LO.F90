

#include "maxdefinitions.inc"
module mod_APW_LO
!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! maximum allowable APW order
integer, parameter :: maxapword=4
! APW order
integer::apword(0:_MAXLAPW_, _MAXSPECIES_)
! maximum of apword over all angular momenta and species
integer::apwordmax
! APW initial linearisation energies
real(8)::apwe0(maxapword, 0:_MAXLAPW_, _MAXSPECIES_)
! APW linearisation energies
real(8), allocatable :: apwe(:, :, :)
! APW derivative order
integer::apwdm(maxapword, 0:_MAXLAPW_, _MAXSPECIES_)
! apwve is .true. if the linearisation energies are allowed to vary
logical::apwve(maxapword, 0:_MAXLAPW_, _MAXSPECIES_)
! APW radial functions
real(8), allocatable :: apwfr(:, :, :, :, :)
! derivate of radial functions at the muffin-tin surface
real(8), allocatable :: apwdfr(:, :, :)
! maximum number of local-orbitals
integer, parameter :: maxlorb=20
! maximum allowable local-orbital order
integer, parameter :: maxlorbord=4
! number of local-orbitals
integer::nlorb(_MAXSPECIES_)
! maximum nlorb over all species
integer::nlomax
! total number of local-orbitals
integer::nlotot
! local-orbital order
integer::lorbord(maxlorb, _MAXSPECIES_)
! local-orbital angular momentum
integer::lorbl(maxlorb, _MAXSPECIES_)
! maximum lorbl over all species
integer::lolmax
! (lolmax+1)^2
integer::lolmmax
! local-orbital initial energies
real(8)::lorbe0(maxlorbord, maxlorb, _MAXSPECIES_)
! local-orbital energies
real(8), allocatable :: lorbe(:, :, :)
! local-orbital derivative order
integer::lorbdm(maxlorbord, maxlorb, _MAXSPECIES_)
! lorbve is .true. if the linearisation energies are allowed to vary
logical::lorbve(maxlorbord, maxlorb, _MAXSPECIES_)
! local-orbital radial functions
real(8), allocatable :: lofr(:, :, :, :)
! energy step size for locating the band energy
!replaced by inputstructurereal(8)::deband
end module

