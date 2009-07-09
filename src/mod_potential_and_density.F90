

#include "maxdefinitions.inc"
module mod_potential_and_density
use modinput
!-----------------------------------------!
!     potential and density variables     !
!-----------------------------------------!
! exchange-correlation functional type
!replaced by inputstructureinteger::xctype
! exchange-correlation functional description
character(256)::xcdescr
! exchange-correlation functional spin treatment
integer::xcspin
! exchange-correlation functional density gradient treatment
integer::xcgrad

integer::maxncv
! muffin-tin charge density
real(8), allocatable :: rhomt(:, :, :)
! interstitial real-space charge density
real(8), allocatable :: rhoir(:)
! muffin-tin magnetisation vector field
real(8), allocatable :: magmt(:, :, :, :)
! interstitial magnetisation vector field
real(8), allocatable :: magir(:, :)
! muffin-tin Coulomb potential
real(8), allocatable :: vclmt(:, :, :)
! interstitial real-space Coulomb potential
real(8), allocatable :: vclir(:)
! order of polynomial for pseudocharge density
!replaced by inputstructureinteger::npsden
! muffin-tin exchange-correlation potential
real(8), allocatable :: vxcmt(:, :, :)
! interstitial real-space exchange-correlation potential
real(8), allocatable :: vxcir(:)
! muffin-tin exchange-correlation magnetic field
real(8), allocatable :: bxcmt(:, :, :, :)
! interstitial exchange-correlation magnetic field
real(8), allocatable :: bxcir(:, :)
! nosource is .true. if the field is to be made source-free
!replaced by inputstructurelogical::nosource
! muffin-tin effective potential
real(8), allocatable :: veffmt(:, :, :)
! interstitial effective potential
real(8), allocatable :: veffir(:)
! G-space interstitial effective potential
complex(8), allocatable :: veffig(:)
! muffin-tin exchange energy density
real(8), allocatable :: exmt(:, :, :)
! interstitial real-space exchange energy density
real(8), allocatable :: exir(:)
! muffin-tin correlation energy density
real(8), allocatable :: ecmt(:, :, :)
! interstitial real-space correlation energy density
real(8), allocatable :: ecir(:)
! type of mixing to use for the potential
!replaced by inputstructureinteger::mixtype
! adaptive mixing parameters
!replaced by inputstructurereal(8)::beta0
!replaced by inputstructurereal(8)::betainc
!replaced by inputstructurereal(8)::betadec
end module

