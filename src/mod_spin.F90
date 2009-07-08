
#include "maxdefinitions.inc"
module mod_spin
!--------------------------------!
!     spin related variables     !
!--------------------------------!
! spinpol is .true. for spin-polarised calculations
logical::spinpol
! spinorb is .true. for spin-orbit coupling
logical::spinorb
! fixspin type: 0 = none, 1 = global, 2 = local, 3 = global + local
integer::fixspin
! dimension of magnetisation and magnetic vector fields (1 or 3)
integer::ndmag
! ncmag is .true. if the magnetisation is non-collinear, i.e. when ndmag = 3
logical::ncmag
! fixed total spin magnetic moment
real(8)::momfix(3)
! fixed spin moment global effective field in Cartesian coordinates
real(8)::bfsmc(3)
! muffin-tin fixed spin moments
real(8)::mommtfix(3, _MAXATOMS_, _MAXSPECIES_)
! muffin-tin fixed spin moment effective fields in Cartesian coordinates
real(8)::bfsmcmt(3, _MAXATOMS_, _MAXSPECIES_)
! fixed spin moment field mixing parameter
real(8)::taufsm
! second-variational spinor dimension (1 or 2)
integer::nspinor
! external magnetic field in each muffin-tin in lattice coordinates
real(8)::bflmt(3, _MAXATOMS_, _MAXSPECIES_)
! external magnetic field in each muffin-tin in Cartesian coordinates
real(8)::bfcmt(3, _MAXATOMS_, _MAXSPECIES_)
! global external magnetic field in lattice coordinates
real(8)::bfieldl(3)
! global external magnetic field in Cartesian coordinates
real(8)::bfieldc(3)
! external magnetic fields are multiplied by reducebf after each iteration
real(8)::reducebf
! spinsprl if .true. if a spin-spiral is to be calculated
logical::spinsprl
! number of spin-dependent first-variational functions per state
integer::nspnfv
! spin-spiral q-vector in lattice coordinates
real(8)::vqlss(3)
! spin-spiral q-vector in Cartesian coordinates
real(8)::vqcss(3)
end module
