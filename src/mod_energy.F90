!
!
#include "maxdefinitions.inc"
Module mod_energy
!--------------------------!
!     energy variables     !
!--------------------------!
! eigenvalue sum
      Real (8) :: evalsum
! electron kinetic energy
      Real (8) :: engykn
! core electron kinetic energy
      Real (8) :: engykncr
! nuclear-nuclear energy
      Real (8) :: engynn
! electron-nuclear energy
      Real (8) :: engyen
! Hartree energy
      Real (8) :: engyhar
! Coulomb energy (E_nn + E_en + E_H)
      Real (8) :: engycl
! electronic Coulomb potential energy
      Real (8) :: engyvcl
! Madelung term
      Real (8) :: engymad
! exchange-correlation potential energy
      Real (8) :: engyvxc
! exchange-correlation effective field energy
      Real (8) :: engybxc
! energy of external global magnetic field
      Real (8) :: engybext
! energy of muffin-tin magnetic fields (non-physical)
      Real (8) :: engybmt
! exchange energy
      Real (8) :: engyx
! correlation energy
      Real (8) :: engyc
! compensating background charge energy
      Real (8) :: engycbc
! total energy
      Real (8) :: engytot
End Module
!
