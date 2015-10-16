
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

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
! Dipole correction energy      
      Real(8)  :: endipc
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
! DFT-D2 dispersion correction
      Real (8) :: e_disp
! DFT-1/2 contribution to total energy
      Real (8) :: engyhalf
! total energy
      Real (8) :: engytot
! kinetic energies of KS states
      Real (8), allocatable :: engyknst(:,:)
End Module
!
