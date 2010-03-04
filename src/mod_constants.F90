
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_constants
!-----------------------------!
!     numerical constants     !
!-----------------------------!
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: twopi = 6.2831853071795864769d0
      Real (8), Parameter :: fourpi = 12.566370614359172954d0
! square root of two
      Real (8), Parameter :: sqtwo = 1.4142135623730950488d0
! spherical harmonic for l=m=0
      Real (8), Parameter :: y00 = 0.28209479177387814347d0
! complex constants
      Complex (8), Parameter :: zzero = (0.d0, 0.d0)
      Complex (8), Parameter :: zhalf = (0.5d0, 0.d0)
      Complex (8), Parameter :: zone = (1.d0, 0.d0)
      Complex (8), Parameter :: zi = (0.d0, 1.d0)
! array of i**l values
      Complex (8), Allocatable :: zil (:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
      Complex (8) sigmat (2, 2, 3)
      Data sigmat / (0.d0, 0.d0), (1.d0, 0.d0), (1.d0, 0.d0), (0.d0, &
     & 0.d0), (0.d0, 0.d0), (0.d0, 1.d0), (0.d0,-1.d0), (0.d0, 0.d0), &
     & (1.d0, 0.d0), (0.d0, 0.d0), (0.d0, 0.d0), (-1.d0, 0.d0) /
! Boltzmann constant in Hartree/kelvin (CODATA 2006)
      Real (8), Parameter :: kboltz = 3.166815343d-6
  ! Kronecker delta
      Integer, Parameter :: krondelta (3, 3) = reshape ( (/ 1, 0, 0, 0, &
     & 1, 0, 0, 0, 1 /), (/ 3, 3 /))
! conversion from Hartrees to electron volts (CODATA 2006):
! 1 Hartree = 27.211 383 86(68) eV
      real(8), parameter :: h2ev = 27.21138386d0
! conversion from Hartrees to cm^{-1} (CODATA 2006):
! 1 Hartree / (hc) = 2.194 746 313 705(15) * 10^7 m^{-1}
      real(8), parameter :: h2cm1 = 2.194746313705d5
! conversion from Hartrees to THz (CODATA 2006):
! 1 Hartree / h = 6.579 683 920 722(44) * 10^{15} Hz
      real(8), parameter :: h2thz = 6.579683920722d3
End Module
!
