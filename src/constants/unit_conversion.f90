! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020 

! For all new conversions, please:
! a) Define using CODATA 2018 (see reference below)
! b) Define using double precision (dp): ha_to_ev = 27.21138386_dp 
! c) Use a naming convention like: ha_to_ev, not hartree_to_ev 
! d) Don't state the value in the comment as comments cannot be tested but are subject to change 
!    - values in comments left from refactor, and require addressing 

! "The 2018 CODATA Recommended Values of the Fundamental Physical Constants"
! (Web Version 8.0). Database developed by J. Baker, M. Douma, and S. Kotochigova.
! Available at http://physics.nist.gov/constants,   

!> Unit conversion 
module unit_conversion
    use precision, only: dp 
    implicit none
    private

    ! TODO(Alex) Issue #20. Update to CODATA 2018 physical units if required
    ! (and NOTE) the last decimal place has been rounded erroneously
    !> Conversion from Hartrees to electron volts (CODATA 2006)
    ! 1 Hartree = 27.211 383 86(68) eV
    real(8), public, parameter :: hartree_to_ev = 27.21138386d0

    ! TODO(Alex) Issue #20. Update to CODATA 2018 physical units if required
    !> Conversion from Hartrees to cm^{-1} (CODATA 2006):
    ! 1 Hartree / (hc) = 2.194 746 313 705(15) * 10^7 m^{-1}
    real(dp), public, parameter :: hartree_to_inv_cm = 2.194746313705e-5_dp

    ! TODO(Alex) Issue #20. Update to CODATA 2018 physical units if required
    !> Conversion from Hartrees to THz (CODATA 2006):
    ! 1 Hartree / h = 6.579 683 920 722(44) * 10^{15} Hz
    real(dp), public, parameter :: hartree_to_thz = 6.579683920722e-3_dp

end module unit_conversion