! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

!> Exciting-specific constants
module constants
  use precision, only: dp, i32
  use iso_c_binding

  implicit none
  private  

  !> Pi 
  real(dp), public, parameter :: pi = 3.1415926535897932385_dp
  !> 2 * Pi
  real(dp), public, parameter :: twopi = 6.2831853071795864769_dp
  !> 4 * Pi 
  real(dp), public, parameter :: fourpi = 12.566370614359172954_dp

  !> Square root of two. 
  real(dp), public, parameter :: sqrt_two = 1.4142135623730950488_dp

  ! Complex initialisation constants
  !> Complex zero
  complex(dp), public, parameter :: zzero = (0._dp, 0._dp)
  !> Complex half 
  complex(dp), public, parameter :: zhalf = (0.5_dp, 0._dp)
  !> Complex one
  complex(dp), public, parameter :: zone = (1._dp, 0._dp)
  !> Complex i 
  complex(dp), public, parameter :: zi = (0._dp, 1._dp)

  ! Real initialisation constants
  !> Real zero
  real(dp), public, parameter :: real_zero = 0._dp
  !> Real one
  real(dp), public, parameter :: real_one = 1._dp

    
  ! TODO(Alex). Issue #24. Replace maxatoms and maxspecies with values from input

  !> Maximum number of different species handled by Exciting. 
  !> This is likely an upper bound originally defined in a legacy version of the code
  !> which was and is used for declaring static arrays
  integer(i32), public, parameter :: maxspecies = 8
  !> Maximum number of atoms per species. 
  !> Upper bound for static array declaration. See mod_atoms.F90, for example 
  integer(i32), public, parameter :: maxatoms = 250
  !> Maximum number of linear augmented plane waves per species
  !> Easy recursive algorithms for generating spherical harmonics are stable up to ~ 50
  !> Upper bound for static array declaration. See mod_APW_LO.F90, for example 
  integer(i32), public, parameter :: maxlapw = 50 

  !> Kronecker delta matrix 
  integer(i32), public, parameter :: krondelta (3, 3) = &
      reshape ( (/ 1, 0, 0, & 
                   0, 1, 0, &
                   0, 0, 1  /), (/ 3, 3 /))

  !> Pauli spin matrix, $\sigma_x$
  !> \f[ sigma_x = ( 0  1 )  
  !>               ( 1  0 )  \f]         
  complex(dp), public, parameter, dimension(2, 2) :: sigma_x = reshape(&
     [(0._dp, 0._dp), (1._dp, 0._dp), &
      (1._dp, 0._dp), (0._dp ,0._dp) ], [2, 2])

  !> Pauli spin matrix, $\sigma_y$
  !> \f[ sigma_y = ( 0 -i )  
  !>               ( i  0 )  \f]
  !> Note, fortran stores data columnwise hence the declaration.
  !> Use of transpose would avoid this but breaks compatibility with Intel 2015 
  complex(dp), public, parameter, dimension(2, 2) :: sigma_y = reshape(&
     [(0._dp,  0._dp), (0._dp, 1._dp),   &
      (0._dp, -1._dp), (0._dp, 0._dp) ], [2, 2])

  !> Pauli spin matrix, $\sigma_z$
  !> \f[ sigma_z = ( 1  0 )  
  !>               ( 0 -1 )  \f]    
  complex(dp), public, parameter, dimension(2, 2) :: sigma_z = reshape(&
    [(1._dp, 0._dp), ( 0._dp, 0._dp), &
     (0._dp, 0._dp), (-1._dp, 0._dp)  ], [2, 2])
    
  !> Pauli spin matrices
  ! Note, no nice constructor for 3D arrays 
  complex(dp), public :: sigmat (2, 2, 3)
  data sigmat / (0._dp, 0._dp), (1._dp, 0._dp), (1._dp, 0._dp), (0._dp, &
  & 0._dp), (0._dp, 0._dp), (0._dp, 1._dp), (0._dp,-1._dp), (0._dp, 0._dp), &
  & (1._dp, 0._dp), (0._dp, 0._dp), (0._dp, 0._dp), (-1._dp, 0._dp) /   

  ! TODO(Alex) This should be moved. Should be with other spherical harmonics  
  !> spherical harmonic for l=m=0. 
  real(dp), public, parameter :: y00 = 0.28209479177387814347_dp  

  ! TODO(Alex) This should be moved. Not initialised with values and not parameter, hence not a constant 
  !> array of i**l values
  complex(dp), public, allocatable :: zil (:)

  !> Set of lower case alphabetic strings.
  character(*), public, parameter :: lower_case_alphabet_set = 'abcdefghijklmnopqrstuvwxyz'
  !> Set of upper case alphabetic strings.
  character(*), public, parameter :: upper_case_alphabet_set = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  !> Set of digid strings.
  character(*), public, parameter :: digit_set = '0123456789'

  !> Value for unintialized c integers.
  integer(c_int), public, parameter :: uninit_c_int = 9999999  

end module  
