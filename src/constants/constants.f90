! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

!> Exciting-specific constants
module constants
  use precision, only: dp 
  implicit none
  private  

  !> Pi 
  Real(8), Public, Parameter :: pi = 3.1415926535897932385d0
  !> 2 * Pi
  Real(8), Public, Parameter :: twopi = 6.2831853071795864769d0
  !> 4 * Pi 
  Real(8), Public, Parameter :: fourpi = 12.566370614359172954d0

  !> Square root of two. 
  Real(8), Public, Parameter :: sqrt_two = 1.4142135623730950488d0

  ! Complex initialisation constants
  !> Complex zero
  Complex (8), Public, Parameter :: zzero = (0.d0, 0.d0)
  !> Complex half 
  Complex (8), Public, Parameter :: zhalf = (0.5d0, 0.d0)
  !> Complex one
  Complex (8), Public, Parameter :: zone = (1.d0, 0.d0)
  !> Complex i 
  Complex (8), Public, Parameter :: zi = (0.d0, 1.d0)
    
  ! TODO(Alex). Issue #24. Replace maxatoms and maxspecies with values from input

  !> Maximum number of different species handled by Exciting. 
  !> This is likely an upper bound originally defined in a legacy version of the code
  !> which was and is used for declaring static arrays
  Integer, Public, Parameter :: maxspecies = 8 
  !> Maximum number of atoms per species. 
  !> Upper bound for static array declaration. See mod_atoms.F90, for example 
  Integer, Public, Parameter :: maxatoms = 200 
  !> Maximum number of linear augmented plane waves per species
  !> Easy recursive algorithms for generating spherical harmonics are stable up to ~ 50
  !> Upper bound for static array declaration. See mod_APW_LO.F90, for example 
  Integer, Public, Parameter :: maxlapw = 50 

  !> Kronecker delta matrix 
  Integer, Public, Parameter :: krondelta (3, 3) = &
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
  complex(8), public :: sigmat (2, 2, 3)
  data sigmat / (0.d0, 0.d0), (1.d0, 0.d0), (1.d0, 0.d0), (0.d0, &
  & 0.d0), (0.d0, 0.d0), (0.d0, 1.d0), (0.d0,-1.d0), (0.d0, 0.d0), &
  & (1.d0, 0.d0), (0.d0, 0.d0), (0.d0, 0.d0), (-1.d0, 0.d0) /   

  ! TODO(Alex) This should be moved. Should be with other spherical harmonics  
  !> spherical harmonic for l=m=0. 
  real(8), public, parameter :: y00 = 0.28209479177387814347d0  

  ! TODO(Alex) This should be moved. Not initialised with values and not parameter, hence not a constant 
  !> array of i**l values
  complex(8), public, allocatable :: zil (:)

end module  
