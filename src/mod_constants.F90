#include "maxdefinitions.inc"
module mod_constants
!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! square root of two
real(8), parameter :: sqtwo=1.4142135623730950488d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
! complex constants
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zhalf=(0.5d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0)
! array of i**l values
complex(8), allocatable :: zil(:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
complex(8) sigmat(2,2,3)
data sigmat / (0.d0,0.d0), (1.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
              (0.d0,0.d0), (0.d0,1.d0),(0.d0,-1.d0), (0.d0,0.d0), &
              (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0),(-1.d0,0.d0) /
! Boltzmann constant in Hartree/kelvin (CODATA 2006)
real(8), parameter :: kboltz=3.166815343d-6
end module
