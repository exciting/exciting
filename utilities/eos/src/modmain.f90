
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmain
! crystal name
character(256) cname
! number of atoms
integer natoms
! EOS type
integer etype
! number of volume points to plot
integer nvplt
! volume plot range
real(8) vplt1,vplt2
! number of energy data points to fit
integer nevpt
! volume and energy data point sets
real(8), allocatable :: vpt(:)
real(8), allocatable :: ept(:)
! maximum number of parameters for an EOS
integer, parameter :: maxparam=100
! number of parameters
integer nparam
! EOS name
character(256) ename(2)
! optimized parameter set
real(8) popt(maxparam)
! parameter names
character(256) pname(maxparam)


!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
! CODATA 2006 constants
! Bohr in SI units
real(8), parameter :: bohr_si=0.52917720859d-10
! electron mass in SI units
real(8), parameter :: emass_si=9.10938215d-31
! atomic unit of time in SI units
real(8), parameter :: autime_si=2.418884326505d-17
! atomic pressure unit in GPa
real(8), parameter :: aupress_gpa=1.d-9*emass_si/(bohr_si*autime_si**2)

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer version(3)
data version /1,4,0/

end module
