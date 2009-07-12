
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOI
! !TITLE: The XS/EXCITING Code (eXited States) Manual \\ Version 0.9
! !AUTHORS: S. Sagmeister and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
!   Welcome to the {\sf XS/EXCITING} code developers' manual.
!   This manual is supposed to collect the routines and modules belonging
!   exclusively to the excited states (TDDFT and BSE) part into one document.
!   The content of this manual is partially taken from the author's PhD thesis.
!   \\\\
!   S. Sagmeister\\
!   Leoben, August 2008
!
!EOI

module modxsmain
use mod_misc
implicit none
integer :: nxstasks
character(256) :: xstasks(maxtasks)
end module

subroutine xsmain
  use modmain
  use modmpi
  use modtetra
  use modxs
  implicit none
  ! initialization
  call xsinit
  ! task selection
  select case(task)
  case(23)
     ! estimate bandgap from regular grid
     call writebandgapgrid
  case(700)
     ! estimate disk-space, cpu-time and memory
     call xsestimate
  case(701)
     ! test timing
     call xstiming
  case(301)
     ! generate eigenvectors, eigenvalues, occupancies and MT-coefficients
     ! for q-point set
     call xsgeneigvec
  case(310)
     ! calculate weights for tetrahedron method
     call tetcalccw
  case(320)
     ! parallel version of momentum matrix elements
     call writepmatxs
  case(321)
     ! ASCII output of momentum matrix elements
     call writepmatasc
  case(322)
     ! convert momentum matrix elements file to old format
     call pmatxs2orig
  case(330)
     ! calculate matrix elements of exponential expression (band combs)
     call writeemat
  case(331)
     ! ASCII output of matrix elements of exponential expression
     call writeematasc
  case(335)
     ! calculate matrix elements of the plane wave (simple version for checking)
     call writepwmat
  case(339)
     ! check relation between matr. el. of exp. and mom. matr. el.
     call emattest
  case(340)
     ! Kohn Sham response function
     call df
  case(341)
     ! ASCII output of Kohn Sham response function
     call x0toasc
  case(342)
     ! binary output of Kohn Sham response function
     call x0tobin
  case(350)
     ! inverse of dielectric function - solve Dyson equation for xc-kernel
     call idf
  case(396)
     ! convolute dielectric function from tetrahedron method with Lorentzian
     call epsconv
  case(398)
     ! check ALDA kernel
     call fxc_alda_check
  case(401)
     ! generate eigenvectors, eigenvalues, occupancies and APW MT coefficients
     ! for screening and BSE(-kernel)
     call scrgeneigvec
  case(410)
     ! calculate weights for tetrahedron method (screening)
     call scrtetcalccw
  case(420)
     ! momentum matrix elements for screening
     call scrwritepmat
  case(430)
     ! RPA screening
     call screen
  case(440)
     ! screened Coulomb interaction
     call scrcoulint
  case(441)
     ! exchange Coulomb interaction
     call exccoulint
  case(445)
     ! Bethe-Salpeter equation
     call bse
  case(450)
     ! BSE-kernel
     call kernxc_bse
  case(451)
     ! BSE-kernel, simple version, very slow
     call kernxc_bse3
  case(499)
     ! call to test xs-routine
     call testxs
  case default
     write(*,*)
     write(*,*) 'Error(xsmain): task not defined:', task
     write(*,*)
     call terminate
  end select
  ! summarize information on run
  call xsfinit
end subroutine xsmain
