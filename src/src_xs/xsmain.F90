
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOI
! !TITLE: The Developers' Guide for the excited-states implementation within
!   the 
!   EXCITING Code \\ Version
! !AUTHORS: S. Sagmeister and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
!   Welcome to the {\sf XS/EXCITING} code developers' manual
!   This is supposed to collect the routines and modules belonging
!   only to the excited states (TDDFT and BSE) part into one document.
!   \\\\
!   S. Sagmeister\\
!   Graz, 2006
!
!EOI

subroutine xsmain
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_getunit
  implicit none
  character(*), parameter :: thisnam = 'xsmain'
!!$  logical :: tskip

  ! save task
  tasktd = task
  ! remember how often this routine is called
  calledxs = calledxs + 1
  ! 
  if (calledxs == 1) call argparse()

  ! basic initialization
  call xsinit

  if (calledxs == 1) then
     ! check verify constraints
     call tdcheck
     write(unitout,'(a)') 'Info('//thisnam//'): initialization done.'
     ! write info
     call writeinfotd
  end if

!!$  if (tresume) then
!!$     tskip=.true.
!!$  else
!!$     tskip=.false.
!!$  end if
!!$  if (tresume.and.(tasktd.eq.resumetask)) tskip=.false.
!!$  if (tskip) goto 10

  ! task selection
  select case(tasktd)
  case(300)
     ! say hello
     write(*,*)
     write(*,'(a)') 'Copyright (C) 2004-2007 by S. Sagmeister and &
          &C. Ambrosch-Draxl'
     write(*,*)
  case(23)
     ! estimate bandgap from regular grid
     call writebandgapgrid
  case(301)
     ! calculate eigenvectors for q-point set
     call tdgeneigvec
  case(310)
     ! calculate weights for tetrahedron method
     call tetcalccw
  case(320)
     ! parallel version of momentum matrix elements
     call writepmattd(.false.)
  case(321)
     ! ASCII output of momentum matrix elements
     call writepmat_ascii
  case(322)
     ! convert momentum matrix elements file to old format
     call pmattd2orig
  case(330)
     ! calculate matrix elements of exponential expression
     call writeemat
  case(335)
     ! calculate matrix elements of the plane wave (new version)
     call writepwmat
  case(331)
     ! ASCII output of matrix elements of exponential expression
     call writeemat_ascii
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
     ! inversion of dielectric function wrt. xc-kernel
     call idf
  case(351)
     ! linear optics of old version
     call linoptold
  case(396)
     ! convolute dielectric function from tetrahedron method with Lorentzian
     call epsconv
  case(397)
     ! estimate disk-space, cpu-time and memory
     call tdestimate
  case(398)
     ! check ALDA kernel
     call fxc_alda_check
  case(399)
     call testmain
  case(400)
     ! RPA screening
     call screen
  case(498,499)
     ! * debug task
     call init0
     call init1
     call writekpts
  case(497)
     ! * debug task
     call init0
     call init1
     call init2xs
     call writekpts
  case(496)
     ! * debug task
     call init0
  case default
     write(*,*) 'Error('//thisnam//'): task not defined:', tasktd
  end select

!!$10 continue

  ! epilog
  call tdepilog

end subroutine xsmain
