! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOI
! !TITLE: The XS/EXCITING Code (eXited States) Manual \\ Version 0.9
! !AUTHORS: S. Sagmeister and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
!   Welcome to the {\sf XS/EXCITING} code developers' manual.
!   This manual is supposed to collect the routines and modules belonging
!   exclusively to the excited states (TDDFT and BSE) part into one document.
!   The content of this manual is partially taken from the author's PhD thesis.
!   \\\&
!   S. Sagmeister\&
!   Leoben, August 2008
!
!EOI

subroutine xsmain(plan, nxstasks)
  use modinput
  use modmpi
  use mod_misc, only: task
  use modxs, only: fhdf5
  use mod_exciton_wf
  use mod_hdf5

  implicit none

  type(plan_type), intent(in) :: plan
  integer(4), intent(in) :: nxstasks

  integer(4) :: i

  ! initialization of hdf5 output
#ifdef _HDF5_
  if(mpiglobal%rank == 0) then
    call hdf5_initialize()
    fhdf5="bse_output.h5"
    call hdf5_create_file(fhdf5)
  end if
#endif
  
  ! task selection, loop over first nxstasks specified in passed plan
  do i = 1, nxstasks

    ! Set task number
    task = plan%doonlyarray(i)%doonly%tasknumber

    ! initialization for xs tasks (dependent on task number)
    call xsinit(i,plan)

    select case(task)

      ! Taskname 'writebandgapgrid'
      case(23)
        ! estimate bandgap from regular grid
        call writebandgapgrid

      ! Taskname 'xsgeneigvec'
      case(301)
        ! generate eigenvectors, eigenvalues, occupancies and MT-coefficients
        ! for q-point set
        call xsgeneigveclauncher

      ! Taksname 'writepmatxs'
      case(320)
        ! parallel version of momentum matrix elements
        call writepmatxs

      ! Taksname 'writepmatasc'
      case(321)
        ! ASCII output of momentum matrix elements
        call writepmatasc

      ! Taskname 'pmatxs2orig'
      case(322)
        ! convert momentum matrix elements file to old format
        call pmatxs2orig

      ! Taskname 'writeemat'
      case(330)
        ! calculate matrix elements of exponential expression (band combs)
        call writeemat

      ! Taskname 'writeematasc'
      case(331)
        ! ASCII output of matrix elements of exponential expression
        call writeematasc

      ! Taskname 'writepwmat'
      case(335)
        ! calculate matrix elements of the plane wave (simple version for checking)
        call writepwmat

      ! Taskname 'emattest'
      case(339)
        ! check relation between matr. el. of exp. and mom. matr. el.
        call emattest

      ! Taskname 'df'
      case(340)
        ! Kohn Sham response function
        call df

      ! Taskname 'x0toasc'
      case(341)
        ! ASCII output of Kohn Sham response function
        call x0toasc

      ! Taskname 'x0tobin'
      case(342)
        ! binary output of Kohn Sham response function
        call x0tobin

      ! Taskname 'idf'
      case(350)
        ! inverse of dielectric function - solve Dyson equation for xc-kernel
        call idf

      ! Taskname 'fxc_alda_check'
      case(398)
        ! check ALDA kernel
        call fxc_alda_check

      ! Taksname 'scrgeneigvec'
      case(401)
        ! generate eigenvectors, eigenvalues, occupancies and APW MT coefficients
        ! for screening and BSE(-kernel)
        call xsgeneigveclauncher

      ! Taskname 'scrwritepmat'
      case(420)
        ! momentum matrix elements for screening
        call scrwritepmat

      ! Taskname 'screen'
      case(430)
        ! RPA screening
        call screenlauncher

      ! Taskname 'scrcoulint'
      case(440)
        ! screened Coulomb interaction
        call scrcoulintlauncher

      ! Taskname 'exccoulint'
      case(441)
        ! exchange Coulomb interaction
        call exccoulintlauncher

      ! Taskname 'bse'
      case(445)
        ! Bethe-Salpeter equation
        call bselauncher

      ! Taskname 'bsegenspec'
      case(446)
        ! regenerate BSE spectrum from exciton output
        call bsegenspec

      ! Taskname 'writebevec'
      case(447)
        ! ASCII output of BSE eigenvectors
        call writeexcevec

      ! Taskname 'writekpathweights'
      case(448)
        ! ASCII output of excitonic weights
        call writekpathweights

      ! Taskname 'bsesurvey'
      case(449)
        ! BSE transition survey
        call bsesurvey

      ! Taskname 'kernxc_bse'
      case(450)
        ! BSE-kernel
        call kernxc_bse

      ! Taskname 'xsestimate'
      case(700)
        ! estimate disk-space, cpu-time and memory
        call xsestimate

      ! Taskname 'xstiming'
      case(701)
        ! test timing
        call xstiming

      ! Taskname 'excitonWavefunction'
      case(710)
        ! Polt of TDA exciton wave function
        call plot_excitonWavefunction

      ! Taskname 'testmain'
      case(999)
        call testmain

      case default
        write(*,*)
        write(*,*) 'Error(xsmain): task not defined:', task
        write(*,*)
        call terminate

     end select

     ! Finalization for each xs task (task dependent)
     call xsfinit

   end do

  ! Finalization of hdf5 output
#ifdef _HDF5_
  if (mpiglobal%rank == 0) then
    call hdf5_finalize()
  end if
#endif

end subroutine xsmain
