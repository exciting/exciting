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
! TODO (Max) Issue 64: Use task names instead of numbers 
!                     (see example task\_screened\_Coulomb)
! TODO (Max) Issue 65 Collect launchers in module(s)  
!              (see example ph\_screening\_launcher)
!EOI

subroutine xsmain(plan, nxstasks)
  use modinput
  use modmpi
  use mod_misc, only: task
  use mod_exciton_wf
  use mod_hdf5, only: fhdf5
  use m_write_hdf5, only: fhdf5_inter
  use xhdf5, only: xhdf5_type

  use mod_write_screen, only: write_screen

  use phonon_screening, only: phonon_screening_launcher
  use expand_add_eps, only: expand_add_eps_launcher
  use write_screening, only: write_screening_launcher
  use xhdf5, only: xhdf5_type
  use xstring, only: validate_filename
  use fastBSE, only: fastBSE_main, fastBSE_sanity_checks
  use fastBSE_write_wfplot, only: fastBSE_write_u
  use fastBSE_transitions, only: fastBSE_setup_transitions
  use fastBSE_isdf, only: fastBSE_isdf_cvt
  use fastBSE_isdf_tests, only: fastBSE_isdf_vexc_test
  use modxs, only: unitout
  use write_screening, only: write_screening_launcher
  
  implicit none

  !> Screening from polar phonons
  integer, parameter :: task_phonon_screening = 431
  !> Expanding dielectric matrix from unit cell to super cell
  integer, parameter :: task_expand_add_eps = 432
  !> Screened Coulomb interaction for BSE
  integer, parameter :: task_screened_coulomb = 440 
  
  !> Writing dielectric matrix for all non-reduced q-vectors
  integer, parameter :: task_write_dielectric_matrix = 442
  !> Writing screened Coulomb matrix for all non-reduced q-vectors
  integer, parameter :: task_write_screened_coulomb = 443

  type(plan_type), intent(in) :: plan
  integer(4), intent(in) :: nxstasks
  integer(4) :: i
  character(:), allocatable :: ghdf5
  type(xhdf5_type) :: h5

  ! initialization of hdf5 output
  fhdf5 = trim( adjustl( input%xs%h5fname ))
  ghdf5 = trim( adjustl( input%xs%h5gname ))

  call terminate_if_false(validate_filename(fhdf5, '.h5'), 'HDF5 file name for bse output is not valid.')

  call h5%initialize(fhdf5, mpiglobal%comm)
  if (ghdf5 /= '/') call h5%initialize_group('/', ghdf5)
  call h5%finalize()

  fhdf5_inter = 'bse_matrix.h5'
  call h5%initialize(fhdf5_inter, mpiglobal%comm)
  call h5%finalize()

  do i = 1, nxstasks
     task = plan%doonlyarray(i)%doonly%tasknumber
  end do

  
  ! task selection, loop over first nxstasks specified in passed plan
  do i = 1, nxstasks

    ! Set task number
    task = plan%doonlyarray(i)%doonly%tasknumber

    ! initialization for xs tasks (dependent on task number)case(321)
        !   
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
      
      ! Taskname 'phonon_screening'
      case(task_phonon_screening) 
        call phonon_screening_launcher

      ! Taskname 'expand_eps'
      case(task_expand_add_eps)
        ! Expanding dielectric matrix
        call expand_add_eps_launcher
 
      ! Taskname 'scrcoulint'
      case(440)
        ! screened Coulomb interaction
        call scrcoulintlauncher

      ! Taskname 'exccoulint'
      case(441)
        ! exchange Coulomb interaction
        call exccoulintlauncher

      ! Taskname 'write_dielectric_matrix'
      case (task_write_dielectric_matrix)
          ! Expanding dielectric matrix
          call write_screening_launcher('write_dielectric_matrix', fhdf5, ghdf5, mpiglobal)

          ! Taskname 'write_screened_coulomb'
      case (task_write_screened_coulomb)
          ! Expanding dielectric matrix
          call write_screening_launcher('write_screened_coulomb', fhdf5, ghdf5, mpiglobal)

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

      ! Taskname 'write_wfplot'
      case(451)
        ! write real-space XS wfcts to file
        call fastBSE_sanity_checks(mpiglobal, input)
        call fastBSE_write_u(fhdf5, ghdf5, input, mpiglobal)

      ! Taskname 'write_screen'
      case(452)
        ! write screened Coulomb potential to file
        call write_screen

      ! Taskname 'fastBSE_main'
      case(501)
        call fastBSE_sanity_checks(mpiglobal, input)
        call fastBSE_main(mpiglobal, input, fhdf5, ghdf5, unitout)

      case(510)
      ! Taskname 'fastBSE_setup_transitions'
        call fastBSE_sanity_checks(mpiglobal, input)
        call fastBSE_setup_transitions(mpiglobal, input, fhdf5, ghdf5, unitout)

      ! Taskname 'fastBSE_isdf_cvt'
      case(512)
        call fastBSE_sanity_checks(mpiglobal, input)
        call fastBSE_isdf_cvt(mpiglobal, input, fhdf5, ghdf5, unitout)  

      ! Taskname 'fastBSE_isdf_vexc_test'
      case(513)
        call fastBSE_sanity_checks(mpiglobal, input)
        call fastBSE_isdf_vexc_test(mpiglobal, input, fhdf5, ghdf5, unitout)

      ! Taskname 'xsestimate'
      case(700)
        ! estimate disk-space, cpu-time and memory
        call xsestimate


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

end subroutine xsmain
