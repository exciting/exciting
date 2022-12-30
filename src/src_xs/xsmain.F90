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
  use mod_hdf5
  use m_write_hdf5, only: fhdf5_inter
  use phonon_screening, only: phonon_screening_launcher
  use write_screening, only: write_screening_launcher

  implicit none

  !> Screening from polar phonons
  integer, parameter :: task_phonon_screening = 431
  !> Screened Coulomb interaction for BSE
  integer, parameter :: task_screened_coulomb = 440 
  !> Writing dielectric matrix for all non-reduced q-vectors
  integer, parameter :: task_write_dielectric_matrix = 442
  !> Writing screened Coulomb matrix for all non-reduced q-vectors
  integer, parameter :: task_write_screened_coulomb = 443

  type(plan_type), intent(in) :: plan
  integer(4), intent(in) :: nxstasks
  logical :: fex
  integer(4) :: i

  ! initialization of hdf5 output
#ifdef _HDF5_
  if(mpiglobal%rank == 0) then
    call hdf5_initialize()
    fhdf5="bse_output.h5"
    ! find out whether file already exists
    inquire(file=trim(fhdf5), exist=fex)
    if (.not. fex) then
      call hdf5_create_file(fhdf5)
    end if
    fhdf5_inter='bse_matrix.h5'
    call hdf5_create_file(fhdf5_inter)
  end if
#endif

  do i = 1, nxstasks
     task = plan%doonlyarray(i)%doonly%tasknumber
  end do

  
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
      
      ! Taskname 'phonon_screening'
      case(task_phonon_screening) 

        call phonon_screening_launcher

      ! Taskname 'scrcoulint'
      case(task_screened_coulomb)
        ! screened Coulomb interaction
        call scrcoulintlauncher

      ! Taskname 'exccoulint'
      case(441)
        ! exchange Coulomb interaction
        call exccoulintlauncher

      ! Taskname 'write_dielectric_matrix'
      case (task_write_dielectric_matrix)
          ! Expanding dielectric matrix
          call write_screening_launcher('write_dielectric_matrix')

          ! Taskname 'write_screened_coulomb'
      case (task_write_screened_coulomb)
          ! Expanding dielectric matrix
          call write_screening_launcher('write_screened_coulomb')

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

      ! Taskname 'screen'
      case(451)
        ! write real-space XS wfcts to file
        call write_wfplot

      ! Taskname 'write_screen'
      case(452)
        ! write screened Coulomb potential to file
        call write_screen

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


  ! Finalization of hdf5 output
#ifdef _HDF5_
  if (mpiglobal%rank == 0) then
    call hdf5_finalize()
  end if
#endif

end subroutine xsmain
