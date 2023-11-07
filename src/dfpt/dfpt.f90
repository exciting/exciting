!> This is the main module for organizing and launching a
!> density-functional perturbation theory (DFPT) calculation.
!> 
!> At the moment, the only perturbations are nuclei displacements
!> which allow for the calculation of lattice vibrations.
!> In the future, also the response w.r.t. an external electric field
!> will be available.
module dfpt
  use dfpt_variables
  use dfpt_density_potential
  use dfpt_eigensystem
  use dfpt_inout

  use precision, only: dp

  implicit none
  private

  public :: dfpt_launcher

  contains

    !> Main launcher for a DFPT calculation.
    !>
    !> This routine is setting up the different tasks that need to be done
    !> according to the input file and starts the execution of these tasks.
    subroutine dfpt_launcher
      use phonons
      use phonons_variables, only: ph_var_init, ph_parts_per_rank
      use efield
      use efield_variables, only: ef_var_init

      use os_utils, only: path_exists
      use modmpi
      use modinput

      integer :: max_num_tasks, num_tasks
      integer :: i
      character(64) :: task

      character(64), allocatable :: task_list(:)

      ! ** CREATE LIST OF TASKS

      ! generate list of tasks that should be done
      max_num_tasks = 16
      if( associated( input%phonons%parts) ) &
        max_num_tasks = size( input%phonons%parts%dopartarray )
      allocate( task_list(max_num_tasks) )

      ! dry run
      if( input%phonons%do == 'dry' ) then
        task_list(1) = 'dry_run'
        num_tasks = 1
      ! read tasks from file
      else if( associated( input%phonons%parts ) ) then
        num_tasks = 0
        do i = 1, max_num_tasks
          task = trim( adjustl( input%phonons%parts%dopartarray(i)%dopart%id ) )
          select case( task )
            ! generate eigenvalues and eigenvectors of unperturbed system
            ! and write files EVALK0.TMP and EVECK0.TMP
            case( 'gen_eig0' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'gen_eig0'
            ! delete temporary files
            case( 'cleanup' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'cleanup'
            ! solve self-consistency cycle for density and potential response
            ! upon a phonon-like perturbation
            case( 'phonons_scf' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'phonons_scf'
            ! compute force response upon a phonon-like perturbation
            case( 'phonons_forces' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'phonons_forces'
            ! compute polarization response upon a phonon-like perturbation
            case( 'phonons_polarization' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'phonons_polarization'
            ! compute dynamical matrices in canonical basis
            case( 'phonons_dynmat' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'phonons_dynmat'
            ! compute Born effective charge tensors in canonical basis
            case( 'phonons_borncharge' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'phonons_borncharge'
            ! write potential response to binary file
            case( 'phonons_write_dpot' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'phonons_write_dpot'
            ! solve self-consistency cycle for density and potential response
            ! upon a perturbing electric field
            case( 'efield_scf' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'efield_scf'
            ! compute polarization response upon a perturbing electric field
            case( 'efield_polarization', 'efield_dielten' )
              num_tasks = num_tasks + 1
              task_list(num_tasks) = 'efield_polarization'
            case default
              write( *, '(*(a))' ) "Warning(dfpt_launcher): Invalid task '", trim( task), "'."
          end select
        end do
      ! set up default task list
      else
        task_list(1) = 'gen_eig0'
        task_list(2) = 'cleanup'
        task_list(3) = 'phonons_scf'
        task_list(4) = 'phonons_forces'
        task_list(5) = 'phonons_dynmat'
        num_tasks = 5
        if( input%phonons%polar ) then
          task_list(6) = 'phonons_polarization'
          task_list(7) = 'phonons_borncharge'
          task_list(8) = 'efield_scf'
          task_list(9) = 'efield_polarization'
          num_tasks = num_tasks + 4
        end if
      end if

      ! ** EXECUTE LIST OF TASKS

      ! initialize global DFPT variables
      call dfpt_var_init

      ! perform phonon dry run
      if( any( task_list == 'dry_run' ) ) then
        call dfpt_io_set_prefix( 'PHONON', 'phonon' )
        call ph_dry_run
      end if

      ! initialize density and potential variables
      if( any( task_list == 'gen_eig0' ) .or. &
          any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) .or. &
          any( task_list == 'efield_scf' ) .or. &
          any( task_list == 'efield_polarization' ) ) &
        call dfpt_rhopot_init

      ! initialize eigensystem variables
      ! and compute eigenvalues and eigenvectors of unperturbed system
      ! and write files EVALK0.TMP and EVECK0.TMP
      if( any( task_list == 'gen_eig0' ) .or. &
          any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) .or. &
          any( task_list == 'efield_scf' ) .or. &
          any( task_list == 'efield_polarization' ) ) then
        call dfpt_eig_init
        call dfpt_prepare
      end if

! -------START EFIELD RESPONSE-------------------------------------------------

      if( any( task_list == 'efield_scf' ) .or. &
          any( task_list == 'efield_polarization' ) ) then
        ! initialize global DFPT electric field variables
        call ef_var_init
        ! prepare electric field calculation
        call dfpt_io_set_prefix( 'EFIELD', 'efield' )
        call ef_prepare
        ! solve self-consistency cycle for density and potential response
        if( any( task_list == 'efield_scf' ) ) then
          if( .not. path_exists( 'EPSINF.OUT' ) ) call ef_scf
        end if
        ! compute polarization response
        ! (dielectric tensor)
        if( any( task_list == 'efield_polarization' ) ) then
          if( .not. path_exists( 'EPSINF.OUT' ) ) call ef_polarization
        end if
        ! finalize electric field calculation
        call ef_finalize( delete_files=input%phonons%delete_eigensystem_response )
      end if

! -------END EFIELD RESPONSE---------------------------------------------------

! -------START PHONON RESPONSE-------------------------------------------------

      ! initialize global DFPT phonons variables
      if( any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) .or. &
          any( task_list == 'phonons_dynmat' ) .or. &
          any( task_list == 'phonons_borncharge' ) .or. &
          any( task_list == 'phonons_write_dpot' ) ) &
        call ph_var_init

      ! prepare phonon calculation
      if( any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) ) then
        call dfpt_io_set_prefix( 'PHONON', 'phonon' )
        call ph_prepare( do_force=any( task_list == 'phonons_forces' ) )
      end if

      ! execute independent phonons parts
      if( any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) ) then
        do i = 1, size( ph_parts_per_rank )
          ! prepare independent part
          call ph_part_prepare( ph_parts_per_rank(i) )
          ! solve self-consistency cycle for density and potential response
          if( any( task_list == 'phonons_scf' ) ) &
            call ph_part_scf( ph_parts_per_rank(i) )
          ! compute force response
          ! (dynamical matrix column in irrep coordinates)
          if( any( task_list == 'phonons_forces' ) ) &
            call ph_part_force( ph_parts_per_rank(i) )
          ! compute polarization response
          ! (Born effective charge column in irrep coordinates)
          if( any( task_list == 'phonons_polarization' ) ) &
            call ph_part_polarization( ph_parts_per_rank(i) )
          ! finalize independent part
          call ph_part_finalize( ph_parts_per_rank(i) )
        end do
      end if
      call barrier( mpicom=mpiglobal )

      ! compute dynamical matrices in canonical basis
      ! and write them to file
      if( any( task_list == 'phonons_dynmat' ) ) &
        call ph_write_dyn_canonical
      call barrier( mpicom=mpiglobal )

      ! compute effective potential response in canonical basis
      ! and write it to file
      if( any( task_list == 'phonons_write_dpot' ) ) &
        call ph_write_dpot_canonical
      call barrier( mpicom=mpiglobal )

      ! compute Born effective charge tensors in canonical basis
      ! and write them to file
      if( any( task_list == 'phonons_borncharge' ) ) &
        call ph_write_borncharge_canonical
      call barrier( mpicom=mpiglobal )

      ! finalize phonon calculation
      if( any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_dynmat' ) .or. &
          any( task_list == 'phonons_polarization' ) ) &
        call ph_finalize

! -------END PHONON RESPONSE---------------------------------------------------

      ! close and delete files EVALK0.TMP and EVECK0.TMP
      if( any( task_list == 'gen_eig0' ) .or. &
          any( task_list == 'cleanup' ) .or. &
          any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) .or. &
          any( task_list == 'efield_scf' ) .or. &
          any( task_list == 'efield_polarization' ) ) &
        call dfpt_finalize( any( task_list == 'cleanup') )

      ! delete eigensystem variables
      if( any( task_list == 'gen_eig0' ) .or. &
          any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) ) &
        call dfpt_eig_free

      ! delete density and potential variables
      if( any( task_list == 'gen_eig0' ) .or. &
          any( task_list == 'phonons_scf' ) .or. &
          any( task_list == 'phonons_forces' ) .or. &
          any( task_list == 'phonons_polarization' ) ) &
        call dfpt_rhopot_free

      ! delete global DFPT variables
      call dfpt_var_free
    end subroutine dfpt_launcher

    !> This subroutine executes preparative tasks that are common to all 
    !> DFPT calculations.
    !>
    !> This inculdes the calculation of all unperturbed Kohn-Sham eigenstates
    !> at the (reduced) set of electronic wavevectors \({\bf k}\). These are 
    !> written to the temporary binary files `EVALK0.TMP` and `EVECK0.TMP` 
    !> which are autSmatically deleted at the end of the calculation.
    subroutine dfpt_prepare
      use modmpi
      use constants, only: zone
      use mod_APW_LO, only: nlotot
      use block_data_file, only: block_data_file_type

      integer :: nmatmax, ik, ik1, ik2
      logical :: exists

      real(dp), allocatable :: evalk(:,:)
      complex(dp), allocatable :: eveck(:,:)

      ! set limits for k-point loops
      ik1 = firstofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
      ik2 = lastofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
      ! * We need all (unoccupied) eigenpairs at each k-point in the k-set.
      ! * We compute them once and store them on disk.
      ! * The file is deleted at the end of the calculation.
      ! allocate local variables
      nmatmax = dfpt_Gkset%ngkmax + nlotot
      allocate( evalk(nmatmax, dfpt_kset%nkpt), source=0.0_dp )
      allocate( eveck(nmatmax, nmatmax) )
      ! open temporary file for eigenvectors and eigenvalues
      fevalk0 = block_data_file_type( 'EVALK0.TMP', [nmatmax], 1.d0 )
      feveck0 = block_data_file_type( 'EVECK0.TMP', [nmatmax, nmatmax], zone )
      exists = fevalk0%exists() .and. feveck0%exists()
      call barrier
      call fevalk0%open( mpiglobal )
      call feveck0%open( mpiglobal )
      ! loop over k-points
      do ik = ik1, ik2
        ! solve KS equation
        if( exists ) then
          call dfpt_eig_ks( ik, dfpt_kset, dfpt_Gset, dfpt_Gkset, nmatmax, evalk(:,ik), eveck, &
                 p0set=dfpt_kset, Gp0set=dfpt_Gkset, feval=fevalk0, fevec=feveck0 )
        else
          call dfpt_eig_ks( ik, dfpt_kset, dfpt_Gset, dfpt_Gkset, nmatmax, evalk(:,ik), eveck )
          ! write eigenvectors and eigenvalues to file
          call fevalk0%write( ik, evalk(:, ik) )
          call feveck0%write( ik, eveck )
        end if
      end do
      deallocate( evalk, eveck )
      call barrier
    end subroutine dfpt_prepare

    !> This subroutine executes finalizing tasks that are common to all 
    !> DFPT calculations.
    !>
    !> This includes freeing memory from unneeded variables, closing files
    !> and deleting temporary files.
    subroutine dfpt_finalize( cleanup )
      use modmpi, only: mpiglobal
      !> delete eigenvector and eigenvalues files
      logical, intent(in) :: cleanup
      ! close temporary files for eigenvectors and eigenvalues and delete them
      call fevalk0%close( mpiglobal )
      call feveck0%close( mpiglobal )
      if( cleanup ) then
        call fevalk0%delete( mpiglobal )
        call feveck0%delete( mpiglobal )
      end if
    end subroutine dfpt_finalize

end module dfpt
