module mod_getoccsv
#include "asserts.fpp"
  use mod_kpoint, only: nkpt, vkl
  use mod_names, only: filetag_occsv
  use modinput, only: input
  use modmpi, only: firstofset, procofindex, splittfile, terminate
  use precision, only: dp, i32, str_256

  implicit none
  private
  integer(i32), parameter :: n_cartesian = 3
  integer(i32), parameter :: max_number_of_read_attempts = 100

  public :: getoccsv

contains

  !> Read the occupation numbers from `OCCSV.OUT`.
  !> For spin-unpolarized systems, the maximum occupancy of a state is 2; for
  !> spin-polarized systems, it is 1.
  !> `OCCSV.OUT` is a direct-access binary file. Its record length can be
  !> determined from the array sizes and data-type specifications.
  !> Each record of this file has the following structure
  !> 
  !> | \( k_{\rm lat} \)  | \( N_{\rm stsv} \) | occupation |
  !> 
  !> The meaning of each component is described in the table below:
  !>
  !> | name | type | shape | description |
  !> | ----------- | ----------- | ----------- | ----------- |
  !> | \( k_{\rm lat} \) | real(dp) | 3 | \( \mathbf{k} \)-point in lattice coordinates |
  !> | \( N_{\rm stsv} \) | integer | 1 | number of (second-variational) states (without core states) |
  !> | occupation | real(dp) | 1 | (second-variational) occupation number |
  subroutine getoccsv( kpt_in_lattice_coordinates, occupations )
    !> \(\mathbf{k}\)-point in lattice coordinates
    real(dp), intent(in) :: kpt_in_lattice_coordinates(:)
    !> Occupation factors of each KS state for the given \(\mathbf{k}\)-point
    real(dp), intent(out) :: occupations(:)
    
    logical :: occupations_file_exist
    integer(i32) :: isym, ik, koffset, i, recl, nstates_from_file, occsv_size, file_unit
    real(dp) :: kpt_from_file(n_cartesian), error_kpt_coordinates
    character(len=:), allocatable :: file_name
    character(len=str_256), external :: outfilenamestring
    real(dp), allocatable :: occupations_from_file(:)
    
    call findkpt( kpt_in_lattice_coordinates, isym, ik )  
    file_name = ( outfilenamestring(trim( filetag_occsv ), ik) )
    occsv_size = size( occupations )
    inquire (IoLength=Recl) kpt_from_file, nstates_from_file
    inquire ( file=file_name, exist=occupations_file_exist )
    if ( .not. occupations_file_exist ) call wait_processes_to_write_occsv( file_name, occupations_file_exist )
    if ( occupations_file_exist ) then
      open( newunit=file_unit, file=file_name, Action='READ', Form='UNFORMATTED', access='DIRECT', Recl=Recl )
    else
      call terminate( "The file containing occupations is missing." )
    end if
    koffset = ik
    if ( splittfile ) koffset = koffset - firstofset(procofindex(ik, nkpt), nkpt) + 1
    read (file_unit, Rec=1) kpt_from_file, nstates_from_file
    close ( file_unit )
    call sanity_check_n_states( occsv_size, nstates_from_file, ik, file_name )
    allocate( occupations_from_file(nstates_from_file) )
    inquire (IoLength=Recl) kpt_from_file, nstates_from_file, occupations_from_file
    open (newunit=file_unit, file=file_name, Action='READ', Form='UNFORMATTED', access='DIRECT', Recl=Recl)
    read ( file_unit, Rec=koffset ) kpt_from_file, nstates_from_file, occupations_from_file
    close ( file_unit )
    call sanity_check_kpt_coordinates( vkl(:, ik), kpt_from_file, input%structure%epslat, ik, file_name )
    occupations (:) = occupations_from_file(:occsv_size)
  end subroutine getoccsv

  !> (private) Wait for other processes to write the occupations file.
  subroutine wait_processes_to_write_occsv( file_name, occupations_file_exist )
    !> Name of the file containing occupations
    character(len=*), intent(in) :: file_name
    !> If .true., the file exsists
    logical, intent(out) :: occupations_file_exist

    integer(i32) :: i

    do i = 1, max_number_of_read_attempts
      inquire ( file=file_name, exist=occupations_file_exist )
      if ( occupations_file_exist ) then
        exit
      else
        call system( 'sync' )
        write (*,*) "Waiting for other process to write" // ":getoccsv:" // file_name
        call sleep(5)
      end if
    end do
  end subroutine

  !> (private) Check whether the occupations file contain enough states.
  subroutine sanity_check_n_states( n_states_requested, max_n_states_expected, ik, file_name )
    !> Number of requested states
    integer(i32), intent(in) :: n_states_requested
    !> Number of states in the occupations file
    integer(i32), intent(in) :: max_n_states_expected
    !> \mathbf{k} point index
    integer(i32), intent(in) :: ik
    !> Name of the file containing occupations
    character(len=*), intent(in) :: file_name

    if ( n_states_requested > max_n_states_expected ) then
      write (*,*)
      write (*, '("Error(getoccsv): invalid nstfv for k-point ", I8)') ik
      write (*, '(" current    : ", I8)') n_states_requested
      write (*, '(" OCCSV.OUT  : ", I8)') max_n_states_expected
      write (*, '(" file	     : ", a)') file_name
      write (*,*)
      call terminate()
    end if
  end subroutine

  !> (private) Check whether the requested \mathbf{k} point corresponds to the one from the occupations file.
  subroutine sanity_check_kpt_coordinates( kpt_requested, kpt_expected, tol, ik, file_name )
    !> Requested \mathbf{k} point in lattice coordinates
    real(dp), intent(in) :: kpt_requested(:)
    !> \mathbf{k} point in lattice coordinates in the occupations file
    real(dp), intent(in) :: kpt_expected(:)
    !> \mathbf{k} point index
    integer(i32), intent(in) :: ik
    !> Tolerance used to compare \mathbf{k} vectors
    real(dp), intent(in) :: tol
    !> Name of the file containing occupations
    character(len=*), intent(in) :: file_name

    real(dp) :: error_kpt_coordinates

    CALL_ASSERT( size( kpt_requested ) == size( kpt_expected ), "Sizes of kpt_requested and kpt_expected don't match." )
    error_kpt_coordinates = sum( abs( kpt_requested - kpt_expected ) )
    if ( error_kpt_coordinates > tol ) then
      write (*,*)
      write (*, '("Error(getoccsv): differing vectors for k-point ", I8)') ik
      write (*, '(" current    : ", 3G18.10)') kpt_expected
      write (*, '(" OCCSV.OUT  : ", 3G18.10)') kpt_requested
      write (*, '(" file	     : ", a)') file_name
      write (*,*)
      call terminate()
    end if
  end subroutine
end module
