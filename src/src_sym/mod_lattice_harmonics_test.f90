!> Unit tests for mod_lattice_harmonics
module mod_lattice_harmonics_test
  use precision, only: dp, i32
  use constants, only: zzero, real_zero, zi, sqrt_two, zone
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, transpose_reshape
  use mod_lattice_harmonics, only: remove_zero_rows, generate_matrix_complex_to_real_spherical_harmonics
  use to_char_conversion, only: to_char

  implicit none

  private
  public :: lattice_harmonics_test_driver

contains

  !> Run tests for [[mod_lattice_harmonics]].
  subroutine lattice_harmonics_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report

    call test_report%init( mpiglobal)

    call test_remove_zero_rows(test_report)
    call test_generate_matrix_complex_to_real_spherical_harmonics(test_report)

    if (present(kill_on_failure)) then
       call test_report%report('lattice_harmonics', kill_on_failure)
    else
       call test_report%report('lattice_harmonics')
    end if

    call test_report%finalise()
  end subroutine lattice_harmonics_test_driver

  !> Test getting removing non-zero rows in a given matrix.
  subroutine test_remove_zero_rows(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: C(4, 3), C_reduced_reference(2, 3)
    real(dp), allocatable :: C_reduced(:, :)
    integer(i32) :: num_non_zero_rows

    C = transpose(reshape([ 1, 2, -3, &
         0, 0, 0, &
         0, 0, 0, &
         -4, 5, 6 ], [3, 4]))

    C_reduced_reference = transpose(reshape([ 1, 2, -3, &
         -4, 5, 6 ], [3, 2]))

    call remove_zero_rows(C, num_non_zero_rows, C_reduced)

    call test_report%assert(num_non_zero_rows == 2, 'Test determining the number of non-zero rows.')
    call test_report%assert(all_close(C_reduced, C_reduced_reference, 1e-8_dp), &
         'Test determining the reduced matrix containing only non-zero rows.')

  end subroutine test_remove_zero_rows

  !> Test getting matrix which transforms complex to real spherical harmonics for \( l \)=2.
  subroutine test_generate_matrix_complex_to_real_spherical_harmonics(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report

    complex(dp) :: A_reference(5, 5)
    complex(dp), allocatable :: A(:, :, :)
    integer(i32) :: lmax, i, j

    lmax = 2

    A_reference = transpose_reshape([ zi / sqrt_two, zzero, zzero, zzero, -zi / sqrt_two, &
         zzero, zi / sqrt_two, zzero, zi / sqrt_two, zzero, &
         zzero, zzero, zone, zzero, zzero, &
         zzero, zone / sqrt_two, zzero, -zone / sqrt_two, zzero, &
         zone / sqrt_two, zzero, zzero, zzero, zone / sqrt_two ], [5, 5])

    call generate_matrix_complex_to_real_spherical_harmonics(lmax, A)

    call test_report%assert(all_close(A(-2:2, -2:2, 2), A_reference, 1e-8_dp), &
         'Test determining matrix which transforms complex to real spherical harmonics for l=2. Error:' // &
         to_char(maxval(abs(A(-2:2, -2:2, 2)-A_reference))) )

  end subroutine test_generate_matrix_complex_to_real_spherical_harmonics

end module mod_lattice_harmonics_test
