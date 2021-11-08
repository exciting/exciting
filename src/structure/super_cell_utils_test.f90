module super_cell_utils_test
   use precision, only: dp
   use modmpi, only: mpiinfo
   use unit_test_framework, only: unit_test_type
   use math_utils, only: all_close
   use super_cell_utils, only: get_translation_vectors, supercell_atomic_positions, &
                               TranslationIntegers_type
   implicit none
   private
   public :: super_cell_utils_test_driver

contains

   !> Run tests for the super_cell_utils module
   subroutine super_cell_utils_test_driver(mpiglobal, kill_on_failure)
      !> mpi environment
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program upon failure of an assertion
      logical, intent(in), optional :: kill_on_failure

      !> Test report object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 9

      call test_report%init(n_assertions, mpiglobal)

      !call functions that test
      call test_get_translation(test_report)
      call test_supercell_atomic_positions(test_report)
      call check_integers_valid(test_report)

      if (present(kill_on_failure)) then
         call test_report%report('super_cell_utils', kill_on_failure)
      else
         call test_report%report('super_cell_utils')
      end if

   end subroutine super_cell_utils_test_driver

   subroutine test_get_translation(test_report)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report

      !> Test input
      !> Lattice vectors
      real(dp) :: lattice_vect(3, 3) = reshape([1.0_dp, 0.0_dp, 1.0_dp, &
                                                0.0_dp, 1.0_dp, 1.0_dp, &
                                                1.0_dp, 1.0_dp, 0.0_dp], &
                                               [3, 3])

      !> Translation integers for zero translation
      type(TranslationIntegers_type) :: n_zero 
      !> Reference for zero translation
      real(dp) :: ref_zero(3, 1)

      !> Translation integers for -1, 1 translation
      type(TranslationIntegers_type) :: n = TranslationIntegers_type([-1, 1], [-1, 1], [-1, 1])
      !> Translation vectors
      real(dp), allocatable :: translation(:, :)
      !> Middle cell ith element of array
      integer :: i_middle_cell
      !> Origin of supercell
      real(dp) :: origin(3)
      
      call test_report%assert(size(get_translation_vectors(lattice_vect, n), 1) == 3, &
                              'Tests translation. Expected result: Dimension one of translation &
                              matrix has a size of 3 for three translation integers.')
      
      call test_report%assert(size(get_translation_vectors(lattice_vect, n), 2) == 27, &
                              'Tests translation. Expected result: Dimension two of translation &
                              matrix has a size of 27 for three translation integers.')
      
      n_zero = TranslationIntegers_type([0, 0], [0, 0], [0, 0])
      ref_zero = 0.0_dp
      call test_report%assert(all_close(get_translation_vectors(lattice_vect, n_zero), ref_zero), &
                              'Tests translation. Expected result: All cells after translation are zero.')
                        
      translation = get_translation_vectors(lattice_vect, n)
      i_middle_cell = int(0.5*size(translation, 2)) + 1
      origin = [0.0_dp, 0.0_dp, 0.0_dp]

      call test_report%assert(all_close(translation(:, i_middle_cell), origin), &
                              'Tests translation. Expected result: Middle cell of translation is zero: [0, 0, 0].')
   end subroutine test_get_translation

   subroutine check_integers_valid(test_report)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report

      !> Test input
      !> Invalid translation integers
      type(TranslationIntegers_type) :: n_integers_invalid = TranslationIntegers_type([1, 5], [1, -2], [1, 5])
      !> Valid translation integers
      type(TranslationIntegers_type) :: n_integers_valid = TranslationIntegers_type([-1, 1], [-1, 1], [-1, 1])

      call test_report%assert(.not. n_integers_invalid%integers_valid(), 'Tests translation integers. &
                              The integers: [1, 5], [1, -2], [1, 5] are not valid.')

      call test_report%assert(n_integers_valid%integers_valid(), 'Tests translation integers. &
                              The integers: [-1, 1], [-1, 1], [-1, 1] are valid.')

   end subroutine check_integers_valid

   subroutine test_supercell_atomic_positions(test_report)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report

      !> Test input
      !> Number of atoms per species
      integer :: natoms(2) = [1, 1]
      !> Translation
      real(dp), parameter :: translation(3, 1) = reshape([-1.0_dp, 2.0_dp, 1.0_dp], &
                                                         [3, 1])
      !> Atomic positions
      real(dp) :: atomic_postion(3, 1, 2) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                     0.5_dp, 0.5_dp, 0.5_dp], [3, 1, 2])
      !> Reference
      real(dp) :: ref(3, 2) = reshape([-1.0_dp, 2.0_dp, 1.0_dp, &
                                       -0.5_dp, 2.5_dp, 1.5_dp], [3, 2])

      call test_report%assert(all_close(supercell_atomic_positions(translation, atomic_postion, natoms), &
                                        ref), 'Tests new positions after translation. Expected: [[-1, 2, 1], [-0.5, 2.5, 1.5]]')

   end subroutine test_supercell_atomic_positions

end module super_cell_utils_test
