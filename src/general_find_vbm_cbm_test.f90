!> Unit tests for the functions in the general_find_vbm_cbm module.

module general_find_vbm_cbm_test
    use precision, only: i32, dp
    use modmpi, only: mpiinfo
    use unit_test_framework, only : unit_test_type
    use general_find_vbm_cbm, only : find_vbm_cbm
  
    implicit none
    private
    public :: general_find_vbm_cbm_test_driver
  
  contains
  
    !> Run tests for find_vbm_cbm procedures
    subroutine general_find_vbm_cbm_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure
      !> Test report object
      type(unit_test_type) :: test_report

      ! Initialize test object
      call test_report%init( mpiglobal)

      ! Run and assert tests
      call test_find_vbm_cbm_efermi(test_report)
      call test_find_vbm_cbm_occupancy(test_report)
  
      ! Report results 
      if (present(kill_on_failure)) then
        call test_report%report('general_find_vbm_cbm', kill_on_failure)
      else
        call test_report%report('general_find_vbm_cbm')
      end if
  
      ! Finalise test object
      call test_report%finalise()
    end subroutine general_find_vbm_cbm_test_driver
  
    !> Test find_vbm_cbm using Fermi energy
    subroutine test_find_vbm_cbm_efermi(test_report)
      type(unit_test_type), intent(inout) :: test_report
  
      integer(i32), parameter :: ib = 1, nb = 4, nk = 3
      real(dp) :: eband(ib:nb, nk)
      real(dp) :: efermi
      integer(i32) :: ibvm, ibcm, ikvm, ikcm, ikvc
      integer(i32) :: expected_ibvm, expected_ibcm, expected_ikvm, expected_ikcm, expected_ikvc
  
      ! Mock data
      eband = reshape([1.0_dp, 2.0_dp, 4.9_dp, 5.3_dp, &  ! k=1
                       1.2_dp, 2.2_dp, 4.7_dp, 5.2_dp, &  ! k=2
                       0.9_dp, 1.9_dp, 4.8_dp, 5.4_dp], & ! k=3
                       [nb, nk])
      efermi = 5.0_dp
  
      call find_vbm_cbm(ib, nb, nk, eband, efermi, ibvm, ibcm, ikvm, ikcm, ikvc)
      expected_ibvm = 3
      expected_ibcm = 4
      expected_ikvm = 1
      expected_ikcm = 2
      expected_ikvc = 1

      call test_report%assert( ibvm == expected_ibvm, 'wrong ibvm!')
      call test_report%assert( ibcm == expected_ibcm, 'wrong ibcm!')
      call test_report%assert( ikvm == expected_ikvm, 'wrong ikvm!')
      call test_report%assert( ikcm == expected_ikcm, 'wrong ikcm!')
      call test_report%assert( ikvc == expected_ikvc, 'wrong ikvc!')
      
    end subroutine test_find_vbm_cbm_efermi
  
    !> Test find_vbm_cbm using occupancy
    subroutine test_find_vbm_cbm_occupancy(test_report)
      type(unit_test_type), intent(inout) :: test_report
  
      integer(i32), parameter :: ib = 1, nb = 4, nk = 3
      real(dp) :: occ(ib:nb, nk), eband(ib:nb, nk)
      integer(i32) :: ibvm, ibcm, ikvm, ikcm, ikvc
      integer(i32) :: expected_ibvm, expected_ibcm, expected_ikvm, expected_ikcm, expected_ikvc
  
      ! Mock data
      occ = reshape([1.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, &
                     1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                     1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], [nb, nk])
      eband = reshape([1.0_dp, 2.0_dp, 5.0_dp, 6.0_dp, &
                       1.5_dp, 2.5_dp, 5.5_dp, 6.5_dp, &
                       0.8_dp, 1.8_dp, 4.8_dp, 5.8_dp], [nb, nk])
  
      call find_vbm_cbm(ib, nb, nk, occ, eband, ibvm, ibcm, ikvm, ikcm, ikvc)

      expected_ibvm = 3
      expected_ibcm = 3
      expected_ikvm = 1
      expected_ikcm = 3
      expected_ikvc = 1

      call test_report%assert( ibvm == expected_ibvm, 'wrong ibvm!')
      call test_report%assert( ibcm == expected_ibcm, 'wrong ibcm!')
      call test_report%assert( ikvm == expected_ikvm, 'wrong ikvm!')
      call test_report%assert( ikcm == expected_ikcm, 'wrong ikcm!')
      call test_report%assert( ikvc == expected_ikvc, 'wrong ikvc!')

    end subroutine test_find_vbm_cbm_occupancy
  
  end module general_find_vbm_cbm_test
  
