!> This module contains the tests for the workspaces
module workspaces_test

    use precision, only: dp, i32
    use modmpi, only: mpiinfo
    use unit_test_framework, only: unit_test_type
    use lapack_workspaces, only: lapack_workspace_complex_dp_t

    private
    public :: workspaces_test_driver

contains
    
    !> Run tests for hermitian eigensolvers
    subroutine workspaces_test_driver(mpiglobal, kill_on_failure)
        !> mpi environment
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program upon failure of an assertion
        logical, intent(in), optional :: kill_on_failure

        !> Test report object
        type(unit_test_type) :: test_report

        call test_report%init( mpiglobal)

        call test_workspace(test_report)
        
        if (present(kill_on_failure)) then
            call test_report%report('lapack_workspaces', kill_on_failure)
        else
            call test_report%report('lapack_workspaces')
        end if

        call test_report%finalise()
    end subroutine workspaces_test_driver


    !> Check workspace class
    subroutine test_workspace(test_report)
        !> Test report
        type(unit_test_type) :: test_report
        
        type(lapack_workspace_complex_dp_t) :: my_space

        ! Check if when inited all is ok
        call test_report%assert(.not. my_space%computed(), &
                        'Test if lapack_workspace_t is properly inited &
                        Expected: Default inited workspace.')
        
        ! Check if workspace sizes are given
        call my_space%initialize(lrwork=42, lwork=42)
        
        call test_report%assert(my_space%computed(), &
                        'Test if lapack_workspace_t%computed and lapack_workspace_t%initialize works (1) &
                        Expected: A valid workspace with sizes.')

        ! Check reset
        call my_space%reset()
        call test_report%assert(.not. my_space%computed(), &
                        'Test if lapack_workspace_t%reset works &
                        Expected: That the workspace is back to the default state (i.e. not computed).')
        
    end subroutine test_workspace


end module workspaces_test
