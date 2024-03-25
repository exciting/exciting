!> Read input file and call main GW routine.
!>
!> Declarations for the cubic routine, as the quartic data is all passed globally.
subroutine gwtasklauncher()
    use modinput
    use inputdom
    use modmpi, only: terminate, mpiglobal 

    use time_freq_grid, only: grid_type
    use groundstate_results, only: groundstate_results_type
    use space_time, only: space_time_init, space_time_main

    !> Ground state results object, defining all quantities required for GW
    type(groundstate_results_type) :: gs_results
    !> Imaginary time grid
    type(grid_type) :: img_time_grid

    call rereadinput()

    if (input%gw%method == 'quartic') then
        call gw_main()
    else if (input%gw%method == 'cubic') then
        call space_time_init(mpiglobal, gs_results, img_time_grid)
        call space_time_main(mpiglobal, gs_results, img_time_grid)
    else
        call terminate('Invalid GW method selected: ' // input%gw%method)
    end if

end subroutine gwtasklauncher
