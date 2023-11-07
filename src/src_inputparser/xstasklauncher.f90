! Copyright (C) 2009-2010 S. Sagmeister, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Modified May 2019 (Ronaldo) - to include RT-TDDFT
! Modified Dec 2020 (Ronaldo) - to improve the interface with RT-TDDFT
module xs_task_launcher
  use modinput, only: input, &
                      getstructtetra, &
                      getstructbse, &
                      getstructfastBSE, &
                      getstructtddft, &
                      getstructscreening, &
                      stringtonumberdoonlytask
  use modxs, only: temat, skipgnd
  use mod_hybrids, only: hybridhf
  use inputdom
  use rttddft_main, only: run_rttddft => coordinate_rttddft_calculation
  use errors_warnings, only: terminate_if_false
  use modmpi, only: mpiglobal
  use xs_plan_manager, only: xs_plan_bse, &
                              xs_plan_tddft, &
                              xs_plan_fastBSE
  use fastBSE, only: fastBSE_sanity_checks

  implicit none

  contains

  subroutine xstasklauncher

    integer(4) :: nxstasks
    character(:), allocatable :: message

    ! Check if RT-TDDFT is desired
    if( trim( input%xs%xstype ) /= "RT-TDDFT") then
      call terminate_if_false( mpiglobal, associated(input%xs%energywindow), &
        & 'ERROR in xs: an energywindow is required!' )
      call terminate_if_false( mpiglobal, associated(input%xs%qpointset), &
        & 'ERROR in xs: a qpointset is required!' )
    end if

    ! Set defaults in the xs part of the input global
    call initialize_xs_input()

    ! User defined plan
    if(associated(input%xs%plan)) then
      nxstasks = size(input%xs%plan%doonlyarray)

    ! RT-TDDFT calculation
    ! TODO(Ronaldo)
    else if ( trim( input%xs%xstype ) == "RT-TDDFT" ) then
      ! We know that xstype is RT-TDDFT
      ! But we need to check if the element input%xs%rt_tddft has been defined
      call terminate_if_false( mpiglobal, associated( input%xs%realTimeTDDFT ), &
        & 'ERROR in RT-TDDFT: you need to add the element rt_tddft inside xs in input.xml!')
      call run_rttddft
      return

    ! TDDFT calculation
    else if (trim(input%xs%xstype) == "TDDFT") then
      call xs_plan_tddft(nxstasks)

    
    else if(trim(input%xs%xstype) == "BSE") then
      ! BSE calculation
      if(trim(input%xs%BSE%solver) == "direct") then
        call xs_plan_bse(nxstasks)

      ! fastBSE calculation
      else if(trim(input%xs%BSE%solver) == "fastBSE") then
        call fastBSE_sanity_checks(mpiglobal, input)
        call xs_plan_fastBSE(nxstasks)
      end if
    else
      ! Stop the code: xstype not recognized
      message = 'error xstasklauncher: ' // trim(input%xs%xstype) // ' is not a valid xstype.'
      call terminate_if_false( mpiglobal, .false., message )

    end if

    call xsmain(input%xs%plan, nxstasks)

  end subroutine xstasklauncher

  !> Set defaults for the xs element of the input instance.
  subroutine initialize_xs_input()
    ! Set the default values if "TDDFT" element is not present
    if( .not. associated(input%xs%tddft) ) then
      input%xs%tddft => getstructtddft(emptynode)
    end if

    ! Set the default values if "screening" element is not present
    if( .not. associated(input%xs%screening) ) then
      input%xs%screening => getstructscreening(emptynode)
    end if

    ! Set the default values if "BSE" element is not present
    if( .not. associated(input%xs%BSE) ) then
      input%xs%BSE => getstructbse(emptynode)
    end if

    ! Set the default values if "fastBSE" is solver but element is not present
    if( .not. associated(input%xs%fastBSE) ) then
      input%xs%fastBSE => getstructfastBSE(emptynode)
    end if


    ! Set the default values if "tetra" element not present
    if( .not. associated(input%xs%tetra) ) then
      input%xs%tetra => getstructtetra(emptynode)
    end if

    ! If "skipgnd" is set to true ks eigenvalues and eigenvectors are not recalculated
    skipgnd=input%xs%skipgnd

    if(associated(input%groundstate%hybrid)) then
      if(input%groundstate%hybrid%exchangetypenumber == 1) then
        hybridhf = .true.
        skipgnd = .true.
      end if
    end if

  end subroutine initialize_xs_input

end module xs_task_launcher