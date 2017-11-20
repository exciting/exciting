! Copyright (C) 2009-2010 S. Sagmeister, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xstasklauncher
  use modinput
  use modxs, only: temat, skipgnd, hybridhf
  use inputdom
  use mod_hdf5

  implicit none

  integer(4) :: nxstasks, nxstasksmax, i
  real(8), parameter :: eps=1.d-7

  ! SET DEFAULTS
  ! Set the default values if "TDDFT" element is not present
  if( .not. associated(input%xs%tddft) ) then
    input%xs%tddft => getstructtddft(emptynode)
  end if

  ! Set the default values if "screening" element is not present
  if( .not. associated(input%xs%screening) ) then
    input%xs%screening => getstructscreening(emptynode)
  end if

  ! Set the default values if "BSE" element is not present
  if( .not. associated(input%xs%bse) ) then
    input%xs%bse => getstructbse(emptynode)
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

  ! Use specified plan or construct default plans

  if(associated(input%xs%plan)) then

    nxstasks = size(input%xs%plan%doonlyarray)
    call xsmain(input%xs%plan, nxstasks)

  else if(trim(input%xs%xstype) .eq. "TDDFT") then

    ! Allocate plan
    nxstasks = 0
    nxstasksmax = 10
    allocate(input%xs%plan)
    allocate(input%xs%plan%doonlyarray(nxstasksmax))
    do i = 1, nxstasksmax
      allocate(input%xs%plan%doonlyarray(i)%doonly)
    end do

    ! Setup default plan
    if(input%xs%tddft%do .eq. "fromscratch") then

      ! Task 301 corresponds to "xsgeneigvec" plan
      ! One shot GS calculation with xs%ngridk, xs%nempty and potential xs%vkloff.
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="xsgeneigvec"

      ! Task 320 corresponds to "writepmatxs" plan
      ! Calculates the momentum matrix elements for the xs GS calculation.
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="writepmatxs"

      ! Task 330 corresponds to "writeemat" plan
      ! Calculates the plane wave matrix elements, is skipped when
      ! gqmax = 0 and only gamma point is considered
      temat=.true.
      if( (size(input%xs%qpointset%qpoint, 2) .eq. 1) .and. (input%xs%gqmax .lt. eps)) then
        if(sum(abs(input%xs%qpointset%qpoint(:, 1))) .lt. eps) then 
          temat = .false.
        end if
      end if
      if(temat) then
        nxstasks = nxstasks+1
        input%xs%plan%doonlyarray(nxstasks)%doonly%task="writeemat"
      end if

      ! BSE derived kernels ?
      ! 7 = "MB1_NLF", 8 = "BO"
      if(input%xs%tddft%fxctypenumber .eq. 7 .or. &
         & input%xs%tddft%fxctypenumber .eq. 8) then

        ! Task 401 corresponds to "scrgeneigvec" plan
        ! One shot GS calculation with more empty states xs%screening%nempty 
        ! but otherwise identical parameters as "xsgeneigvec".
        nxstasks = nxstasks+1
        input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrgeneigvec"

        ! Task 420 corresponds to "scrwritepmat" plan
        nxstasks = nxstasks+1
        input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrwritepmat"

        if(input%xs%screening%do .eq. "fromscratch") then

          ! Task 430 corresponds to "screen" plan
          ! Generate KS RPA screening
          nxstasks = nxstasks+1
          input%xs%plan%doonlyarray(nxstasks)%doonly%task="screen"

          ! Task 440 corresponds to "scrcoulint" plan
          ! Generate screened Coulomb interaction matrix
          nxstasks = nxstasks+1
          input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrcoulint"

        end if

        ! Task 450 corresponds to "kernxs_bse" plan
        nxstasks = nxstasks+1
        input%xs%plan%doonlyarray(nxstasks)%doonly%task="kernxc_bse"

      end if

      ! Task 340 corresponds to "df" plan
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="df"

      ! Task 350 corresponds to "idf" plan
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="idf"

    else

      ! Task 350 corresponds to "idf" plan
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="idf"

    end if

    ! Set associated task numbers
    do i = 1, nxstasks
     input%xs%plan%doonlyarray(i)%doonly%tasknumber =&
      stringtonumberdoonlytask(input%xs%plan%doonlyarray(i)%doonly%task)
    end do

    ! Execute plan
    call xsmain(input%xs%plan, nxstasks)

  else if(trim(input%xs%xstype) .eq. "BSE") then

    ! Allocate plan
    nxstasksmax = 10
    allocate(input%xs%plan)
    allocate(input%xs%plan%doonlyarray(nxstasksmax))
    do i = 1, nxstasksmax
      allocate(input%xs%plan%doonlyarray(i)%doonly)
    end do

    ! Setup default plan

    ! Task 301 corresponds to "xsgeneigvec" plan
    ! One shot GS calculation with xs%ngridk, xs%nempty and potential xs%vkloff.
    nxstasks=1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="xsgeneigvec"
    ! Task 320 corresponds to "writepmatxs" plan
    ! Calculates the momentum matrix elements for the xs GS calculation.
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="writepmatxs"

    ! Task 401 corresponds to "scrgeneigvec" plan
    ! One shot GS calculation with more empty states xs%screening%nempty 
    ! but otherwise identical parameters as "xsgeneigvec".
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrgeneigvec"
    ! Task 420 corresponds to "scrwritepmat" plan
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrwritepmat"

    if(input%xs%screening%do .eq. "fromscratch") then
      ! Task 430 corresponds to "screen" plan
      ! Generate KS RPA screening
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="screen"
      ! Task 440 corresponds to "scrcoulint" plan
      ! Generate screened Coulomb interaction matrix
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrcoulint"
    end if

    ! Task 441 corresponds to "exccoulint" plan
    ! Generate unscreened Coulomb exchange interaction matrix
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="exccoulint"

    ! Task 445 corresponds to "bse" plan
    ! Set up and solve BSE
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="bse"

    ! Set associated taks numbers
    do i = 1, nxstasks
     input%xs%plan%doonlyarray(i)%doonly%tasknumber =&
      stringtonumberdoonlytask(input%xs%plan%doonlyarray(i)%doonly%task)
    end do

    ! Execute plan
    call xsmain(input%xs%plan, nxstasks)

  else

     write (*,*) "error xstasklauncher"
     write (*,*) trim (input%xs%xstype), "no valid xstype"
     stop

  end if

end subroutine
