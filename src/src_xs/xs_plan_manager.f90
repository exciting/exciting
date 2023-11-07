module xs_plan_manager
  use inputdom
  use modinput, only: input, stringtonumberdoonlytask
  use modxs, only: temat
  use errors_warnings, only: terminate_if_false
  use modmpi, only: mpiglobal
  use bse_utils, only: bse_type_to_bool
  

  implicit none

  !character(*), dimension(:), parameter ::  


  contains

  !> Setup default plan for BSE calculation in the input global.
  subroutine xs_plan_bse(nxstasks)
    !> Number of tasks
    integer, intent(out) :: nxstasks

    integer :: nxstasksmax, i

     ! Allocate plan
    nxstasksmax = 10
    allocate(input%xs%plan)
    allocate(input%xs%plan%doonlyarray(nxstasksmax))
    do i = 1, nxstasksmax
      allocate(input%xs%plan%doonlyarray(i)%doonly)
    end do

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
      ! Task 431 corresponds to "phonon_screening" plan
      ! Generate phonon screening
      if( associated(input%xs%phonon_screening) ) then
        nxstasks = nxstasks+1
        input%xs%plan%doonlyarray(nxstasks)%doonly%task="phonon_screening"
      end if
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
      input%xs%plan%doonlyarray(i)%doonly%tasknumber = &
              stringtonumberdoonlytask(input%xs%plan%doonlyarray(i)%doonly%task)
    end do
    
  end subroutine xs_plan_bse

  !> Setup default plan for TDDFT calculation in the input global.
  subroutine xs_plan_tddft(nxstasks)
    !> Number of tasks
    integer, intent(out) :: nxstasks

    integer :: nxstasksmax, i
    real(8), parameter :: eps=1.d-7

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
  
  end subroutine xs_plan_tddft

  
  !> Setup default plan for fast BSE calculation in the input global.
  subroutine xs_plan_fastBSE(nxstasks)
    !> Number of tasks
    integer, intent(out) :: nxstasks

    integer :: nxstasksmax, i

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


    ! Task 321 corresponds to "write_pmatrix_hdf5" plan
    ! ASCII output of momentum matrix elements
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="writepmatxs"

    ! Task 451 corresponds to "write_wfplot" plan
    ! Write wavfunctions to hdf5 file
    nxstasks = nxstasks + 1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="write_wfplot"


    ! Allows for skipping screening calculation
    if(input%xs%screening%do .eq. "fromscratch") then
      ! Task 401 corresponds to "scrgeneigvec" plan
      ! One shot GS calculation with more empty states xs%screening%nempty
      ! but otherwise identical parameters as "xsgeneigvec".
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrgeneigvec"

      ! Task 420 corresponds to "scrwritepmat" plan
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="scrwritepmat"

      ! Task 430 corresponds to "screen" plan
      ! Generate KS RPA screening
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="screen"

      ! Task 452 corresponds to "write_screened_coulomb" plan
      ! Write screened Coulomb potential to file
      nxstasks = nxstasks+1
      input%xs%plan%doonlyarray(nxstasks)%doonly%task="write_screened_coulomb"
    end if

     ! Task 510 corresponds to "isdf_lanczos_bse" plan
    ! Set up and solve BSE with ISDF compression of the BSH and Lanczos solver
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="fastBSE_setup_transitions"

     ! Task 512 corresponds to "isdf_lanczos_bse" plan
    ! Set up and solve BSE with ISDF compression of the BSH and Lanczos solver
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="fastBSE_isdf_cvt"

    ! Task 501 corresponds to "isdf_lanczos_bse" plan
    ! Set up and solve BSE with ISDF compression of the BSH and Lanczos solver
    nxstasks = nxstasks+1
    input%xs%plan%doonlyarray(nxstasks)%doonly%task="fastBSE_main"

    ! Set associated taks numbers
    do i = 1, nxstasks
     input%xs%plan%doonlyarray(i)%doonly%tasknumber = stringtonumberdoonlytask(input%xs%plan%doonlyarray(i)%doonly%task)
    end do
    
  end subroutine xs_plan_fastBSE


  

end module xs_plan_manager