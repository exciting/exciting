! Copyright (C) 2009-2010 S. Sagmeister, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xstasklauncher
  use modinput
  use modmain, only: task
  use modxs, only: dgrid, nksubpt, iksubpt, temat, doscreen0, vkloff_xs_b, &
    & hybridhf,skipgnd, fhdf5
  use inputdom
  use mod_hdf5
  implicit none

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

  call backup0
  call backup1
  call backup2

  if(associated(input%xs%plan)) then

    call xsmain(input%xs%plan)

  else if(trim(input%xs%xstype) .eq. "TDDFT") then

    if(input%xs%tddft%do .eq. "fromkernel") go to 10

    task = 301
    call xsinit
    call xsgeneigvec
    call xsfinit

    if(input%xs%tetra%tetradf) then
#ifdef TETRA         
      task = 310
      call xsinit
      call tetcalccw
      call xsfinit
#else
      ! added by din
      write(*,*) 'tetrahedron method for xs is disabled!'
      write(*,*) 'check -dtetra option in make.inc' 
      stop
#endif            
    end if

    task = 320
    call xsinit
    call writepmatxs
    call xsfinit

    task = 330
    call xsinit
    if(temat) then
      call writeemat
    end if
    call xsfinit

    ! What is fxctypenumber ? - Benjmain Aurich 2016-08-01
    if(input%xs%tddft%fxctypenumber .eq. 7 .or. &
       & input%xs%tddft%fxctypenumber .eq. 8) then
      task = 400
      call xsinit
      call scrgeneigvec
      call xsfinit

      task = 420
      call xsinit
      call scrwritepmat
      call xsfinit

      if ((input%xs%tetra%tetradf)) then
#ifdef TETRA            
        task = 410
        call xsinit
        call scrtetcalccw
        call xsfinit
#else
        ! added by din
        write(*,*) 'tetrahedron method for xs is disabled!'
        write(*,*) 'check -dtetra option in make.inc' 
        stop
#endif
      end if

      if (input%xs%screening%do .eq. "fromscratch") then
        task = 430
        call xsinit
        call screen
        call xsfinit

        task = 440
        call xsinit
        call scrcoulint
        call xsfinit
      end if

      task = 450
      call xsinit
      call kernxc_bse
      call xsfinit

    end if

    task = 340
    call xsinit
    call df
    call xsfinit

10  continue

    task = 350
    call xsinit
    call idf
    call xsfinit

  else if(trim(input%xs%xstype)=="BSE" .and. .not. (input%xs%bse%xas .or. input%xs%bse%beyond)) then

    ! STK
    ! Apply double grid technique if requested
    if(any(input%xs%bse%ngridksub .gt. 1)) then
      dgrid = .true.
    else
      dgrid = .false.
      nksubpt = 1
    endif

    if(dgrid) then
      ! append xs output
      input%xs%tappinfo = .true.
      ! backup input xs vkloff (it will be added to all generated grids)
      vkloff_xs_b(:) = input%xs%vkloff(:)
      ! save screening status for later use
      doscreen0 = input%xs%screening%do
      ! generate subgrid
      call genksubpts
    endif

    ! Double grid related loop, if no dgrid then nksubpt=1
    do iksubpt = 1, nksubpt

      if(dgrid) call bsedgridinit

      task = 301
      call xsinit
      call xsgeneigvec
      call xsfinit

      task = 320
      call xsinit
      call writepmatxs
      call xsfinit

      task = 401
      call xsinit
      call scrgeneigvec
      call xsfinit

      if((input%xs%tetra%tetradf)) then
#ifdef TETRA            
        task = 410
        call xsinit
        call scrtetcalccw
        call xsfinit
#else
        ! added by din
        write(*,*) 'tetrahedron method for xs is disabled!'
        write(*,*) 'check -dtetra option in make.inc' 
        stop
#endif
      end if

      task = 420
      call xsinit
      call scrwritepmat
      call xsfinit

      if(input%xs%screening%do .eq. "fromscratch") then
        task = 430
        call xsinit
        call screen
        call xsfinit

        task = 440
        call xsinit
        call scrcoulint
        call xsfinit
      end if

      task = 441
      call xsinit
      call exccoulint
      call xsfinit

      task = 445
      call xsinit
      call bse
      call xsfinit

      if(dgrid) input%xs%screening%do = "skip"

    enddo

    if (dgrid) then
       call bsedgrid
       ! restore input settings
       input%xs%vkloff(:) = vkloff_xs_b(:)
       input%xs%screening%do = doscreen0
    endif

  else if(trim(input%xs%xstype)=="BSE" .and. input%xs%bse%xas .and. input%xs%bse%beyond) then

    !! Removed dubble grid code, since no-one knows how it works.
    !! Removed tetra code, since no-one knows if it works.
    
    ! HDF5 Initialzation
#ifdef _HDF5_
    call hdf5_initialize()
    fhdf5="bse_output.h5"
    call hdf5_create_file(fhdf5)
    write(*,*) 'file created'
#endif
    ! Task 301 corresponds to "xsgeneigvec" plan
    ! One shot GS calculation with xs%ngridk, xs%nempty and potential xs%vkloff.
    task = 301
    call xsinit
    call b_xsgeneigveclauncher
    call xsfinit

    ! Task 320 corresponds to "writepmatxs" plan
    ! Calculates the momentum matrix elements for the xs GS calculation.
    task = 320
    call xsinit
    call xasinit
    call writepmatxs
    call xasfinit
    call xsfinit

    ! Task 401 corresponds to "scrgeneigvec" plan
    ! One shot GS calculation with more empty states xs%screening%nempty 
    ! but otherwise identical parameters as "xsgeneigvec".
    task = 401
    call xsinit
    call b_xsgeneigveclauncher ! Does one shot GS runs with screening GS parameters
    call xsfinit

    ! Task 420 corresponds to "scrwritepmat" plan
    task = 420
    call xsinit
    call scrwritepmat ! Calls writepmatxs to write PMAT_XS.OUT
    call xsfinit

    if(input%xs%screening%do .eq. "fromscratch") then
      ! Task 430 corresponds to "screen" plan
      task = 430
      call xsinit
      call b_screenlauncher
      call xsfinit

      ! Task 440 corresponds to "scrcoulint" plan
      task = 440
      call xsinit
      call xasinit
      call b_scrcoulintlauncher
      call xasfinit
      call xsfinit
    end if

    ! Task 441 corresponds to "exccoulint" plan
    task = 441
    call xsinit
    call xasinit
    call b_exccoulintlauncher
    call xasfinit
    call xsfinit

    ! Task 445 corresponds to "bse" plan
    task = 445
    call xsinit
    call xasinit
    call b_bselauncher
    call xasfinit
    call xsfinit

    
    ! HDF5 Finalization
#ifdef _HDF5_
    call hdf5_finalize()
#endif
  else if(trim(input%xs%xstype)=="BSE" .and. input%xs%bse%beyond .eqv. .true.) then


    !! Removed dubble grid code, since no-one knows how it works.
    !! Removed tetra code, since no-one knows if it works.

    ! HDF5 Initialzation
#ifdef _HDF5_
    call hdf5_initialize()
    fhdf5="bse_output.h5"
    call hdf5_create_file(fhdf5)
#endif

    ! Task 301 corresponds to "xsgeneigvec" plan
    ! One shot GS calculation with xs%ngridk, xs%nempty and potential xs%vkloff.
    task = 301
    call xsinit
    call b_xsgeneigveclauncher
    call xsfinit

    ! Task 320 corresponds to "writepmatxs" plan
    ! Calculates the momentum matrix elements for the xs GS calculation.
    task = 320
    call xsinit
    call writepmatxs
    call xsfinit

    ! Task 401 corresponds to "scrgeneigvec" plan
    ! One shot GS calculation with more empty states xs%screening%nempty 
    ! but otherwise identical parameters as "xsgeneigvec".
    task = 401
    call xsinit
    call b_xsgeneigveclauncher ! Does one shot GS runs with screening GS parameters
    call xsfinit

    ! Task 420 corresponds to "scrwritepmat" plan
    task = 420
    call xsinit
    call scrwritepmat ! Calls writepmatxs to write PMAT_XS.OUT
    call xsfinit

    if(input%xs%screening%do .eq. "fromscratch") then
      ! Task 430 corresponds to "screen" plan
      task = 430
      call xsinit
      call b_screenlauncher
      call xsfinit

      ! Task 440 corresponds to "scrcoulint" plan
      task = 440
      call xsinit
      call b_scrcoulintlauncher
      call xsfinit
    end if

    ! Task 441 corresponds to "exccoulint" plan
    task = 441
    call xsinit
    call b_exccoulintlauncher
    call xsfinit

    ! Task 445 corresponds to "bse" plan
    task = 445
    call xsinit
    call b_bselauncher
    call xsfinit

    ! HDF5 Finalization
#ifdef _HDF5_
    call hdf5_finalize()
#endif
  else

     write (*,*) "error xstasklauncher"
     write (*,*) trim (input%xs%xstype), "no valid xstype"
     stop

  end if

end subroutine
