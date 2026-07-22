!> Manage allocation, dispatch, and cleanup for the ground-state mixer state.
Module modmixer_lifecycle
  Use mod_LDA_LU, Only: ldapu
  Use modinput, Only: input
  Use modmain, Only: ngrtot, natmtot, nrmtmax, nspecies, rmt, rhomt, rhoir, veffmt, veffir
  Use modmpi, Only: terminate_if_false
  Use modmixer, Only: Amix, beta_ir, beta_mt, densitymixing, history_ir, history_mt, history_size, &
 &    lmmaxmix, lmaxmix, npsd, pickmixer, residual_history_ir, residual_history_mt
  implicit none

Contains

  !> Initialize mixer state, history storage, and compatibility checks for the selected mixer.
  Subroutine initmixer
    implicit none

    densitymixing = (input%groundstate%mixerswitch .eq. 2)
    lmaxmix = input%groundstate%lmaxvr
    lmmaxmix = (lmaxmix + 1)**2
    npsd = nint(0.25d0 * input%groundstate%gmaxvr * maxval(rmt(1:nspecies)))

    call terminate_if_false((input%groundstate%mixerswitch .eq. 1) .or. &
 &   (input%groundstate%mixerswitch .eq. 2), &
 &   'Error(initmixer): mixerswitch must be 1 (potential) or 2 (density)')

    call terminate_if_false(.not.associated(input%groundstate%mgga), &
 &   'Error(initmixer): mixer=' // trim(input%groundstate%mixer) // &
 &   ' is not supported with meta-GGA calculations')
    call terminate_if_false(.not.associated(input%groundstate%spin), &
 &   'Error(initmixer): mixer=' // trim(input%groundstate%mixer) // &
 &   ' is not supported for spin-polarized calculations')
    call terminate_if_false(ldapu .eq. 0, &
 &   'Error(initmixer): mixer=' // trim(input%groundstate%mixer) // &
 &   ' is not supported with LDA+U calculations')
    if (input%groundstate%mixer .eq. 'kerker') then
      call terminate_if_false(input%groundstate%mixerswitch .eq. 2, &
 &     'Error(initmixer): mixer=kerker requires mixerswitch=2 for density mixing')
    end if
    if ((input%groundstate%mixer .eq. 'kerker') .or. &
 &      ((input%groundstate%mixer .eq. 'pulay') .and. densitymixing)) then
      call terminate_if_false(input%groundstate%lambda .gt. 0.d0, &
 &     'Error(initmixer): lambda must be positive for screened density mixing')
    end if

    if (input%groundstate%mixer .eq. 'lin') then
      history_size = 1
      allocate(beta_ir(ngrtot))
      allocate(beta_mt(lmmaxmix, nrmtmax, natmtot))
      beta_ir = input%groundstate%beta0
      beta_mt = input%groundstate%beta0
    elseif (input%groundstate%mixer .eq. 'msec') then
      history_size = input%groundstate%msecStoredSteps
    elseif (input%groundstate%mixer .eq. 'pulay') then
      history_size = input%groundstate%pulayStoredSteps
      call terminate_if_false(history_size .gt. 0, &
 &     'Error(initmixer): pulayStoredSteps must be positive')
      allocate(Amix(history_size, history_size))
      Amix = 0d0
    elseif ((input%groundstate%mixer .eq. 'kerker') .or. &
 &          (input%groundstate%mixer .eq. 'simplelinear')) then
      history_size = 1
    else
      call terminate_if_false(.false., &
 &     'Error(initmixer): unknown mixer ' // trim(input%groundstate%mixer))
    end if

    allocate(history_ir(ngrtot, history_size))
    allocate(residual_history_ir(ngrtot, history_size))
    allocate(history_mt(lmmaxmix, nrmtmax, natmtot, history_size))
    allocate(residual_history_mt(lmmaxmix, nrmtmax, natmtot, history_size))

    residual_history_ir = 0d0
    residual_history_mt = 0d0

    if (densitymixing) then
      history_ir(:, 1) = rhoir
      history_mt(:, :, :, 1) = rhomt
    else
      history_ir(:, 1) = veffir
      history_mt(:, :, :, 1) = veffmt
    end if
  End Subroutine

  !> Run the active density or potential mixer on SCF iteration `sclstep`.
  Subroutine runmixer(sclstep)
    Use modmain, Only: rhomt, rhoir, veffmt, veffir
    implicit none
    !> Current self-consistent-field iteration index.
    integer, intent(in) :: sclstep

    if (densitymixing) then
      call pickmixer(rhoir, rhomt, sclstep)
    else
      call pickmixer(veffir, veffmt, sclstep)
    end if
  End Subroutine

  !> Release all mixer-owned arrays allocated during `initmixer`.
  Subroutine finish_mixer
    implicit none

    if (allocated(Amix)) deallocate(Amix)
    if (allocated(history_ir)) deallocate(history_ir)
    if (allocated(residual_history_ir)) deallocate(residual_history_ir)
    if (allocated(beta_ir)) deallocate(beta_ir)
    if (allocated(history_mt)) deallocate(history_mt)
    if (allocated(residual_history_mt)) deallocate(residual_history_mt)
    if (allocated(beta_mt)) deallocate(beta_mt)
  End Subroutine

End Module
