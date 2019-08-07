subroutine calcbarcmb_lr(iq)

    use modmain
    use modgw
    use mod_hybrids, only: barc_lr
    Use modinput, only: input
    implicit none
    integer(4), intent(in) :: iq ! index of the q-point

    integer(4) :: i, ipw, ipw0, npw
    real(8) :: vc
    real(8) :: omega_hyb, omega2, exp_omega
    real(8) :: gpq(3), gpq2
    complex(8), allocatable :: tmat(:,:)

    omega_hyb = input%groundstate%Hybrid%omega
    omega2 = omega_hyb*omega_hyb

    npw = Gqbarc%ngk(1,iq)
    allocate(tmat(matsiz,npw))
    tmat(:,:) = zzero

    !------------------
    ! Loop over G+q
    !------------------
    ipw0 = 1
    if (Gamma) then
      ipw0 = 2
      ! tmat(:,1) = -(pi/omega2)*mpwmix(:,1) !!(?)
    end if

    do ipw = ipw0, npw
      gpq(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      gpq2 = gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3)
      ! is missing 1/2, NO LOOK at zgemm
      exp_omega = exp(-gpq2/(4.0d0*omega2))
      vc = (4.0d0*pi/gpq2)*exp_omega
      tmat(:,ipw) = vc*mpwmix(:,ipw)
    end do ! ipw

    call zgemm('n', 'c', matsiz, matsiz, npw, &
    &          zone, &
    &          tmat, matsiz, &
    &          mpwmix, matsiz, &
    &          zzero, &
    &          barc_lr, matsiz)

    deallocate(tmat)

    return
end subroutine
