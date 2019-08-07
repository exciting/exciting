
subroutine calcbarcmb_pw(iq)

    use modmain
    use modgw
    use mod_coulomb_potential

    implicit none
    integer(4), intent(in) :: iq ! index of the q-point

    integer(4) :: i, ipw, ipw0, npw
    real(8) :: omega_hyb, omega2, exp_omega
    real(8) :: gpq(3), gpq2
    real(8),    allocatable :: vc(:)
    complex(8), allocatable :: tmat(:,:)

    if ((task == 7) .and. (xctype(1)==408)) then
      ! hybrid SCF
      omega_hyb = input%groundstate%Hybrid%omega
      omega2 = omega_hyb*omega_hyb
    end if

    npw = Gqbarc%ngk(1,iq)
    allocate(vc(npw))

    select case (trim(input%gw%barecoul%cutofftype))

      case('0d')
        call vcoul_0d(Gamma, iq, Gqbarc, vc)

      case('1d')
        call vcoul_1d(Gamma, iq, Gqbarc, vc)
        ! call vcoul_1d_Rozzi(Gamma, iq, Gqbarc, vc)
        ! if (Gamma) singc2 = vc(1)/(4.d0*pi)/dble(kqset%nkpt)

      case('2d')
        call vcoul_2d(Gamma, iq, Gqbarc, vc)

      case('none')

        if ((task == 7) .and. (xctype(1)==408)) then
          ! hybrid SCF
          exp_omega = exp(-gpq2/(4.d0*omega2))
          vc = ((4.0d0*pi)/gpq2)*(1.d0-exp_omega)
        else
          ! GW
          if (trim(input%gw%selfenergy%singularity) == 'rim') then
              call vcoul_3d_RIM(Gamma, input%gw%ngridq, iq, Gqbarc, vc)
              if (Gamma) singc2 = vc(1)/(4.d0*pi)/dble(kqset%nkpt)
          else
              call vcoul_3d(Gamma, iq, Gqbarc, vc)
          end if
       end if

    end select

    if (Gamma) then
        ipw0 = 2
    else
        ipw0 = 1
    end if

    allocate(tmat(matsiz,npw))
    do ipw = ipw0, npw
      tmat(:,ipw) = vc(ipw)*mpwmix(:,ipw)
    end do ! ipw
    deallocate(vc)

    call zgemm('n', 'c', matsiz, matsiz, npw, &
               zone, &
               tmat, matsiz, &
               mpwmix, matsiz, &
               zzero, &
               barc, matsiz)
    deallocate(tmat)

    return
end subroutine
