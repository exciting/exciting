
subroutine calcbarcmb_pw(iq)
      
    use modmain
    use modgw
    use mod_coulomb_potential

    implicit none
    integer(4), intent(in) :: iq ! index of the q-point

    integer(4) :: i, ipw, ipw0, npw
    real(8) :: gpq(3), gpq2
    real(8),    allocatable :: vc(:)
    complex(8), allocatable :: tmat(:,:)
    
    npw = Gqbarc%ngk(1,iq)
    allocate(vc(npw))
       
    select case (trim(input%gw%barecoul%cutofftype))

      case('0d')
        call vcoul_0d(Gamma, iq, Gqbarc, vc)

      case('1d') 
        if (trim(input%gw%selfenergy%singularity) == 'rim') then
            call vcoul_1d_RIM(Gamma, input%gw%ngridq, iq, Gqbarc, vc)
            ! if (Gamma) singc2 = vc(1)/(4.d0*pi)/dble(kqset%nkpt)
        else
            call vcoul_1d(Gamma, iq, Gqbarc, vc)
            ! call vcoul_1d_Rozzi(Gamma, iq, Gqbarc, vc)
            ! if (Gamma) singc2 = vc(1)/(4.d0*pi)/dble(kqset%nkpt)
        end if

      case('2d')
        call vcoul_2d(Gamma, iq, Gqbarc, vc)

      case('none')       
        if (trim(input%gw%selfenergy%singularity) == 'rim') then
            call vcoul_3d_RIM(Gamma, input%gw%ngridq, iq, Gqbarc, vc)
            if (Gamma) singc2 = vc(1)/(4.d0*pi)/dble(kqset%nkpt)
        else
            call vcoul_3d(Gamma, iq, Gqbarc, vc)
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
