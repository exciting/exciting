
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
    
    real(8) :: kxy, kz, ab_plane, ab_norm(3), q0_vol
    
    npw = Gqbarc%ngk(1,iq)
    
    allocate(vc(npw))
    vc(:) = 0.d0
        
    select case (trim(input%gw%barecoul%cutofftype))

      case('none')
        ! no cutoff
        if (Gamma) then
          ipw0 = 2
        else
          ipw0 = 1
        end if          
        do ipw = ipw0, npw
          vc(ipw) = 4.0d0*pi / Gqbarc%gkc(ipw,1,iq)**2
        end do

      case('0d')
        call vcoul_0d(iq, Gqbarc, vc)

      case('1d')
        call vcoul_1d(iq, Gqbarc, vc)
        if (Gamma) vc(1) = vcq0

      case('2d')
        call vcoul_2d(iq, Gqbarc, vc)
        if (Gamma) vc(1) = vcq0

    end select

    allocate(tmat(matsiz,npw))
    do ipw = 1, npw
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
