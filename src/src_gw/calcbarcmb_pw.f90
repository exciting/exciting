
subroutine calcbarcmb_pw(iq)
      
    use modmain
    use modgw

    implicit none
    integer(4), intent(in) :: iq ! index of the q-point

    integer(4) :: i, ipw, ipw0, npw
    real(8) :: vc
    real(8) :: gpq(3), gpq2
    complex(8), allocatable :: tmat(:,:)
    
    real(8) :: kxy, kz, ab_plane, ab_norm(3), q0_vol
    
    npw = Gqbarc%ngk(1,iq)
    allocate(tmat(matsiz,npw))
    tmat(:,:) = zzero
    
    ipw0 = 1
    if (Gamma) then
      ipw0 = 2
      if (vccut) tmat(:,1) = i_sz*mpwmix(:,1)
    end if ! Gamma
    
    !------------------
    ! Loop over G+q
    !------------------
    do ipw = ipw0, npw
    
      gpq(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      gpq2 = gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3)
      
      if (vccut) then
      
        select case (trim(input%gw%barecoul%cutofftype))
        
          case('0d')
            vc = 4.0d0*pi/gpq2
            vc = vc*(1.d0-dcos(dsqrt(gpq2)*rccut))
      
          case('2d')
            ! version by Ismail-Beigi (fixed rc = L_z/2)
            kxy = dsqrt(gpq(1)*gpq(1)+gpq(2)*gpq(2))
            kz = dabs(gpq(3))
            vc = 4.0d0*pi/gpq2
            vc = vc*(1.d0-dexp(-kxy*rccut)*dcos(kz*rccut))
            
          case default
            write(*,*) 'ERROR(calcbarcmb_pw): Specified cutoff type is not implemented!'
          
        end select
        
      else
        ! no cutoff
        vc = 4.0d0*pi/gpq2
        
      end if
      
      tmat(:,ipw) = vc*mpwmix(:,ipw)
      
    end do ! ipw
    
    call zgemm('n','c',matsiz,matsiz,npw, &
    &          zone, &
    &          tmat,matsiz, &
    &          mpwmix,matsiz, &
    &          zzero, &
    &          barc,matsiz)
    
    deallocate(tmat)
      
    return  
end subroutine
