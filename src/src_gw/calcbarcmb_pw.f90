
subroutine calcbarcmb_pw(iq)
      
    use modmain
    use modgw
    Use modinput, only: input

    implicit none
    integer(4), intent(in) :: iq ! index of the q-point

    integer(4) :: i, ipw, ipw0, npw
    real(8) :: vc
    real(8) :: omega_hyb, omega2, exp_omega
    real(8) :: gpq(3), gpq2
    complex(8), allocatable :: tmat(:,:)
    
    real(8) :: kxy, kz, ab_plane, ab_norm(3), q0_vol
    
    npw = Gqbarc%ngk(1,iq)
    allocate(tmat(matsiz,npw))
    tmat(:,:) = zzero
    
    if (xctype(1)==408) then
       omega_hyb=input%groundstate%Hybrid%omega
       omega2=omega_hyb*omega_hyb
    end if
    ipw0 = 1
    if (Gamma) then
      ipw0 = 2
      if (xctype(1)==408) then
         tmat(:,1)=(pi/omega2)*mpwmix(:,1)
      else
         if (vccut) tmat(:,1) = i_sz*mpwmix(:,1)
      endif
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
        !CECI this is the 3d
        if (xctype(1)==408) then
           exp_omega=exp(-gpq2/(4.d0*omega2)) 
           vc = ((4.0d0*pi)/gpq2)*(1.d0-exp_omega)
        else
           vc = 4.0d0*pi/gpq2
        endif
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
