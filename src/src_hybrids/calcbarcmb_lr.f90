subroutine calcbarcmb_lr(iq)
      
    use modmain
    use modgw
    use mod_hybrids, only: barc_lr
    Use modinput, only: input
    implicit none
    integer(4), intent(in) :: iq ! index of the q-point

    integer(4) :: i, ipw, ipw0, npw
    real(8) :: vc, exp_omega
    real(8) :: omega_hyb
    real(8) :: gpq(3), gpq2
    complex(8), allocatable :: tmat(:,:)
       
    real(8) :: kxy, kz, ab_plane, ab_norm(3), q0_vol
    
    omega_hyb=input%groundstate%Hybrid%omega
     
    npw = Gqbarc%ngk(1,iq)
    allocate(tmat(matsiz,npw))
    tmat(:,:) = zzero
    ipw0 = 1
    if (Gamma) then
      ipw0 = 2
      tmat(:,1) = (pi/omega_hyb)*mpwmix(:,1) !!(?)
    end if ! Gamma
    
    !------------------
    ! Loop over G+q
    !------------------
    do ipw = ipw0, npw
    
      gpq(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      gpq2 = gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3)
!      write(*,*) "deb4", gpq2     
      
        ! no cutoff
        !CECI this is the 3d
        !is missing 1/2, NO LOOK at zgemm
        exp_omega=exp(-gpq2/(4*omega_hyb**2)) 
        vc = (4.0d0*pi/gpq2)*exp_omega
!    write(*,*) "deb5"     
        
      
      tmat(:,ipw) = vc*mpwmix(:,ipw)
      
    end do ! ipw
!    write(*,*) "deb6", tmat(:,1)  
    
    call zgemm('n','c',matsiz,matsiz,npw, &
    &          zone, &
    &          tmat,matsiz, &
    &          mpwmix,matsiz, &
    &          zzero, &
    &          barc_lr,matsiz)
    
    deallocate(tmat)
      
    return  
end subroutine

