!
subroutine calcbarcmb_ipw_ipw(iq)
!
    use modinput
    use modmain,               only : zzero, zone, pi
    use modgw,                 only : Gset, Gamma, kqset, Gqset, Gqbarc
    use mod_product_basis,     only : locmatsiz, mpwipw
    use mod_coulomb_potential, only : barc, rccut, i_sz, vccut
    implicit none
    ! input variables
    integer, intent(in) :: iq
    ! local variables 
    integer :: npw, ipw, ipw0  ! PW
    integer :: ngq, igq, jgq   ! IPW
    real(8) :: gqvec(3), gqlen
    real(8) :: vc, kxy, kz
    complex(8), allocatable :: tmat1(:,:), tmat2(:,:)
    ! external routine 
    external zgemm    
    
    npw = Gqbarc%ngk(1,iq)
    ngq = Gqset%ngk(1,iq)
    
    ! local array
    allocate(tmat1(1:ngq,1:npw))
    tmat1(:,:) = zzero
    
    ipw0 = 1
    if (Gamma) then 
      ipw0 = 2
      !if (vccut) tmat1(1:ngq,1) = i_sz*mpwipw(1:ngq,1)
    end if
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,gqvec,gqlen,vc,kxy,kz,igq)
!$OMP DO
#endif    
    do ipw = ipw0, npw
      gqvec(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      gqlen = gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2)+gqvec(3)*gqvec(3)
      
      if (vccut) then
      
        ! apply cutoff
        select case (trim(input%gw%barecoul%cutofftype))
      
          case('0d')
            vc = 1.d0/gqlen
            vc = vc*(1.d0-dcos(dsqrt(gqlen)*rccut))
      
          case ('2d')
            ! version by Ismail-Beigi (fixed rc = L_z/2)
            kxy = dsqrt(gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2))
            kz = dabs(gqvec(3))
            vc = 4.0d0*pi/gqlen
            vc = vc*(1.d0-dexp(-kxy*rccut)*dcos(kz*rccut))
          
          case default
            write(*,*) 'ERROR(calcbarcmb_ipw_ipw): Specified cutoff type is not implemented!'
          
          end select
          
      else
      
        ! no cutoff
        vc = 4.0d0*pi/gqlen
        
      end if
      
      do igq = 1, ngq
        tmat1(igq,ipw) = vc*mpwipw(igq,ipw)
      end do
      
    end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
    
    allocate(tmat2(1:ngq,1:ngq))  
    call zgemm( 'n','c',ngq,ngq,npw, &
    &           zone,tmat1,ngq, &
    &           mpwipw,ngq, &
    &           zzero,tmat2,ngq)
    deallocate(tmat1)
      
    do jgq = 1, ngq
      do igq = 1, ngq
        barc(locmatsiz+igq,locmatsiz+jgq) = tmat2(igq,jgq)
      end do 
    end do
    deallocate(tmat2)
    
    !write(*,*) 'IPW-IPW'
    !do igq = 1, ngq, ngq/10
    !do jgq = 1, ngq, ngq/10
    ! write(*,*) igq, jgq, barc(locmatsiz+igq,locmatsiz+jgq)
    !end do
    !end do
    
end subroutine
    
