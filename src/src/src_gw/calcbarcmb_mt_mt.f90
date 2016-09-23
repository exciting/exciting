!
subroutine calcbarcmb_mt_mt(iq)
!
    use modinput
    use modmain,               only : idxas, pi
    use modgw,                 only : Gamma
    use mod_product_basis,     only : locmatsiz, mbindex, rtl, rrint
    use mod_coulomb_potential, only : barc, sgm
    use mod_misc_gw,           only : vi
    implicit none
    ! input variable
    integer, intent(in) :: iq
    ! local variables
    integer :: imix, ia, is, ias, irm, l1, m1
    integer :: jmix, ja, js, jas, jrm, l2, m2
    integer :: lm12, ijrm
    real(8) :: tg, minu
    complex(8) :: stc
    ! external function      
    real(8), external :: gettildeg

    !---------------------------------------------------------
    ! Calculate the matrix \Sigma for the structure constants
    !---------------------------------------------------------
    call sigma(iq,4*(input%gw%MixBasis%lmaxmb+1))
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(imix,is,ia,ias,irm,l1,m1,jmix,js,ja,jas,jrm,l2,m2,tg,lm12,minu,stc,ijrm)
!$OMP DO
#endif
    do imix = 1, locmatsiz
      is  = mbindex(imix,1)
      ia  = mbindex(imix,2)
      ias = idxas(ia,is)
      irm = mbindex(imix,3)
      l1  = mbindex(imix,4)
      m1  = mbindex(imix,5)
      
      do jmix = 1, locmatsiz
        js  = mbindex(jmix,1)
        ja  = mbindex(jmix,2)
        jas = idxas(ja,js)
        jrm = mbindex(jmix,3)
        l2  = mbindex(jmix,4)
        m2  = mbindex(jmix,5)
                      
        if ((.not.Gamma).or.(l1.ne.0).or.(l2.ne.0)) then

          if (jas >= ias) then
            tg = gettildeg(l1,l2,-m1,m2)
            lm12 = (l1+l2)*(l2+l1+1)+m2-m1+1
            minu = (-1.0d0)**m1
            stc = tg*sgm(ias,jas,lm12)
          else
            tg = gettildeg(l1,l2,m1,-m2)
            lm12 = (l1+l2)*(l2+l1+1)-m2+m1+1
            minu = (-1.0d0)**(m2+l1+l2)
            stc = tg*conjg(sgm(ias,jas,lm12))
          end if
          barc(imix,jmix) = rtl(irm,ias)*rtl(jrm,jas)*minu*stc
            
          ! Second term of Eq.(56)
          if ((m1==m2).and.(l1==l2).and.(ias==jas)) then
            if (jrm>=irm) then
              ijrm = irm+(jrm*(jrm-1))/2
            else
              ijrm = jrm+(irm*(irm-1))/2
            end if
            barc(imix,jmix) = barc(imix,jmix)+ &
            &                 cmplx(rrint(ijrm,ias)*4.0d0*pi/dble(2*l1+1),0.0d0,8)
          end if
          
        end if
        
      end do ! jmix
      
    end do ! imix
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

    ! deallocate (global) structure factors
    deallocate(sgm)
    
    !write(*,*) 'MT-MT'
    !do imix = 1, locmatsiz, locmatsiz/10
    !do jmix = 1, locmatsiz, locmatsiz/10
    !  write(*,*) imix, jmix, barc(imix,jmix)
    !end do
    !end do

end subroutine
