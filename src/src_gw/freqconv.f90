!==================================================================
! Perform frequency convolution integral 
! $ S(i\omega; \epsilon ) = 1/\pi \int_0^\infty 
!   \frac{ \epsilon - i \omega ) W(i\omega')}  
!   { (\epsilon - i \omega)^2 + \omega'^2} d\omega' $
!
! Created by Hong Jiang
!==================================================================
complex(8) function freqconv(iom,nomg,om,enk,mwm,omg,womg)

    implicit none
    ! input variables
    integer(4),intent(in) :: iom   !! control how to do convolution 
                                   !!  iom.eq.0 -- direct integraion 
                                   !!  iom.gt.0 -- indicating that om = omg(iom) 
                                   !!    so that special treatment as explained 
                                   !!    in Developer Guide (H.4) is needed 
    integer(4),intent(in) :: nomg  !! number of frequency points used for the convolution
    real(8),   intent(in) :: om    !! the frequency \omega 
    real(8),   intent(in) :: enk   !! \epsilon 
    complex(8),intent(in) :: mwm(nomg)  !! W(i\omega')
    real(8),   intent(in) :: omg(nomg)  !! the array for \omega'
    real(8),   intent(in) :: womg(nomg) !! integraion weight 
    ! local variables
    integer(4) :: jom
    complex(8) :: ffac   ! the frequecy dependent factor in each term
    complex(8) :: comsq
    complex(8) :: ommek  ! i*omega-e_nk
    complex(8) :: sc
    real(8), parameter :: pi = 3.14159265358979d0

    ommek = cmplx(enk,-1.0d0*om,8)
    sc = cmplx(0.0d0,0.0d0,8)
    if (iom>0) then 
      do jom = 1, nomg    ! integration frequencies
        comsq = cmplx(omg(jom)*omg(jom),0.0d0,8)
        ffac = cmplx(womg(jom),0.0d0,8)/(ommek*ommek+comsq)
        sc = sc+ffac*(mwm(jom)-mwm(iom))
      end do
      sc = ommek*sc/pi+mwm(iom)*sign(0.5d0,enk)
    else 
      do jom = 1, nomg   
        comsq = cmplx(omg(jom)*omg(jom),0.0d0,8)
        ffac = cmplx(womg(jom),0.0d0,8)/(ommek*ommek+comsq)
        sc = sc+ffac*mwm(jom)
      end do
      sc = ommek*sc/pi
    end if 
    freqconv = sc
  
    return
end function freqconv
