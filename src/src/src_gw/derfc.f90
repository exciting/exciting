!BOP
! !ROUTINE: derfc
! 
! !INTERFACE:
      real(8) function derfc(x)
!
! !INPUT PARAMETERS:      
!    
!  real(8) :: x
!
! !DESCRIPTION:
!
! Error function in double precision. Some compilers include it as an 
! internal procedure, but since some of them don't, it is better to have
! the external one
!
!
! !LOCAL VARIABLES:

      implicit none
      
      real(8) :: t
      real(8) :: u
      real(8) :: x
      real(8) :: y
      

! !DEFINED PARAMETERS:
!
      real(8), parameter :: &
     &    pa = 3.97886080735226000d+00, &
     &    p0 = 2.75374741597376782d-01, &
     &    p1 = 4.90165080585318424d-01, &
     &    p2 = 7.74368199119538609d-01, &
     &    p3 = 1.07925515155856677d+00, &
     &    p4 = 1.31314653831023098d+00, &
     &    p5 = 1.37040217682338167d+00, &
     &    p6 = 1.18902982909273333d+00, &
     &    p7 = 8.05276408752910567d-01, &
     &    p8 = 3.57524274449531043d-01, &
     &    p9 = 1.66207924969367356d-02, &
     &    p10 = -1.19463959964325415d-01, &
     &    p11 = -8.38864557023001992d-02
      real(8), parameter :: &
     &    p12 = 2.49367200053503304d-03, &
     &    p13 = 3.90976845588484035d-02, &
     &    p14 = 1.61315329733252248d-02, &
     &    p15 = -1.33823644533460069d-02, &
     &    p16 = -1.27223813782122755d-02, &
     &    p17 = 3.83335126264887303d-03, &
     &    p18 = 7.73672528313526668d-03, &
     &    p19 = -8.70779635317295828d-04, &
     &    p20 = -3.96385097360513500d-03, &
     &    p21 = 1.19314022838340944d-04, &
     &    p22 = 1.27109764952614092d-03
!EOP
!
!BOC
      t = pa / (pa + abs(x))
      u = t - 0.5d0
      y = (((((((((p22 * u + p21) * u + p20) * u + &
     &    p19) * u + p18) * u + p17) * u + p16) * u + &
     &    p15) * u + p14) * u + p13) * u + p12
      y = ((((((((((((y * u + p11) * u + p10) * u + &
     &    p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + &
     &    p4) * u + p3) * u + p2) * u + p1) * u + p0) * t * &
     &    exp(-x * x)
      if (x .lt. 0) y = 2 - y
      derfc = y
      end function derfc
!EOC
