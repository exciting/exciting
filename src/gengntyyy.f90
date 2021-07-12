!
!BOP
! !ROUTINE: init1
! !INTERFACE:
!
!
Subroutine gengntyyy !(gntyyy)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!  Gaunt coefficients <Y3 | Y2 | Y1 >.  
!  
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3

! external functions
      real(8), external :: gaunt
!
! allocate and generate complex Gaunt coefficient array
!      If (allocated(gntyry)) deallocate (gntyry)
      Allocate (gntyyy(lmmaxapw, lmmaxapw, lmmaxapw))
      Do l1 = 0, input%groundstate%lmaxapw
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do l2 = 0, input%groundstate%lmaxapw
               Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        gntyyy (lm3, lm2, lm1) = gaunt (l1, l2, l3, m1, m2, m3)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

      Return
End Subroutine
!EOC
