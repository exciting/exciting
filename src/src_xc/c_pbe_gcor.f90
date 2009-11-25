!
!
!
! This routine is based on code written by K. Burke.
!
!
Subroutine c_pbe_gcor (a, a1, b1, b2, b3, b4, rtrs, gg, ggrs)
      Implicit None
! arguments
      Real (8), Intent (In) :: a
      Real (8), Intent (In) :: a1
      Real (8), Intent (In) :: b1
      Real (8), Intent (In) :: b2
      Real (8), Intent (In) :: b3
      Real (8), Intent (In) :: b4
      Real (8), Intent (In) :: rtrs
      Real (8), Intent (Out) :: gg
      Real (8), Intent (Out) :: ggrs
! local variables
      Real (8) :: q0, q1, q2, q3
      q0 = - 2.d0 * a * (1.d0+a1*rtrs*rtrs)
      q1 = 2.d0 * a * rtrs * (b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = Log (1.d0+1.d0/q1)
      gg = q0 * q2
      q3 = a * (b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
      ggrs = - 2.d0 * a * a1 * q2 - q0 * q3 / (q1*(1.d0+q1))
      Return
End Subroutine
