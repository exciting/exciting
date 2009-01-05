
! This routine is based on code written by K. Burke.

subroutine c_pbe_gcor(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
implicit none
! arguments
real(8), intent(in) :: a
real(8), intent(in) :: a1
real(8), intent(in) :: b1
real(8), intent(in) :: b2
real(8), intent(in) :: b3
real(8), intent(in) :: b4
real(8), intent(in) :: rtrs
real(8), intent(out) :: gg
real(8), intent(out) :: ggrs
! local variables
real(8) q0,q1,q2,q3
q0=-2.d0*a*(1.d0+a1*rtrs*rtrs)
q1=2.d0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
q2=log(1.d0+1.d0/q1)
gg=q0*q2
q3=a*(b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
ggrs=-2.d0*a*a1*q2-q0*q3/(q1*(1.d0+q1))
return
end subroutine

