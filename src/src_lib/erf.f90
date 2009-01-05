
! This code is based on a routine from the NSWC Library of Mathematics
! Subroutines and is in the public domain.

!BOP
! !ROUTINE: erf
! !INTERFACE:
real(8) function erf(x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the error function ${\rm erf}(x)$ using a rational function
!   approximation. This procedure is numerically stable and accurate to near
!   machine precision.
!
! !REVISION HISTORY:
!   Modified version of a NSWC routine, April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x
! local variables
real(8) ax,x2,t,bot,top
real(8), parameter :: c=0.564189583547756d0
real(8), parameter :: &
 a(5)=(/ 0.771058495001320d-04,-0.133733772997339d-02, &
         0.323076579225834d-01, 0.479137145607681d-01, &
         0.128379167095513d+00 /), &
 b(3)=(/ 0.301048631703895d-02, 0.538971687740286d-01, &
         0.375795757275549d+00 /), &
 p(8)=(/-1.36864857382717d-07,  5.64195517478974d-01, &
         7.21175825088309d+00,  4.31622272220567d+01, &
         1.52989285046940d+02,  3.39320816734344d+02, &
         4.51918953711873d+02,  3.00459261020162d+02 /), &
 q(8)=(/ 1.00000000000000d+00,  1.27827273196294d+01, &
         7.70001529352295d+01,  2.77585444743988d+02, &
         6.38980264465631d+02,  9.31354094850610d+02, &
         7.90950925327898d+02,  3.00459260956983d+02 /), &
 r(5)=(/ 2.10144126479064d+00,  2.62370141675169d+01, &
         2.13688200555087d+01,  4.65807828718470d+00, &
         2.82094791773523d-01 /), &
 s(4)=(/ 9.41537750555460d+01,  1.87114811799590d+02, &
         9.90191814623914d+01,  1.80124575948747d+01 /)
ax=abs(x)
if (ax.lt.0.5d0) then
  t=x**2
  top=((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5))+1.d0
  bot=((b(1)*t+b(2))*t+b(3))*t+1.d0
  erf=x*(top/bot)
  return
end if
if (ax.lt.4.d0) then
  top=((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+p(6))*ax+p(7))*ax+p(8)
  bot=((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+q(6))*ax+q(7))*ax+q(8)
  erf=0.5d0+(0.5d0-exp(-x**2)*top/bot)
  if (x.lt.0.d0) erf=-erf
  return
end if
if (ax.lt.5.8d0) then
  x2=x**2
  t=1.d0/x2
  top=(((r(1)*t+r(2))*t+r(3))*t+r(4))*t+r(5)
  bot=(((s(1)*t+s(2))*t+s(3))*t+s(4))*t+1.d0
  erf=(c-top/(x2*bot))/ax
  erf=0.5d0+(0.5d0-exp(-x2)*erf)
  if (x.lt.0.d0) erf=-erf
  return
end if
erf=sign(1.d0,x)
return
end function
!EOC
