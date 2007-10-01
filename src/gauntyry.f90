
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gauntyry
! !INTERFACE:
complex(8) function gauntyry(l1,l2,l3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   l1, l2, l3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the complex Gaunt-like coefficient given by
!   $\langle Y^{l_1}_{m_1}|R^{l_2}_{m_2}|Y^{l_3}_{m_3}\rangle$, where $Y_{lm}$
!   and $R_{lm}$ are the complex and real spherical harmonics, respectively.
!   Suitable for $l_i$ less than 50. See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: l1
integer, intent(in) :: l2
integer, intent(in) :: l3
integer, intent(in) :: m1
integer, intent(in) :: m2
integer, intent(in) :: m3
! local variables
! real constant sqrt(2)/2
real(8), parameter :: c1=0.7071067811865475244d0
real(8) t1,t2
! external functions
real(8) gaunt
external gaunt
if (m2.gt.0) then
  if (mod(m2,2).eq.0) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  t2=c1*(gaunt(l1,l2,l3,m1,m2,m3)+t1*gaunt(l1,l2,l3,m1,-m2,m3))
  gauntyry=cmplx(t2,0.d0,8)
else if (m2.lt.0) then
  if (mod(m2,2).eq.0) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  t2=c1*(gaunt(l1,l2,l3,m1,m2,m3)-t1*gaunt(l1,l2,l3,m1,-m2,m3))
  gauntyry=cmplx(0.d0,-t2,8)
else
  gauntyry=cmplx(gaunt(l1,l2,l3,m1,m2,m3),0.d0,8)
end if
return
end function
!EOC
