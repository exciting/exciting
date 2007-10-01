
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: clebgor
! !INTERFACE:
real(8) function clebgor(j1,j2,j3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Clebsch-Gordon coefficients using the Wigner $3j$-symbols
!   $$ C(J_1 J_2 J_3 | m_1 m_2 m_3)=(-1)^{J_1-J_2+m_3}\sqrt{2J_3+1}
!    \begin{pmatrix} J_1 & J_2 & J_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}. $$
!   Suitable for $J_i\le 50$.
!
! !REVISION HISTORY:
!   Created September 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j1
integer, intent(in) :: j2
integer, intent(in) :: j3
integer, intent(in) :: m1
integer, intent(in) :: m2
integer, intent(in) :: m3
! external functions
real(8) wigner3j
external wigner3j
if ((j1.lt.0).or.(j2.lt.0).or.(j3.lt.0).or.(abs(m1).gt.j1).or.(abs(m2).gt.j2) &
 .or.(abs(m2).gt.j2).or.(abs(m3).gt.j3)) then
  write(*,*)
  write(*,'("Error(clebgor): non-physical arguments :")')
  write(*,'("j1 = ",I8," j2 = ",I8," j3 = ",I8)') j1,j2,j3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((j1.eq.0).and.(j2.eq.0).and.(j3.eq.0)) then
  clebgor=1.d0
  return
end if
if ((j1.gt.50).or.(j2.gt.50).or.(j3.gt.50)) then
  write(*,*)
  write(*,'("Error(clebgor): angular momenta out of range : ",3I8)') j1,j2,j3
  write(*,*)
  stop
end if
if ((m1+m2-m3.ne.0).or.(j2+j3.lt.j1).or.(j1+j3.lt.j2).or.(j1+j2.lt.j3)) then
  clebgor=0.d0
  return
end if
clebgor=sqrt(dble(2*j3+1))*wigner3j(j1,j2,j3,m1,m2,-m3)
if (mod(j1-j2+m3,2).ne.0) clebgor=-clebgor
return
end function
!EOC
