
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: factnm
! !INTERFACE:
real(8) function factnm(n,m)
! !INPUT/OUTPUT PARAMETERS:
!   n : input (in,integer)
!   m : order of multifactorial (in,integer)
! !DESCRIPTION:
!   Returns the multifactorial
!   $$ n\underbrace{!!\,...\,!}_{m\,{\rm times}}=\prod_{i\ge 0,\,n-im>0}n-im $$
!   for $n,\,m \ge 0$. $n$ should be less than 150.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
integer, intent(in) :: m
! local variables
integer i,j
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(factnm): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
if (m.le.0) then
  write(*,*)
  write(*,'("Error(factnm): m <= 0 : ",I8)') m
  write(*,*)
  stop
end if
if (n.gt.150) then
  write(*,*)
  write(*,'("Error(factnm): n out of range : ",I8)') n
  write(*,*)
  stop
end if
if (n.eq.0) then
  factnm=1.d0
  return
end if
if (m.eq.1) then
  factnm=1.d0
  do i=2,n
    factnm=factnm*dble(i)
  end do
else
  j=n/m
  if (mod(n,m).eq.0) j=j-1
  factnm=dble(n)
  do i=1,j
    factnm=factnm*dble(n-i*m)
  end do
end if
return
end function
!EOC
