
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
real(8) f1(24),f2(38)
data f1 / &
                       1.d0,                        2.d0,  &
                       6.d0,                       24.d0,  &
                     120.d0,                      720.d0,  &
                    5040.d0,                    40320.d0,  &
                  362880.d0,                  3628800.d0,  &
                39916800.d0,                479001600.d0,  &
              6227020800.d0,              87178291200.d0,  &
           1307674368000.d0,           20922789888000.d0,  &
         355687428096000.d0,         6402373705728000.d0,  &
      121645100408832000.d0,      2432902008176640000.d0,  &
    51090942171709440000.d0,   1124000727777607680000.d0,  &
 25852016738884976640000.d0, 620448401733239439360000.d0 /
data f2 / &
                       1.d0,                        2.d0,  &
                       3.d0,                        8.d0,  &
                      15.d0,                       48.d0,  &
                     105.d0,                      384.d0,  &
                     945.d0,                     3840.d0,  &
                   10395.d0,                    46080.d0,  &
                  135135.d0,                   645120.d0,  &
                 2027025.d0,                 10321920.d0,  &
                34459425.d0,                185794560.d0,  &
               654729075.d0,               3715891200.d0,  &
             13749310575.d0,              81749606400.d0,  &
            316234143225.d0,            1961990553600.d0,  &
           7905853580625.d0,           51011754393600.d0,  &
         213458046676875.d0,         1428329123020800.d0,  &
        6190283353629375.d0,        42849873690624000.d0,  &
      191898783962510625.d0,      1371195958099968000.d0,  &
     6332659870762850625.d0,     46620662575398912000.d0,  &
   221643095476699771875.d0,   1678343852714360832000.d0,  &
  8200794532637891559375.d0,  63777066403145711616000.d0 /
! fast return if possible
if (n.eq.0) then
  factnm=1.d0
  return
end if
if (m.eq.1) then
  if ((n.ge.1).and.(n.le.24)) then
    factnm=f1(n)
    return
  end if
end if
if (m.eq.2) then
  if ((n.ge.1).and.(n.le.38)) then
    factnm=f2(n)
    return
  end if
end if
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
if (m.eq.1) then
  factnm=f1(24)
  do i=25,n
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

