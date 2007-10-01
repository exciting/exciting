
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sdelta
! !INTERFACE:
real(8) function sdelta(stype,x)
! !INPUT/OUTPUT PARAMETERS:
!   stype : smearing type (in,integer)
!   x     : real argument (in,real)
! !DESCRIPTION:
!   Returns a normalised smooth approximation to the Dirac delta function. These
!   functions are defined such that
!   $$ \int\tilde{\delta}(x)dx=1. $$
!   The effective width, $w$, of the delta function may be varied by using the
!   normalising transformation
!   $$ \tilde{\delta}_w(x)\equiv\frac{\tilde{\delta}(x/w)}{w}. $$
!   Currently implimented are:
!   \begin{list}{}{\itemsep -2pt}
!    \item[0.] Gaussian
!    \item[1.] Methfessel-Paxton order 1
!    \item[2.] Methfessel-Paxton order 2
!    \item[3.] Fermi-Dirac
!   \end{list}
!   See routines {\tt stheta}, {\tt sdelta\_mp} and {\tt sdelta\_fd}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: stype
real(8), intent(in) :: x
! external functions
real(8) sdelta_mp,sdelta_fd
external sdelta_mp,sdelta_fd
sdelta=0.d0
select case(stype)
case(0)
  sdelta=sdelta_mp(0,x)
  return
case(1)
  sdelta=sdelta_mp(1,x)
  return
case(2)
  sdelta=sdelta_mp(2,x)
  return
case(3)
  sdelta=sdelta_fd(x)
  return
case default
  write(*,*)
  write(*,'("Error(sdelta): sytpe not defined : ",I8)') stype
  write(*,*)
  stop
end select
end function
!EOC

!BOP
! !ROUTINE: getsdata
! !INTERFACE:
subroutine getsdata(stype,sdescr)
! !INPUT/OUTPUT PARAMETERS:
!   stype  : smearing type (in,integer)
!   sdescr : smearing scheme description (out,character(256))
! !DESCRIPTION:
!   Returns a description of the smearing scheme as string {\tt sdescr} up
!   to 256 characters long.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: stype
character(256), intent(out) :: sdescr
select case(stype)
case(0)
  sdescr='Gaussian'
  return
case(1)
  sdescr='Methfessel-Paxton order 1, Phys. Rev. B 40, 3616 (1989)'
  return
case(2)
  sdescr='Methfessel-Paxton order 2, Phys. Rev. B 40, 3616 (1989)'
  return
case(3)
  sdescr='Fermi-Dirac'
  return
case default
  write(*,*)
  write(*,'("Error(getsdata): sytpe not defined : ",I8)') stype
  write(*,*)
  stop
end select
end subroutine
!EOC
