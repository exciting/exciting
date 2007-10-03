
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findband
! !INTERFACE:
subroutine findband(l,k,np,nr,r,vr,de0,e)
! !INPUT/OUTPUT PARAMETERS:
!   l   : angular momentum quantum number (in,integer)
!   k   : quantum number k, zero if Dirac eqn. is not to be used (in,integer)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   de0 : default energy step size (in,real)
!   e   : input energy and returned band energy (inout,real)
! !DESCRIPTION:
!   Finds the band energies for a given radial potential and angular momentum.
!   This is done by first searching upwards in energy until the radial
!   wavefunction at the muffin-tin radius is zero. This is the energy at the top
!   of the band, denoted $E_{\rm t}$. A downward search is now performed from
!   $E_{\rm t}$ until the slope of the radial wavefunction at the muffin-tin
!   radius is zero. This energy, $E_{\rm b}$, is at the bottom of the band. The
!   band energy is taken as $(E_{\rm t}+E_{\rm b})/2$. If either $E_{\rm t}$ or
!   $E_{\rm b}$ cannot be found then the band energy is set to the input value.
!
! !REVISION HISTORY:
!   Created September 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: l
integer, intent(in) :: k
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
real(8), intent(in) :: de0
real(8), intent(inout) :: e
! local variables
! maximum number of steps
integer, parameter :: maxstp=1000
integer ie,nn
! energy search tolerance
real(8), parameter :: eps=1.d-5
real(8) de,et,eb,t,tp
! automatic arrays
real(8) p0(nr),p1(nr),q0(nr),q1(nr)
tp=0.d0
! find the top of the band
de=abs(de0)
et=e
do ie=1,maxstp
  et=et+de
  call rschroddme(0,l,k,et,np,nr,r,vr,nn,p0,p1,q0,q1)
  t=p0(nr)
  if (ie.gt.1) then
    if (t*tp.le.0.d0) then
      if (abs(de).lt.eps) goto 10
      de=-0.5d0*de
    end if
  end if
  tp=t
end do
goto 30
10 continue
! find the bottom of the band
de=-abs(de0)
eb=et+5.d0*abs(de0)
do ie=1,maxstp
  eb=eb+de
  call rschroddme(0,l,k,eb,np,nr,r,vr,nn,p0,p1,q0,q1)
  t=p1(nr)
  if (ie.gt.1) then
    if (t*tp.le.0.d0) then
      if (abs(de).lt.eps) goto 20
      de=-0.5d0*de
    end if
  end if
  tp=t
end do
goto 30
20 continue
! set the band energy to the mid-point
e=(et+eb)/2.d0
30 continue
return
end subroutine
!EOC

