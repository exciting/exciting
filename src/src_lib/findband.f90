
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findband
! !INTERFACE:

subroutine findband(findlinentype,l, k, np, nr, r, vr, de0, e)
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
  character(*), intent(in) :: findlinentype
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
  character*4, parameter :: emain='continue'
  integer::ie, nn
  ! energy search tolerance
  real(8), parameter :: eps=1.d-5
  real(8)::de, et, eb, t, tp
  real(8)::u, uup, udn, upold, ddnold, dnold, du, dudn, dupold, duup, utst, dutst
  real(8)::e1, e2, eidn, eiup
  ! automatic arrays
  real(8)::p0(nr), p1(nr), q0(nr), q1(nr)
  select case(trim(findlinentype))
  case('simple')
     tp=0.d0
     ! find the top of the band
     de=abs(de0)
     et=e
     do ie=1, maxstp
        et=et+de
        call rschroddme(0, l, k, et, np, nr, r, vr, nn, p0, p1, q0, q1)
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
10   continue
     ! find the bottom of the band
     de=-abs(de0)
     eb=et+5.d0*abs(de0)
     do ie=1, maxstp
        eb=eb+de
        call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
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
20   continue
     ! set the band energy to the mid-point
     e=(et+eb)/2.d0
30   continue
  case('advanced')
     de=de0
     call rschroddme(0,l,k,e,np,nr,r,vr,nn,p0,p1,q0,q1)
     u=p0(nr)
     du=p1(nr)
     eiup=e
     eidn=e
     e1=-1000.0d0
     e2=-1000.0d0
     upold=u
     dupold=du
     dnold=u
     ddnold=du
101  continue
     eiup=eiup+de
     call rschroddme(0,l,k,eiup,np,nr,r,vr,nn,p0,p1,q0,q1)
     uup=p0(nr)
     duup=p1(nr)
     utst=upold*uup
     dutst=dupold*duup
     upold=uup
     dupold=duup
     if (utst.lt.0.d0) then
        e2=eiup
        if (e1.gt.-30.d0) goto 301
     end if
     if (dutst.lt.0.d0) then
        e1=eiup
        if (e2.gt.-30.d0) goto 301
     end if
     if ((e1.lt.-30.d0).and.(eiup.lt.2.5d0)) then
201     continue
        eidn=eidn-de
        call rschroddme(0,l,k,eidn,np,nr,r,vr,nn,p0,p1,q0,q1)
        udn=p0(nr)
        dudn=p1(nr)
        utst=dnold*udn
        dutst=ddnold*dudn
        dnold=udn
        ddnold=dudn
        if (utst.lt.0.d0) then
           e2=eidn
           if (e1 .gt. -30.d0) goto 301
        end if
        if (dutst .lt. 0.d0) then
           e1=eidn
           if (e2 .gt. -30.d0) goto 301
        end if
        if (e2 .lt. -30.d0) then
           goto 101
        else if (eidn .gt. (e-3.d0)) then
           goto 201
        end if
     else if (eiup .lt. 0.5d0) then
        goto 101
     end if
301  continue
     if ((e1 .lt. -30.d0) .and. (e2 .lt. -30.d0)) then
        if (emain .eq. 'stop') then
           goto 900
        else
           e=1.d0
           write(*,*) 'l,ei,e1,e2', l, e, e1, e2
        end if
     else if (e2 .lt. -30.d0) then
        if (emain .eq. 'stop') then
           goto 900
        else
           e=max(e1,e)
           write(*,*) 'l,ei,e1,e2', l, e, e1, e2
        end if
     else if (e1 .lt. -30.d0) then
        goto 900
     else
        e=(e1+e2)*0.5d0
        write(*,*) 'l,ei,e1,e2', l, e, e1, e2
     end if
     return
900  continue
     write(*,*)
     write(*,'("Error(findband): no energy limits found for l=",i2)') l
     write(*,'("E-bottom ",g18.10,4x,"E-top ",g18.10)') e1,e2
     stop
  case default
     write(*,*)
     write(*,'("Error(findband): No such method for search of linearization energies: ",a)') &
          trim(findlinentype)
     write(*,*)
     stop
  end select
end subroutine findband
!EOC
