


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
integer::ie, nn
! energy search tolerance
real(8), parameter :: eps=1.d-5
real(8)::de, et, eb, t, tp
! automatic arrays
real(8)::p0(nr), p1(nr), q0(nr), q1(nr)


character*4, parameter ::        EMAIN='CONT'

real(8)   DDNOLD, DNOLD, DU, DUDN, DUPOLD, DUTST, DUUP
real(8)   E1, E2, EIDN, EIUP, U, UDN, UPOLD, UTST
real(8)   UUP

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
10 continue
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
20 continue
! set the band energy to the mid-point
e=(et+eb)/2.d0
30 continue

case('advanced')
 de=de0
      call rschroddme(0,l,k,E,np,nr,r,vr,nn,p0,p1,q0,q1)
      U=p0(nr)
      DU=p1(nr)
      EIUP = E
      EIDN = E
      E1 = -100.0D+0
      E2 = -100.0D+0
      UPOLD = U
      DUPOLD = DU
      DNOLD = U
      DDNOLD = DU
   101 CONTINUE
      EIUP = EIUP + DE
      call rschroddme(0,l,k,EIUP,np,nr,r,vr,nn,p0,p1,q0,q1)
      UUP=p0(nr)
      DUUP=p1(nr)
      UTST = UPOLD*UUP
      DUTST = DUPOLD*DUUP
      UPOLD = UUP
      DUPOLD = DUUP
      IF (UTST .LT. 0.0D+0) THEN
         E2 = EIUP
         IF (E1 .GT. -30.0D+0) GOTO 301
      ENDIF
      IF (DUTST .LT. 0.0D+0) THEN
         E1 = EIUP
         IF (E2 .GT. -30.0D+0) GOTO 301
      ENDIF
      IF ((E1 .LT. -30.0D+0) .AND. (EIUP .LT. 2.5D+0)) THEN
   201    CONTINUE
         EIDN = EIDN - DE
         call rschroddme(0,l,k,EIDN,np,nr,r,vr,nn,p0,p1,q0,q1)
         UDN=p0(nr)
         DUDN=p1(nr)
         UTST = DNOLD*UDN
         DUTST = DDNOLD*DUDN
         DNOLD = UDN
         DDNOLD = DUDN
         IF (UTST .LT. 0.0D+0) THEN
            E2 = EIDN
            IF (E1 .GT. -30.0D+0) GOTO 301
         ENDIF
         IF (DUTST .LT. 0.) THEN
            E1 = EIDN
            IF (E2 .GT. -30.0D+0) GOTO 301
         ENDIF
         IF (E2 .LT. -30.0D+0) THEN
            GOTO 101
         ELSEIF (EIDN .GT. (E - 3.0D+0)) THEN
            GOTO 201
         ENDIF
      ELSEIF (EIUP .LT. 0.5D+0) THEN
         GOTO 101
      ENDIF
   301 CONTINUE
      IF ((E1 .LT. -30.0D+0) .AND. (E2 .LT. -30.0D+0)) THEN
         IF (EMAIN .EQ. 'STOP') THEN
            GOTO 900
         ELSE
            E = 1.0D+0
            write(*,*) 'l,ei,e1,e2', l, E, E1, E2
         ENDIF
      ELSEIF (E2 .LT. -30.0D+0) THEN
         IF (EMAIN .EQ. 'STOP') THEN
            GOTO 900
         ELSE
            E = MAX(E1,E)
            write(*,*) 'l,ei,e1,e2', l, E, E1, E2
         ENDIF
      ELSEIF (E1 .LT. -30.0D+0) THEN
         GOTO 900
      ELSE
         E = (E1+E2)*0.5D+0
         write(*,*) 'l,ei,e1,e2', l, E, E1, E2
      ENDIF
      RETURN
  900 continue
      write(*,'("no energy limits found for L=",i2)') l
      write(*,'("E-bottom ",g18.10,3x,"E-top ",g18.10)') e1,e2

case default
  write(*,*)
  write(*,'("Error(findband): No such method for search of linearization energies: ",a)') &
    trim(findlinentype)
  write(*,*)
  stop
end select

end subroutine
!EOC
