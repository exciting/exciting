!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: findband
! !INTERFACE:
!
Subroutine findband (findlinentype, l, k, np, nr, r, vr, de0, etol, e, tfnd)
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
      Implicit None
  ! arguments
      Character (*), Intent (In) :: findlinentype
      Integer, Intent (In) :: l
      Integer, Intent (In) :: k
      Integer, Intent (In) :: np
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Real (8), Intent (In) :: de0
      real(8), intent(in) :: etol
      Real (8), Intent (Inout) :: e
      logical, intent(out) :: tfnd
  ! local variables
  ! maximum number of steps
      Integer, Parameter :: maxstp = 1000
      Character * 4, Parameter :: emain = 'continue'
      Real (8), Parameter :: etoolow = - 1000.d0
      Real (8), Parameter :: ecutlow = - 30.d0
      Real (8), Parameter :: efermibands = 0.5d0
      Real (8), Parameter :: erangebands = 2.d0
      Real (8), Parameter :: ediffusebands = efermibands + erangebands
      Real (8), Parameter :: erange1 = 3.d0
      Real (8), Parameter :: edefault1 = 1.d0
      Integer :: ie, nn
      Real (8) :: de, et, eb, t, tp
      Real (8) :: u, uup, udn, upold, ddnold, dnold, du, dudn, dupold, &
     & duup, utst, dutst
      Real (8) :: e1, e2, eidn, eiup
  ! automatic arrays
      Real (8) :: p0 (nr), p1 (nr), q0 (nr), q1 (nr)
      tfnd=.false.
      Select Case (trim(findlinentype))
      Case ('simple')
         tp = 0.d0
     ! find the top of the band
         de = Abs (de0)
         et = e
         Do ie = 1, maxstp
            et = et + de
            Call rschroddme (0, l, k, et, np, nr, r, vr, nn, p0, p1, &
           & q0, q1)
            t = p0 (nr)
            If (ie .Gt. 1) Then
               If (t*tp .Le. 0.d0) Then
                  If (Abs(de) .Lt. etol) Go To 10
                  de = - 0.5d0 * de
               End If
            End If
            tp = t
         End Do
         Go To 30
10       Continue
     ! find the bottom of the band
         de = - Abs (de0)
         eb = et + 5.d0 * Abs (de0)
         Do ie = 1, maxstp
            eb = eb + de
            Call rschroddme (0, l, k, eb, np, nr, r, vr, nn, p0, p1, &
           & q0, q1)
            t = p1 (nr)
            If (ie .Gt. 1) Then
               If (t*tp .Le. 0.d0) Then
                  If (Abs(de) .Lt. etol) Go To 20
                  de = - 0.5d0 * de
               End If
            End If
            tp = t
         End Do
         Go To 30
20       Continue
     ! set the band energy to the mid-point
         e = (et+eb) / 2.d0
         tfnd=.true.
30       Continue
      Case ('advanced')
         de = de0
         Call rschroddme (0, l, k, e, np, nr, r, vr, nn, p0, p1, q0, &
        & q1)
         u = p0 (nr)
         du = p1 (nr)
         eiup = e
         eidn = e
         e1 = etoolow
         e2 = etoolow
         upold = u
         dupold = du
         dnold = u
         ddnold = du
11       Continue
         eiup = eiup + de
         Call rschroddme (0, l, k, eiup, np, nr, r, vr, nn, p0, p1, q0, &
        & q1)
         uup = p0 (nr)
         duup = p1 (nr)
         utst = upold * uup
         dutst = dupold * duup
         upold = uup
         dupold = duup
         If (utst .Lt. 0.d0) Then
            e2 = eiup
            If (e1 .Gt. ecutlow) Go To 31
         End If
         If (dutst .Lt. 0.d0) Then
            e1 = eiup
            If (e2 .Gt. ecutlow) Go To 31
         End If
         If ((e1 .Lt. ecutlow) .And. (eiup .Lt. ediffusebands)) Then
21          Continue
            eidn = eidn - de
            Call rschroddme (0, l, k, eidn, np, nr, r, vr, nn, p0, p1, &
           & q0, q1)
            udn = p0 (nr)
            dudn = p1 (nr)
            utst = dnold * udn
            dutst = ddnold * dudn
            dnold = udn
            ddnold = dudn
            If (utst .Lt. 0.d0) Then
               e2 = eidn
               If (e1 .Gt. ecutlow) Go To 31
            End If
            If (dutst .Lt. 0.d0) Then
               e1 = eidn
               If (e2 .Gt. ecutlow) Go To 31
            End If
            If (e2 .Lt. ecutlow) Then
               Go To 11
            Else If (eidn .Gt. (e-erange1)) Then
               Go To 21
            End If
         Else If (eiup .Lt. efermibands) Then
            Go To 11
         End If
31       Continue
         If ((e1 .Lt. ecutlow) .And. (e2 .Lt. ecutlow)) Then
            If (emain .Eq. 'stop') Then
               Go To 40
            Else
               e = edefault1
            End If
         Else If (e2 .Lt. ecutlow) Then
            If (emain .Eq. 'stop') Then
               Go To 40
            Else
               e = Max (e1, e)
            End If
         Else If (e1 .Lt. ecutlow) Then
            Go To 40
         Else
            e = (e1+e2) * 0.5d0
         End If
         tfnd=.true.
         Return
40       Continue
         Write (*,*)
         Write (*, '("Warning(findband): no energy limits found for l=",i&
        &2)') l
         Write (*, '("E-bottom ",g18.10,4x,"E-top ",g18.10)') e1, e2
         tfnd=.false.
         return
      Case Default
         Write (*,*)
         Write (*, '("Error(findband): No such method for search of lin&
        &earization energies: ",a)') trim (findlinentype)
         Write (*,*)
         Stop
      End Select
End Subroutine findband
!EOC
