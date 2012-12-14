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
      Real (8), Parameter :: ecutlow = - 100.d0
      Real (8), Parameter :: efermibands = 0.5d0
      Real (8), Parameter :: erangebands = 2.d0
      Real (8), Parameter :: ediffusebands = efermibands + erangebands
      Real (8), Parameter :: erange1 = 3.d0
      Real (8), Parameter :: edefault1 = 1.d0
      Integer :: ie, nn
      Real (8) :: de, et, eb, t0, t1, pt, pb, dl, dl0
      Real (8) :: u, uup, udn, upold, ddnold, dnold, du, dudn, dupold, &
     & duup, utst, dutst
      Real (8) :: e1, e2, eidn, eiup
  ! automatic arrays
      Real (8) :: p0 (nr), p1 (nr), q0 (nr), q1 (nr)
      character(10) :: fname

      tfnd=.false.

      de = Abs(de0)

!     visualize the logarithmic derivative D_l
      write(fname,'("dl_l=",i1,".dat")') l
      write(*,*) fname
      open(777,file=fname,action='write')
      eb = -10.d0
      do while (eb.le.30.d0)
        call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
        t0 = p0(nr)
        t1 = p1(nr)
        if (dabs(t0)>1.0d-6) then
          ! shifted value of the logarithmic derivative
          dl = r(nr)*t1/t0+(l+1)/r(nr)
          write(777,*) eb, dl
        end if
        eb = eb + de
      end do ! eb
      close(777)


      Select Case (trim(findlinentype))

!-------------------------------------------
      Case ('simple')
!-------------------------------------------
         de = Abs(de0)

!        find the Linearization Energy from the equation D_{l}(E)=-(l+1)/R_{MT}
         eb = e
         
         call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
         t0 = p0(nr)
         t1 = p1(nr)
         if (dabs(t0)>1.0d-6) then
           dl0 = r(nr)*t1/t0+(l+1)/r(nr)
         else
           dl0 = 1000.d0*dsign(1.d0,t0)*dsign(1.d0,t1)
         end if
        
         do while (eb < 30.0d0)
           eb = eb + de
           call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
           t0 = p0(nr)
           t1 = p1(nr)
           if (dabs(t0)>1.0d-6) then
             ! shifted value of the logarithmic derivative
             dl = r(nr)*t1/t0+(l+1)/r(nr)
             if (dl*dl0.le.0.d0) then
               dl0 = dl
               if (dabs(dl).le.1.0d0) then
                 e = eb
                 exit
               end if
             end if
           end if
         end do

!        eb = -10.0d0
!        call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
!        pt = p0(nr)
!        pb = p1(nr)
!        
!        ! search up in the energy
!        do while (eb < 30.0d0)
!           
!           eb = eb + de
!           Call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
!          
!           t0 = p0(nr)
!           t1 = p1(nr)
!           if (t1*pb .Le. 0.d0) then
!             e1 = eb
!             write(*,*) 'Node DOWN found: Ebottom=', e1
!             pb = t1
!             pt = t0
!             ! search down in the energy
!             do while (eb < 30.0d0)
!              eb = eb + de
!              Call rschroddme(0, l, k, eb, np, nr, r, vr, nn, p0, p1, q0, q1)
!              t0 = p0(nr)
!              if (t0*pt .Le. 0.d0) then
!                e2 = eb
!                write(*,*) 'Node UP found: Etop=', e2
!                write(*,*) 'l=', l, '    LE=', 0.5d0*(e1+e2)
!                exit
!              end if
!             end do
!           end if
!            
!        end do ! eb
         
         tfnd=.true.
30       Continue

!-------------------------------------------
      Case ('advanced')
!-------------------------------------------

         de = de0
         Call rschroddme(0, l, k, e, np, nr, r, vr, nn, p0, p1, q0, q1)
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
         Call rschroddme(0, l, k, eiup, np, nr, r, vr, nn, p0, p1, q0, q1)
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
            Call rschroddme(0, l, k, eidn, np, nr, r, vr, nn, p0, p1, q0, q1)
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
         Write (100,*)
         Write (100, '("Warning(findband): no energy limits found for l=",i&
        &2)') l
         Write (100, '("E-bottom ",g18.10,4x,"E-top ",g18.10)') e1, e2
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
