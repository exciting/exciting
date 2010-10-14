!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phdos
      Use modmain
      Use modinput
      Use FoX_wxml
      Implicit None
! local variables
! number of temperature values
      Integer, Parameter :: ntemp = 10
      Integer :: n, iq, i, iw
      Integer :: i1, i2, i3
      Real (8) :: wmin, wmax, wd, dw
      Real (8) :: tmax, temp (ntemp), s (ntemp)
      Real (8) :: v (3), t1, t2
      Type (xmlf_t), Save :: xf
      Character (256) :: buffer
! allocatable arrays
      Real (8), Allocatable :: wp (:)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: gw (:)
      Real (8), Allocatable :: f (:), g (:), cf (:, :)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: dynr (:, :, :)
      Complex (8), Allocatable :: dynp (:, :)
      Complex (8), Allocatable :: ev (:, :)
! initialise universal variables
      Call init0
      Call init2
      n = 3 * natmtot
      Allocate (wp(n))
      Allocate (w(input%phonons%phonondos%nwdos))
      Allocate (gw(input%phonons%phonondos%nwdos))
      Allocate (f(input%phonons%phonondos%nwdos), &
     & g(input%phonons%phonondos%nwdos), cf(3, &
     & input%phonons%phonondos%nwdos))
      Allocate (dynq(n, n, nqpt))
      Allocate (dynr(n, n, ngridq(1)*ngridq(2)*ngridq(3)))
      Allocate (dynp(n, n))
      Allocate (ev(n, n))
! read in the dynamical matrices
      Call readdyn (.true.,dynq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! Fourier transform the dynamical matrices to real-space
      Call dynqtor (dynq, dynr)
! find the minimum and maximum frequencies
      wmin = 0.d0
      wmax = 0.d0
      Do iq = 1, nqpt
         Call dyndiag (dynq(:, :, iq), wp, ev)
         wmin = Min (wmin, wp(1))
         wmax = Max (wmax, wp(n))
      End Do
      wmax = wmax + (wmax-wmin) * 0.1d0
      wmin = wmin - (wmax-wmin) * 0.1d0
      wd = wmax - wmin
      If (wd .Lt. 1.d-8) wd = 1.d0
      dw = wd / dble (input%phonons%phonondos%nwdos)
! generate energy grid
      Do iw = 1, input%phonons%phonondos%nwdos
         w (iw) = dw * dble (iw-1) + wmin
      End Do
      gw (:) = 0.d0
      Do i1 = 0, input%phonons%phonondos%ngrdos - 1
         v (1) = dble (i1) / dble (input%phonons%phonondos%ngrdos)
         Do i2 = 0, input%phonons%phonondos%ngrdos - 1
            v (2) = dble (i2) / dble (input%phonons%phonondos%ngrdos)
            Do i3 = 0, input%phonons%phonondos%ngrdos - 1
               v (3) = dble (i3) / dble (input%phonons%phonondos%ngrdos)
! compute the dynamical matrix at this particular q-point
               Call dynrtoq (v, dynr, dynp)
! find the phonon frequencies
               Call dyndiag (dynp, wp, ev)
               Do i = 1, n
                  t1 = (wp(i)-wmin) / dw + 1.d0
                  iw = Nint (t1)
                  If ((iw .Ge. 1) .And. (iw .Le. &
                 & input%phonons%phonondos%nwdos)) Then
                     gw (iw) = gw (iw) + 1.d0
                  End If
               End Do
            End Do
         End Do
      End Do
      t1 = 1.d0 / (dw*dble(input%phonons%phonondos%ngrdos)**3)
      gw (:) = t1 * gw (:)
! smooth phonon DOS if required
      If (input%phonons%phonondos%nsmdos .Gt. 0) Call fsmooth &
     & (input%phonons%phonondos%nsmdos, input%phonons%phonondos%nwdos, 1, gw)
! write phonon DOS to file
      Open (50, File='PHDOS.OUT', Action='WRITE', Form='FORMATTED')
      Do iw = 1, input%phonons%phonondos%nwdos
         Write (50, '(2G18.10)') w (iw), gw (iw)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(phdos):")')
      Write (*, '(" phonon density of states written to PHDOS.OUT")')
!-------------------------------------------!
!     thermodynamic properties from DOS     !
!-------------------------------------------!
! maximum temperature
      tmax = wmax / kboltz
! temperature grid
      Do i = 1, ntemp
         temp (i) = tmax * dble (i) / dble (ntemp)
      End Do
      Open (50, File='THERMO.OUT', Action='WRITE', Form='FORMATTED')
      Call xml_OpenFile ("thermo.xml", xf, replace=.True., pretty_print=.True.)
      Call xml_NewElement (xf, "thermodynamicproperties")
! zero point energy
      Do iw = 1, input%phonons%phonondos%nwdos
         f (iw) = gw (iw) * w (iw)
      End Do
      Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
      t1 = 0.5d0 * dble (natmtot) * g (input%phonons%phonondos%nwdos)
      Write (50,*)
      Write (50, '("Zero-point energy : ", G18.10)') t1
      Call xml_NewElement (xf, "zeropointenergy")
      Write (buffer, '(g18.10)') t1
      Call xml_AddAttribute (xf, "name", "zero-point energy")
      Call xml_AddAttribute (xf, "value", trim(adjustl(buffer)))
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "zeropointenergy")
! vibrational energy
      Write (50,*)
      Write (50, '("Vibrational energy vs. temperature :")')
      Call xml_NewElement (xf, "vibrationalenergy")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "vibrational energy")
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, ntemp
         Do iw = 1, input%phonons%phonondos%nwdos
            t1 = w (iw) / (2.d0*kboltz*temp(i))
            If (t1 .Gt. 0.d0) Then
               f (iw) = gw (iw) * w (iw) * Cosh (t1) / Sinh (t1)
            Else
               f (iw) = 0.d0
            End If
         End Do
         Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
         t1 = 0.5d0 * dble (natmtot) * g (input%phonons%phonondos%nwdos)
         Write (50, '(2G18.10)') temp (i), t1
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') t1
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
         s (i) = t1
      End Do
      Call xml_endElement (xf, "vibrationalenergy")
! free energy
      Write (50,*)
      Write (50, '("Free energy vs. temperature :")')
      Call xml_NewElement (xf, "vibrationalfreeenergy")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "vibrational free energy")
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, ntemp
         Do iw = 1, input%phonons%phonondos%nwdos
            t1 = 2.d0 * Sinh (w(iw)/(2.d0*kboltz*temp(i)))
            If (t1 .Gt. 0.d0) Then
               f (iw) = gw (iw) * Log (t1)
            Else
               f (iw) = 0.d0
            End If
         End Do
         Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
         t1 = dble (natmtot) * kboltz * temp (i) * g &
        & (input%phonons%phonondos%nwdos)
         Write (50, '(2G18.10)') temp (i), t1
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') t1
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
! compute entropy from S = (F-E)/T
         s (i) = (s(i)-t1) / temp (i)
      End Do
      Call xml_endElement (xf, "vibrationalfreeenergy")
! entropy
      Write (50,*)
      Write (50, '("Entropy vs. temperature :")')
      Call xml_NewElement (xf, "vibrationalentropy")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "vibrational entropy")
      Call xml_AddAttribute (xf, "unit", "Hartree/Kelvin")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, ntemp
         Write (50, '(2G18.10)') temp (i), s (i)
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') s(i)
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
      End Do
      Call xml_endElement (xf, "vibrationalentropy")
! heat capacity
      Write (50,*)
      Write (50, '("Heat capacity vs. temperature :")')
      Call xml_NewElement (xf, "heatcapacity")
      Call xml_NewElement (xf, "mapdef")
      Call xml_NewElement (xf, "variable1")
      Call xml_AddAttribute (xf, "name", "temperature")
      Call xml_AddAttribute (xf, "unit", "Kelvin")
      Call xml_endElement (xf, "variable1")
      Call xml_NewElement (xf, "function1")
      Call xml_AddAttribute (xf, "name", "heat capacity")
      Call xml_AddAttribute (xf, "unit", "Hartree")
      Call xml_endElement (xf, "function1")
      Call xml_endElement (xf, "mapdef")
      Do i = 1, ntemp
         Do iw = 1, input%phonons%phonondos%nwdos
            t1 = w (iw) / (kboltz*temp(i))
            t2 = Exp (t1) - 1.d0
            If (t2 .Ne. 0.d0) Then
               f (iw) = gw (iw) * (t1**2) * (t2+1.d0) / t2 ** 2
            Else
               f (iw) = 0.d0
            End If
         End Do
         Call fderiv (-1, input%phonons%phonondos%nwdos, w, f, g, cf)
         t1 = dble (natmtot) * kboltz * g (input%phonons%phonondos%nwdos)
         Write (50, '(2G18.10)') temp (i), t1
         Call xml_NewElement (xf, "map")
         Write (buffer, '(4g18.10)') temp(i)
         Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
         Write (buffer, '(4g18.10)') t1
         Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
         Call xml_endElement (xf, "map")
      End Do
      Call xml_endElement (xf, "heatcapacity")
      Close (50)
      Call xml_endElement (xf, "thermodynamicproperties")
      Call xml_Close (xf)
      Write (*, '(" thermodynamic properties written to THERMO.OUT")')
      Write (*,*)
      Deallocate (wp, w, gw, f, g, cf, dynq, dynr, dynp, ev)
      Return
End Subroutine
