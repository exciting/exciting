!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine writeengy (fnum)
      Use modmain
      Use modinput
      Use mod_hybrids, only: ihyb, exnl
      Implicit None
      ! arguments
      Integer, Intent (In) :: fnum
      character(len=77)    :: string

      write (fnum, '(" Total energy", T45, ": ", F18.8)') engytot
      write (fnum, '(" _______________________________________________________________")')
      write (fnum, '(" Fermi energy", T45, ": ", F18.8)') efermi
      If ( input%groundstate%outputlevelnumber > 0 ) Then
         write (fnum, '(" Kinetic energy"    , T45, ": ", F18.8)') engykn
         write (fnum, '(" Coulomb energy"    , T45, ": ", F18.8)') engycl
         if (associated(input%groundstate%Hybrid)) then
           if ((input%groundstate%Hybrid%exchangetypenumber==1).and.(ihyb>0)) then
             write (fnum, '(" Non-local exchange energy"   , T45, ": ", F18.8)') exnl
           end if
         end if
         Write (fnum, '(" Exchange energy"   , T45, ": ", F18.8)') engyx
         Write (fnum, '(" Correlation energy", T45, ": ", F18.8)') engyc
         If (input%groundstate%chgexs .Ne. 0) Then
            write (fnum, '(" Correction to the background charge", T45, ": ", F18.8)') engycbc
         End If
         If(engylu .Ne. 0) then
            write (fnum, '(" LDA+U correction"                   , T45, ": ", F18.8)') engylu
         End if
         If (associated(input%groundstate%spin)) Then
            write (fnum, '(" No external B-field energy in total energy")')
         End If
      End If
      If ( input%groundstate%outputlevelnumber > 1 ) Then
         write (fnum, '(" Sum of eigenvalues"        , T45, ": ", F18.8)') evalsum
         write (fnum, '(" Effective potential energy", T45, ": ", F18.8)') &
         & engyvcl + engyvxc + engybxc + engybext + engybmt
         write (fnum, '(" Coulomb potential energy"      , T45, ": ", F18.8)') engyvcl
         write (fnum, '(" xc potential energy"           , T45, ": ", F18.8)') engyvxc
         If (associated(input%groundstate%spin)) Then
            Write (fnum, '(" xc effective B-field energy", T45, ": ", F18.8)') engybxc
            Write (fnum, '(" External B-field energy"    , T45, ": ", F18.8)') engybext
         End If
         If(engybmt .Ne. 0) then
            write (fnum, '(" Non-physical muffin-tin B-field energy", T45, ": ", F18.8)') engybmt
         End If
         write (fnum, '(" Hartree energy"         , T45, ": ", F18.8)') engyhar
         write (fnum, '(" Electron-nuclear energy", T45, ": ", F18.8)') engyen
         write (fnum, '(" Nuclear-nuclear energy" , T45, ": ", F18.8)') engynn
         Write (fnum, '(" Madelung energy"             , T45, ": ", F18.8)') engymad
         Write (fnum, '(" Core-electron kinetic energy", T45, ": ", F18.8)') engykncr
         if (associated(input%groundstate%dfthalf)) then
           Write (fnum, '(" DFT-1/2 contribution to total energy", T45, ": ", F18.8)') engyhalf
         endif
      End If
      If ( tlast .And. input%groundstate%vdWcorrection .Ne. "none" ) Then
         If ( input%groundstate%vdWcorrection .Eq. "DFTD2" ) Then
            Write (fnum, '(" DFT-D2 dispersion correction", T45, ": ", F18.8)') e_disp
         Else If ( input%groundstate%vdWcorrection .Eq. "TSvdW" ) Then
            Write (fnum, '(" TS-vdW dispersion correction", T45, ": ", F18.8)') e_disp
         End If
      End If

      call flushifc(fnum)
      Return
End Subroutine
