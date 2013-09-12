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
      Implicit None
      ! arguments
      Integer, Intent (In) :: fnum
      character(len=77)    :: string

      write (fnum, '(" Total energy", T45, ": ", F22.12)') engytot
      write (fnum, '(" ___________________________________________________________________")')
      write (fnum, '(" Fermi energy", T45, ": ", F22.12)') efermi

      If ( input%groundstate%outputlevelnumber > 0 ) Then
         write (fnum, '(" Kinetic energy"    , T45, ": ", F22.12)') engykn
         write (fnum, '(" Coulomb energy"    , T45, ": ", F22.12)') engycl
         Write (fnum, '(" Exchange energy"   , T45, ": ", F22.12)') engyx
         Write (fnum, '(" Correlation energy", T45, ": ", F22.12)') engyc
         If (input%groundstate%chgexs .Ne. 0) Then
            write (fnum, '(" Correction to the background charge", T45, ": ", F22.12)') engycbc
         End If
         If(engylu .Ne. 0) then
            write (fnum, '(" LDA+U correction"                   , T45, ": ", F22.12)') engylu
         End if
         If (associated(input%groundstate%spin)) Then
            write (fnum, '(" (External B-field energy excluded from the total energy)")')
         End If
      End If
      If ( input%groundstate%outputlevelnumber > 1 ) Then
         write (fnum, '(" Sum of eigenvalues"        , T45, ": ", F22.12)') evalsum
         write (fnum, '(" Effective potential energy", T45, ": ", F22.12)') &
         & engyvcl + engyvxc + engybxc + engybext + engybmt
         write (fnum, '(" Coulomb potential energy"      , T45, ": ", F22.12)') engyvcl
         write (fnum, '(" xc potential energy"           , T45, ": ", F22.12)') engyvxc
         If (associated(input%groundstate%spin)) Then
            Write (fnum, '(" xc effective B-field energy", T45, ": ", F22.12)') engybxc
            Write (fnum, '(" External B-field energy"    , T45, ": ", F22.12)') engybext
         End If
         If(engybmt .Ne. 0) then
            write (fnum, '(" Non-physical muffin-tin B-field energy", T45, ": ", F22.12)') engybmt
         End If
         write (fnum, '(" Hartree energy"         , T45, ": ", F22.12)') engyhar
         write (fnum, '(" Electron-nuclear energy", T45, ": ", F22.12)') engyen
         write (fnum, '(" Nuclear-nuclear energy" , T45, ": ", F22.12)') engynn
         Write (fnum, '(" Madelung energy"             , T45, ": ", F22.12)') engymad
         Write (fnum, '(" Core-electron kinetic energy", T45, ": ", F22.12)') engykncr
      End If
      call flushifc(fnum)
      Return
End Subroutine