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
      Write (fnum,*)
      Write (fnum, '("Energies :")')
      Write (fnum, '(" Fermi", T30, ": ", G22.12)') efermi
      Write (fnum, '(" sum of eigenvalues", T30, ": ", G22.12)') &
     & evalsum
      Write (fnum, '(" electronic kinetic", T30, ": ", G22.12)') engykn
      Write (fnum, '(" core electron kinetic", T30, ": ", G22.12)') &
     & engykncr
      Write (fnum, '(" Coulomb", T30, ": ", G22.12)') engycl
      Write (fnum, '(" Coulomb potential", T30, ": ", G22.12)') engyvcl
      Write (fnum, '(" nuclear-nuclear", T30, ": ", G22.12)') engynn
      Write (fnum, '(" electron-nuclear", T30, ": ", G22.12)') engyen
      Write (fnum, '(" Hartree", T30, ": ", G22.12)') engyhar
      Write (fnum, '(" Madelung", T30, ": ", G22.12)') engymad
      If (input%groundstate%chgexs .Ne. 0.d0) Then
         Write (fnum, '(" comp. background charge", T30, ": ", G22.12)') engycbc
      End If
      Write (fnum, '(" xc potential", T30, ": ", G22.12)') engyvxc
      If (associated(input%groundstate%spin)) Then
         Write (fnum, '(" xc effective B-field", T30, ": ", G22.12)') &
        & engybxc
         Write (fnum, '(" external B-field", T30, ": ", G22.12)') &
        & engybext
      End If
      Write (fnum, '(" exchange", T30, ": ", G22.12)') engyx
      Write (fnum, '(" correlation", T30, ": ", G22.12)') engyc
      If (ldapu .Ne. 0) Then
         Write (fnum, '(" LDA+U", T30, ": ", G22.12)') engylu
      End If
      Write (fnum, '(" total energy", T30, ": ", G22.12)') engytot
      If (associated(input%groundstate%spin)) Then
         Write (fnum, '(" (external B-field energy excluded from total)&
        &")')
      End If
!call flushifc(fnum)
      Return
End Subroutine
