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
      Use mod_hartreefock, only: exnl
      Implicit None
      ! arguments
      Integer, Intent (In) :: fnum
      CHARACTER(LEN=6)     :: input_groundstate_verbosity

      input_groundstate_verbosity = 'high'
      If ( input_groundstate_verbosity .eq. 'low'    .or. &
         & input_groundstate_verbosity .eq. 'normal' .or. &
         & input_groundstate_verbosity .eq. 'high') Then
         write (fnum,*)
         write (fnum, '("% START ENERGY BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
         write (fnum,*)
         write (fnum, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
         write (fnum, '("  Fermi energy", T45, ": ", F22.12)') efermi
         write (fnum, '("  Total energy", T45, ": ", F22.12)') engytot
         write (fnum, '("--------------------------------------------------------------------")')
      End If
      If ( input_groundstate_verbosity .eq. 'normal' .or. &
         & input_groundstate_verbosity .eq. 'high') Then
         write (fnum,*)
         write (fnum, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
         write (fnum, '("Contributions to the total energy")')
         write (fnum, '("--------------------------------------------------------------------")')
         write (fnum, '("  Kinetic energy"                        , T45, ": ", F22.12)') engykn
         write (fnum, '("  Coulomb energy"                        , T45, ": ", F22.12)') engycl
         Write (fnum, '("  Exchange energy"                       , T45, ": ", F22.12)') engyx
         if (associated(input%groundstate%Hybrid)) then
           if (input%groundstate%Hybrid%exchangetypenumber == 1) Then         
             write (fnum, '("  Hartree-Fock energy"                       , T45, ": ", F22.12)') exnl
           end if 
         end if
         Write (fnum, '("  Correlation energy"                    , T45, ": ", F22.12)') engyc
         If (input%groundstate%chgexs .Ne. 0) Then
            write (fnum, '("  Correction to the background charge", T45, ": ", F22.12)') engycbc
         End If
         If(engylu .Ne. 0) then
            write (fnum, '("  LDA+U correction"                   , T45, ": ", F22.12)') engylu
         End if
         If (associated(input%groundstate%spin)) Then
            write (fnum, '("(External B-field energy excluded from the total energy)")')
         End If
         write (fnum, '("--------------------------------------------------------------------")')
      End If
      If (input_groundstate_verbosity .eq. 'high') Then
         write (fnum,*)
         write (fnum, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
         write (fnum, '("The kinetic energy is the difference of the following two terms")')
         write (fnum, '("--------------------------------------------------------------------")')
         write (fnum, '("  Sum of eigenvalues"        , T45, ": ", F22.12)') evalsum
         write (fnum, '("  Effective potential energy", T45, ": ", F22.12)') &
         & engyvcl + engyvxc + engybxc + engybext + engybmt
         write (fnum, '("--------------------------------------------------------------------")')

         write (fnum,*)
         write (fnum, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
         write (fnum, '("Contributions to the effective potential energy")')
         write (fnum, '("--------------------------------------------------------------------")')
         write (fnum, '("  Coulomb potential energy"                 , T45, ": ", F22.12)') engyvcl
         write (fnum, '("  xc potential energy"                      , T45, ": ", F22.12)') engyvxc
         If (associated(input%groundstate%spin)) Then
            Write (fnum, '("  xc effective B-field energy"           , T45, ": ", F22.12)') engybxc
            Write (fnum, '("  External B-field energy"               , T45, ": ", F22.12)') engybext
         End If
         If(engybmt .Ne. 0) then
            write (fnum, '("  Non-physical muffin-tin B-field energy", T45, ": ", F22.12)') engybmt
         End If
         write (fnum, '("--------------------------------------------------------------------")')

         write (fnum,*)
         write (fnum, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
         write (fnum, '("Contributions to the Coulomb energy")')
         write (fnum, '("--------------------------------------------------------------------")')
         write (fnum, '("  Hartree energy"         , T45, ": ", F22.12)') engyhar
         write (fnum, '("  Electron-nuclear energy", T45, ": ", F22.12)') engyen
         write (fnum, '("  Nuclear-nuclear energy" , T45, ": ", F22.12)') engynn
         write (fnum, '("--------------------------------------------------------------------")')

         write (fnum,*)
         write (fnum, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
         write (fnum, '("Further energies")')
         write (fnum, '("--------------------------------------------------------------------")')
         Write (fnum, '("  Madelung energy"             , T45, ": ", F22.12)') engymad
         Write (fnum, '("  Core-electron kinetic energy", T45, ": ", F22.12)') engykncr
         write (fnum, '("--------------------------------------------------------------------")')
      End If
      write (fnum,*)
      write (fnum, '("% END ENERGY BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
      !call flushifc(fnum)
      Return
End Subroutine
