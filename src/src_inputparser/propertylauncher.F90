! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine propertylauncher
      Use modinput
      Use inputdom
      Use modmain, Only: task
      Use modmpi, Only: rank
      Implicit None
      integer :: l, a, b, c
      
      call delete_warnings

! properties which depend on the ground state only

      If (associated(input%properties%bandstructure)) Then
         call rereadinput
         ! tasks are: 20, 21
         task = 20
         Call bandstr
      End If
      If (associated(input%properties%fermisurfaceplot) .And. rank .Eq. 0) Then
         call rereadinput
         if (associated(input%properties%fermisurfaceplot%plot2d)) then
            task = 101
            Call fermisurf
         Else
            task = 100
            Call fermisurf
         End If
      End If
      If (associated(input%properties%dos) .And. rank .Eq. 0) Then
         call rereadinput
         task = 10
         Call dos
      End If
      If (associated(input%properties%wfplot) .And. rank .Eq. 0) Then
#define NOTSTM .false.
         Call wfplot (NOTSTM)
      End If
      If (associated(input%properties%stm) .And. rank .Eq. 0) Then
#define STM .true.
         Call stm
      End If
      If (associated(input%properties%LSJ) .And. rank .Eq. 0) Then
         call rereadinput
! read in input again to reset the magnetic moments for proper symmetry after
! a possible run of the groundstate
         task=15
         if (associated(input%properties%LSJ%kstlist)) task=16
         Call writelsj
      End If
      If (associated(input%properties%masstensor) .And. rank .Eq. 0) Then
         call rereadinput
         task=25
         Call effmass
      End If

      If (associated(input%properties%chargedensityplot) .And. rank .Eq. 0) Then
         Call rhoplot
      End If
      If (associated(input%properties%exccplot) .And. rank .Eq. 0) Then
         Call potplot
      End If
      If (associated(input%properties%elfplot) .And. rank .Eq. 0) Then
         Call elfplot
      End If
      If (associated(input%properties%xcmvecfield) .And. rank .Eq. 0) Then
         If (associated(input%properties%xcmvecfield%plot2d) .And. rank .Eq. 0) Then
            task = 82
            Call vecplot
         End If
         If (associated(input%properties%xcmvecfield%plot3d) .And. rank .Eq. 0) Then
            task = 83
            Call vecplot
         End If
      End If
      If (associated(input%properties%mvecfield)) Then
         If (associated(input%properties%mvecfield%plot2d) .And. rank .Eq. 0) Then
            task = 72
            Call vecplot
         End If
         If (associated(input%properties%mvecfield%plot3d) .And. rank .Eq. 0) Then
            task = 73
            Call vecplot
         End If
      End If
      If (associated(input%properties%electricfield)) Then
         If (associated(input%properties%electricfield%plot2d) .And. rank .Eq. 0) Then
            task = 142
            Call vecplot
         End If
         If (associated(input%properties%electricfield%plot3d) .And. rank .Eq. 0) Then
            task = 143
            Call vecplot
         End If
      End If
      If (associated(input%properties%gradmvecfield) .And. rank .Eq. 0) Then
         Call dbxcplot
      End If
      If (associated(input%properties%EFG) .And. rank .Eq. 0) Then
         Call writeefg
      End If
      If (associated(input%properties%mossbauer) .And. rank .Eq. 0) Then
         Call mossbauer
      End If
      If (associated(input%properties%expiqr) .And. rank .Eq. 0) Then
         call rereadinput
         Call writeexpiqr
      End If
      If (associated(input%properties%elnes) .And. rank .Eq. 0) Then
         Call elnes
      End If

! calculate and print the momentum matrix elements
      If (associated(input%properties%momentummatrix)) Then
         call rereadinput
         Call writepmat
      End If

! IP-RPA dielectric tensor      
      If (associated(input%properties%dielmat)) Then
         call rereadinput
         call dielmat
      End If    

! MOKE effect      
      If (associated(input%properties%moke)) Then
         call rereadinput
         call moke
      End If

! Nonlinear optics: Second Harmonic Generation 
      If (associated(input%properties%shg)) Then
         call rereadinput
         do l = 1, size(input%properties%shg%chicomp,2)
           a = input%properties%shg%chicomp(1,l)
           b = input%properties%shg%chicomp(2,l)
           c = input%properties%shg%chicomp(3,l)
           call shg(a,b,c)
         end do
      End If

      If (associated(input%properties%eliashberg) .And. rank .Eq. 0) Then
         ! this task depends on the results triggered by
         ! "input%properties%phonon"
         ! set the default values if dos element not present
         if (.not.associated(input%properties%dos)) &
           input%properties%dos => getstructdos (emptynode)
         ! electron-phonon coupling
         task=240
         call epcouple
         ! phonon linewidths
         task=245
         call phlwidth
         ! Eliashberg function
         task=250
         Call alpha2f
      End If

      ! Raman scattering      
      ! the subroutine raman triggers a phonon calculation, if requested there, and
      ! requires the input of element xs
      If (associated(input%properties%raman)) Then
         call raman
      End If

End Subroutine propertylauncher
