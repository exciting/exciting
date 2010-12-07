
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine propertylauncher
      Use modinput
      Use inputdom
      Use modmain, Only: task
      Use modmpi, Only: rank
      Implicit None

! properties which depend on the ground state only

      If (associated(input%properties%bandstructure) .And. rank .Eq. 0) Then
         ! tasks are: 20, 21
         task = 20
         Call bandstr
      End If
      If (associated(input%properties%dos) .And. rank .Eq. 0) Then
         task = 10
         Call dos
      End If
      If (associated(input%properties%wfplot) .And. rank .Eq. 0) Then
#define NOTSTM .false.
         Call wfplot (NOTSTM)
      End If
      If (associated(input%properties%STM) .And. rank .Eq. 0) Then
#define STM .true.
         Call wfplot (STM)
      End If
      If (associated(input%properties%LSJ) .And. rank .Eq. 0) Then
         task=15
         if (associated(input%properties%LSJ%kstlist)) task=16
         Call writelsj
      End If
      If (associated(input%properties%masstensor) .And. rank .Eq. 0) Then
         task=25
         Call effmass
      End If
      If (associated(input%properties%fermisurfaceplot) .And. rank .Eq. 0) Then
         If (input%properties%fermisurfaceplot%separate) Then
            task = 101
            Call fermisurf
         Else
            task = 100
            Call fermisurf
         End If
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
         Call writeexpiqr
      End If
      If (associated(input%properties%elnes) .And. rank .Eq. 0) Then
         Call elnes
      End If
      If (associated(input%properties%momentummatrix) .And. rank .Eq. 0) Then
         task = 120
         Call writepmatxs
      End If

! properties which depend on the ground state and/or on the outputs of other properties

      If (associated(input%properties%dielectric) .And. rank .Eq. 0) Then
         ! set the default values if dos element not present
         if (.not.associated(input%properties%dos)) &
           input%properties%dos => getstructdos (emptynode)
         ! this task depends on the results triggered by
         ! "input%properties%momentummatrix"
         task = 121
         Call dielectric
      End If
      If (associated(input%properties%moke) .And. rank .Eq. 0) Then
         ! set the default values if dos element not present
         if (.not.associated(input%properties%dos)) &
           input%properties%dos => getstructdos (emptynode)
         ! set the default values if dielectric element not present
         if (.not.associated(input%properties%dielectric)) &
           input%properties%dielectric => getstructdielectric (emptynode)
         ! this task depends on the results triggered by
         ! "input%properties%momentummatrix"
         Call moke
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

End Subroutine propertylauncher
