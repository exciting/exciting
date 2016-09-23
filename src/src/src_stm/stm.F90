!
!
!
!
! Copyright (C) 2014 S. Rigamonti and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine stm
  Use modmain
  Use modinput
  use modplotlabels
  Implicit None
  ! local variables
  Integer :: ik
  Real (8) ::  bias
  type(plotlabels),pointer ::labels
  ! allocatable arrays
  Complex (8), Allocatable :: evecfv (:, :)
  Complex (8), Allocatable :: evecsv (:, :)
  Character(256) :: string
  ! initialise universal variables
  Call init0
  Call init1
  Allocate (evecfv(nmatmax, nstfv))
  Allocate (evecsv(nstsv, nstsv))
  ! read the density and potentials from file
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
  ! in case of HF hybrids use PBE potential
            string=filext
            filext='_PBE.OUT'
            Call readstate
            filext=string
           Else
               Call readstate
           End If
        Else         
           Call readstate
        End If 
  ! read Fermi energy from file
  Call readfermi
  ! find the new linearisation energies
  Call linengy
  ! generate the APW radial functions
  Call genapwfr
  ! generate the local-orbital radial functions
  Call genlofr
  ! update potential in case if HF Hybrids
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
               Call readstate
           End If
        End If 
  ! verify consistency of input for stm
  Call checkinput
  ! set the occupancies

  bias = input%properties%stm%bias

  ! Setting occupations for STM imaging in the Tersoff-Hamann approx (PRB 31,805 (1985)).
  If ( input%properties%stm%stmtype.Eq.'differentialConductance') Then
     Write(*,*)
     Write (*, '("Info(stm):")')
     Write (*, '("Setting occupations for constant-height/differential conductance.")')
     Call genocc(efermi+bias,efermi+bias)
  Else If ( input%properties%stm%stmtype.Eq.'integratedLDOS') Then
     Write(*,*)
     Write (*, '("Info(stm):")')
     Write (*, '("Setting occupations for constant-height/Integrated LDOS.")')
     If (bias.lt.0.0d0) Then
        Write(*,*) 'negative bias'
        Call genocc(efermi+bias,efermi)
     Else
        Call genocc(efermi,efermi+bias)
     End If
  Else
     call warning('Error(stm): unknown stmtype value.')
     stop 'Error(stm): unknown stmtype value.'
  End If
  ! set the charge density to zero
  rhomt (:, :, :) = 0.d0
  rhoir (:) = 0.d0
  ! compute the charge density with the new occupancies
  Do ik = 1, nkpt
     ! get the eigenvectors from file
     Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
     Call getevecsv (vkl(:, ik), evecsv)
     Call rhovalk (ik, evecfv, evecsv)
     Call genrhoir (ik, evecfv, evecsv)
  End Do
  ! symmetrise the density for the STM plot
  Call symrf (input%groundstate%lradstep, rhomt, rhoir)
  ! convert the density from a coarse to a fine radial mesh
  Call rfmtctof (rhomt)

  If( input%properties%stm%stmmode.Eq. 'constantHeight') Then
     Call genplot2d("STM2d")
  Else If ( input%properties%stm%stmmode.Eq. 'topographic') Then
     Call genplot3d("STM3d")
  End If

  Write (*,*)
  Write (*, '("Info(stm):")')
  Write (*, '(" STM image written ")')
  Write (*,*)

  Deallocate (evecfv, evecsv)
  Return
End Subroutine stm
