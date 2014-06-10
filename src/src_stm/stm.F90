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
  Integer :: ik, ist
  Real (8) :: x, y, t1, bias
  type(plotlabels),pointer ::labels
  ! allocatable arrays
  Complex (8), Allocatable :: evecfv (:, :)
  Complex (8), Allocatable :: evecsv (:, :)
  ! external functions
  External :: occstm
  ! initialise universal variables
  Call init0
  Call init1
  Allocate (evecfv(nmatmax, nstfv))
  Allocate (evecsv(nstsv, nstsv))
  ! read the density and potentials from file
  Call readstate
  ! read Fermi energy from file
  Call readfermi
  ! find the new linearisation energies
  Call linengy
  ! generate the APW radial functions
  Call genapwfr
  ! generate the local-orbital radial functions
  Call genlofr
  ! set the occupancies

  bias = input%properties%stm%bias

  ! Setting occupations for STM imaging in the Tersoff-Hamann approx (PRB 31,805 (1985)).
  If ( input%properties%stm%stmtype.Eq.'differentialConductance' &
       & .And. input%properties%stm%stmmode.Eq. 'constantHeight') Then
     Write(*,*)
     Write (*, '("Info(stm):")')
     Write (*, '("Setting occupations for constant-height/differential conductance.")')
     Call genocc(efermi+bias,efermi+bias)
  Else If ( input%properties%stm%stmtype.Eq.'integratedLDOS' &
       & .And. input%properties%stm%stmmode.Eq. 'constantHeight') Then
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
     call warning('Error(stm): STM still not implemented for direct topographic plot.')
     stop 'Error(stm): STM still not implemented for direct topographic plot.'
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

  ! write the wavefunction modulus squared plot to file

  ! plotdef%parallelogram%grid(1)
  ! plotdef%parallelogram%grid(2)
  ! plotdef%parallelogram%pointarray(1)%point%coord
  ! plotdef%parallelogram%pointarray(2)%point%coord
  ! plotdef%parallelogram%origin%coord

  labels=>create_plotlablels("2D STM image","STM2d",2)
  call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,3,"STM","","graceunit")
  Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
       & rhomt, rhoir, input%properties%stm%plot2d)
  call destroy_plotlablels(labels)
  If ( input%properties%stm%stmmode.Eq.'topographic') Then
     Write (*,*)
     Write (*, '("Info(stm):")')
     Write (*, '("STM still not implemented for direct topographic plot.")')
     Write (*, '("For topographic plot generation consider to make a series &
          constant-height calculations at different heights and postprocess the &
          output to find the iso-surface. ")')
  Else
     Write (*,*)
     Write (*, '("Info(stm):")')
     Write (*, '(" 2D STM image written to STM2d.xml")')
  End If

  Write (*,*)
  Deallocate (evecfv, evecsv)
  Return
End Subroutine stm
