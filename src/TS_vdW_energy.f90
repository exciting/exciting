Subroutine TS_vdW_energy
  !Tkatchenko-Scheffler van der Waals correction 
  Use mod_energy, Only: e_disp
  Use TS_vdW_module, Only: get_TS_parameters, s6, rs6, damping_const, cutoff, C6ab, R0_eff_ab
  Use vdw_general_routines, Only: vdw_energy_pairwiseC6
  Use mod_atoms, Only: natmtot
  Use modinput
  Implicit None
!  Real(8), Allocatable :: C6ab(:,:), R0_eff_ab(:,:)
  If ( input%groundstate%do .Eq. "skip" ) Then
     Call init0
     ! read density from file
     Call readstate
  End If

  If (.Not. Allocated(C6ab)) Then
     Allocate(C6ab(natmtot,natmtot), R0_eff_ab(natmtot,natmtot))
     Call get_TS_parameters()
  End If

  e_disp = vdw_energy_pairwiseC6(s6, rs6, damping_const, cutoff, C6ab, R0_eff_ab)
  If (associated(input%properties)) Then
     If (associated(input%properties%TS_vdW)) Then
        Open(246, File="E_TS_vdW.OUT")
        Write(246,'(F18.8)')e_disp
        Close(246)
     End If
  End If
End Subroutine TS_vdW_energy
