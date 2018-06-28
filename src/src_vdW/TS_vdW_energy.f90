Subroutine TS_vdW_energy
  !Tkatchenko-Scheffler van der Waals correction 
  Use mod_energy, Only: e_disp
  Use TS_vdW_module, Only: get_TS_parameters, C6ab, R0_eff_ab!, s6, sr6, damping_const, cutoff
  Use vdw_general_routines, Only: vdw_energy_pairwiseC6, set_default_vdW_parameters
  Use mod_atoms, Only: natmtot
  Use modinput
  Implicit None
  !  Real(8), Allocatable :: C6ab(:,:), R0_eff_ab(:,:)

  Call set_default_vdW_parameters
  
  If ( input%groundstate%do .Eq. "skip" ) Then
     Call init0
     ! read density from file
     Call readstate
  End If

  If (.Not. Allocated(C6ab)) Then
     Allocate(C6ab(natmtot,natmtot), R0_eff_ab(natmtot,natmtot))
     Call get_TS_parameters()
  End If

  e_disp = vdw_energy_pairwiseC6(input%groundstate%TSvdWparameters%s6, input%groundstate%TSvdWparameters%sr6, input%groundstate%TSvdWparameters%d, input%groundstate%TSvdWparameters%cutoff, C6ab, R0_eff_ab)
  If (associated(input%properties)) Then
     If (associated(input%properties%TSvdW)) Then
        Open(246, File="TSvdW.OUT")
        Write(246,'(F18.8)')e_disp
        Close(246)
     End If
  End If
End Subroutine TS_vdW_energy
