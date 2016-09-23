Subroutine TS_vdW_force
  !Tkatchenko-Scheffler van der Waals correction 
  Use TS_vdW_module, Only: get_TS_parameters, C6ab, R0_eff_ab!, s6, sr6, damping_const, cutoff
  Use mod_atoms, Only: natmtot
  Use vdw_general_routines, Only: vdw_force_pairwiseC6
  Use mod_force, Only: force_disp
  Use modinput
  Implicit None

  If (.Not. Allocated(C6ab)) Then
     Allocate(C6ab(natmtot,natmtot), R0_eff_ab(natmtot,natmtot))
     Call get_TS_parameters()
  End If

  force_disp=vdw_force_pairwiseC6(input%groundstate%TSvdWparameters%s6, input%groundstate%TSvdWparameters%sr6, input%groundstate%TSvdWparameters%d, input%groundstate%TSvdWparameters%cutoff, C6ab, R0_eff_ab)
End Subroutine TS_vdW_force
