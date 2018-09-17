Subroutine DFT_D2_energy
  Use mod_energy, Only: e_disp
  Use DFT_D2_module, Only : loadoldpar!, s6, sr6, damping_const, cutoff
  Use vdw_general_routines, Only: vdw_energy_pairwiseC6, set_default_vdW_parameters
  Use mod_atoms, Only: natmtot
  Use modinput
  Implicit None
  Real(8), Allocatable :: C6ab(:,:), R0ab(:,:)

  Call set_default_vdW_parameters
  
  Allocate(C6ab(natmtot, natmtot), R0ab(natmtot, natmtot))

  Call loadoldpar(C6ab,R0ab)

  e_disp = vdw_energy_pairwiseC6(input%groundstate%DFTD2parameters%s6, input%groundstate%DFTD2parameters%sr6, input%groundstate%DFTD2parameters%d, input%groundstate%DFTD2parameters%cutoff, C6ab, R0ab)
  If (associated(input%properties)) Then
     If (associated(input%properties%DFTD2)) Then
        Open(246, File="DFTD2.OUT")
        Write(246,'(F18.8)')e_disp
        Close(246)
     End If
  End If
End Subroutine DFT_D2_energy
