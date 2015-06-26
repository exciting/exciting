Subroutine DFT_D2_energy
  Use mod_energy, Only: e_disp
  Use DFT_D2_mod, Only : s6, rs6, damping_const, cutoff, loadoldpar
  Use vdw_general_routines, Only: vdw_energy_pairwiseC6
  Use mod_atoms, Only: natmtot
  Use modinput
  Implicit None
  Real(8), Allocatable :: C6ab(:,:), R0ab(:,:)

  If ( input%groundstate%do .Eq. "skip" ) Then
     Call init0
  End If
  Allocate(C6ab(natmtot, natmtot), R0ab(natmtot, natmtot))

  Call loadoldpar(C6ab,R0ab)

  e_disp = vdw_energy_pairwiseC6(s6, rs6, damping_const, cutoff, C6ab, R0ab)
  If (associated(input%properties)) Then
     If (associated(input%properties%DFT_D2)) Then
        Open(246, File="E_DFT_D2.OUT")
        Write(246,'(F18.8)')e_disp
        Close(246)
     End If
  End If
End Subroutine DFT_D2_energy
