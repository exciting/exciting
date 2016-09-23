Subroutine DFT_D2_force
  Use DFT_D2_module, Only: loadoldpar!, s6, sr6, damping_const, cutoff
  Use mod_atoms, Only: natmtot
  Use vdw_general_routines, Only: vdw_force_pairwiseC6
  Use mod_force, Only: force_disp
  Use modinput
  Implicit None
  Real(8) :: C6ab(natmtot,natmtot), R0ab(natmtot,natmtot)
  Call loadoldpar(C6ab,R0ab)
  force_disp=vdw_force_pairwiseC6(input%groundstate%DFTD2parameters%s6, input%groundstate%DFTD2parameters%sr6, input%groundstate%DFTD2parameters%d, input%groundstate%DFTD2parameters%cutoff, C6ab, R0ab)
End Subroutine DFT_D2_force
