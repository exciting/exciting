Subroutine DFT_D2_force
  Use DFT_D2_mod, Only: s6, rs6, damping_const, cutoff, loadoldpar
  Use mod_atoms, Only: natmtot
  Use vdw_general_routines, Only: vdw_force_pairwiseC6
  Use mod_force, Only: force_disp
  Implicit None
  Real(8) :: C6ab(natmtot,natmtot), R0ab(natmtot,natmtot)
  Call loadoldpar(C6ab,R0ab)
  force_disp=vdw_force_pairwiseC6(s6, rs6, damping_const, cutoff, C6ab, R0ab)
End Subroutine DFT_D2_force
