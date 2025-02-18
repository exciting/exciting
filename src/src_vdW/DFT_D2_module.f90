Module DFT_D2_module
  use precision, only: i32, dp
  Implicit None
  Real(dp), Parameter :: au_to_ang = 0.52917726_dp
  Real(dp), Parameter :: J_to_au = 4.35974417e-18_dp
  Real(dp), Parameter :: N_A = 6.02214129d23 !1/mol
  !  Real(8), Parameter :: damping_const = 20_dp
  Integer(i32), Parameter :: max_elem = 86
  !  Real(8), Parameter :: cutoff = 95 !au
  !  Real(8), Parameter :: s6 = 0.75 !for pbe: s6 = 0.75
  !  Real(8), Parameter :: sr6 = 1.1 !for pbe: sr6 = 1.1

  private
  public  :: loadoldpar

Contains

  Subroutine loadoldpar(C6ab_idxas,R0ab_idxas)
    Use mod_atoms, Only: nspecies, natoms, spzn, idxas, natmtot
    Implicit None
    Real(dp) :: r0(max_elem),c6(max_elem)
    Integer(i32) :: is1, ia1, is2, ia2
    Real(dp) :: C6ab_idxas(natmtot, natmtot), R0ab_idxas(natmtot, natmtot)

    ! the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799 (tab 1)
    ! refer to the following values multiplied by 1.1 (sr6 in this code)

    r0(1:max_elem) = (/0.91_dp , & ! H
         0.92_dp , & ! He
         0.75_dp , & ! Li
         1.28_dp , & ! Be
         1.35_dp , & ! B
         1.32_dp , & ! C
         1.27_dp , & ! N
         1.22_dp , & ! O
         1.17_dp , & ! F
         1.13_dp , & ! Ne
         1.04_dp , & ! Na
         1.24_dp , & ! Mg
         1.49_dp , & ! Al
         1.56_dp , & ! Si
         1.55_dp , & ! P
         1.53_dp , & ! S
         1.49_dp , & ! Cl
         1.45_dp , & ! Ar
         1.35_dp , & ! K
         1.34_dp , & ! Ca
         1.42_dp , & ! Sc
         1.42_dp , & ! Ti
         1.42_dp , & ! V
         1.42_dp , & ! Cr
         1.42_dp , & ! Mn
         1.42_dp , & ! Fe
         1.42_dp , & ! Co
         1.42_dp , & ! Ni
         1.42_dp , & ! Cu
         1.42_dp , & ! Zn
         1.50_dp , & ! Ga
         1.57_dp , & ! Ge
         1.60_dp , & ! As
         1.61_dp , & ! Se
         1.59_dp , & ! Br
         1.57_dp , & ! Kr
         1.48_dp , & ! Rb
         1.46_dp , & ! Sr
         1.49_dp , & ! Y
         1.49_dp , & ! Zr
         1.49_dp , & ! Nb
         1.49_dp , & ! Mo
         1.49_dp , & ! Tc
         1.49_dp , & ! Ru
         1.49_dp , & ! Rh
         1.49_dp , & ! Pd
         1.49_dp , & ! Ag
         1.49_dp , & ! Cd
         1.52_dp , & ! In
         1.64_dp , & ! Sn
         1.71_dp , & ! Sb
         1.72_dp , & ! Te
         1.72_dp , & ! I
         1.71_dp , & ! Xe
         1.638_dp , & ! Cs
         1.602_dp , & ! Ba
         1.564_dp , & ! La
         1.594_dp , & ! Ce
         1.594_dp , & ! Pr
         1.594_dp , & ! Nd
         1.594_dp , & ! Pm
         1.594_dp , & ! Sm
         1.594_dp , & ! Eu
         1.594_dp , & ! Gd
         1.594_dp , & ! Tb
         1.594_dp , & ! Dy
         1.594_dp , & ! Ho
         1.594_dp , & ! Er
         1.594_dp , & ! Tm
         1.594_dp , & ! Yb
         1.594_dp , & ! Lu
         1.625_dp , & ! Hf
         1.611_dp , & ! Ta
         1.611_dp , & ! W
         1.611_dp , & ! Re
         1.611_dp , & ! Os
         1.611_dp , & ! Ir
         1.611_dp , & ! Pt
         1.611_dp , & ! Au
         1.598_dp , & ! Hg
         1.805_dp , & ! Tl
         1.767_dp , & ! Pb
         1.725_dp , & ! Bi
         1.823_dp , & ! Po
         1.810_dp , & ! At
         1.749_dp/)  ! Rn
    ! (/ 0.91_dp,0.92_dp,&
    !         0.75_dp,1.28_dp,1.35_dp,1.32_dp,1.27_dp,1.22_dp,1.17_dp,1.13_dp, &
    !         1.04_dp,1.24_dp,1.49_dp,1.56_dp,1.55_dp,1.53_dp,1.49_dp,1.45_dp,&
    !         1.35_dp,1.34_dp,&
    !         1.42_dp,1.42_dp,1.42_dp,1.42_dp,1.42_dp,&
    !         1.42_dp,1.42_dp,1.42_dp,1.42_dp,1.42_dp,&
    !         1.50_dp,1.57_dp,1.60_dp,1.61_dp,1.59_dp,1.57_dp,&
    !         1.48_dp,1.46_dp,&
    !         1.49_dp,1.49_dp,1.49_dp,1.49_dp,1.49_dp,&
    !         1.49_dp,1.49_dp,1.49_dp,1.49_dp,1.49_dp,&
    !         1.52_dp,1.64_dp,1.71_dp,1.72_dp,1.72_dp,1.71_dp,&
    !         1.638_dp,1.602_dp,1.564_dp,1.594_dp,1.594_dp,1.594_dp,1.594_dp,&
    !         1.594_dp,1.594_dp,1.594_dp,1.594_dp,1.594_dp,1.594_dp,1.594_dp,&
    !         1.594_dp,1.594_dp,1.594_dp,&
    !         1.625_dp,1.611_dp,1.611_dp,1.611_dp,1.611_dp,1.611_dp,1.611_dp,&
    !         1.611_dp,&
    !         1.598_dp,1.805_dp,1.767_dp,1.725_dp,1.823_dp,1.810_dp,1.749_dp/)


    c6(1:max_elem) = (/0.14_dp , & ! H
         0.08_dp , & ! He
         1.61_dp , & ! Li
         1.61_dp , & ! Be
         3.13_dp , & ! B
         1.75_dp , & ! C
         1.23_dp , & ! N
         0.70_dp , & ! O
         0.75_dp , & ! F
         0.63_dp , & ! Ne
         5.71_dp , & ! Na
         5.71_dp , & ! Mg
         10.79_dp , & ! Al
         9.23_dp , & ! Si
         7.84_dp , & ! P
         5.57_dp , & ! S
         5.07_dp , & ! Cl
         4.61_dp , & ! Ar
         10.8_dp , & ! K
         10.8_dp , & ! Ca
         10.8_dp , & ! Sc
         10.8_dp , & ! Ti
         10.8_dp , & ! V
         10.8_dp , & ! Cr
         10.8_dp , & ! Mn
         10.8_dp , & ! Fe
         10.8_dp , & ! Co
         10.8_dp , & ! Ni
         10.8_dp , & ! Cu
         10.8_dp , & ! Zn
         16.99_dp , & ! Ga
         17.10_dp , & ! Ge
         16.37_dp , & ! As
         12.64_dp , & ! Se
         12.47_dp , & ! Br
         12.01_dp , & ! Kr
         24.67_dp , & ! Rb
         24.67_dp , & ! Sr
         24.67_dp , & ! Y
         24.67_dp , & ! Zr
         24.67_dp , & ! Nb
         24.67_dp , & ! Mo
         24.67_dp , & ! Tc
         24.67_dp , & ! Ru
         24.67_dp , & ! Rh
         24.67_dp , & ! Pd
         24.67_dp , & ! Ag
         24.67_dp , & ! Cd
         37.32_dp , & ! In
         38.71_dp , & ! Sn
         38.44_dp , & ! Sb
         31.74_dp , & ! Te
         31.50_dp , & ! I
         29.99_dp , & ! Xe
         315.275_dp , & ! Cs
         226.994_dp , & ! Ba
         176.252_dp , & ! La
         140.68_dp , & ! Ce
         140.68_dp , & ! Pr
         140.68_dp , & ! Nd
         140.68_dp , & ! Pm
         140.68_dp , & ! Sm
         140.68_dp , & ! Eu
         140.68_dp , & ! Gd
         140.68_dp , & ! Tb
         140.68_dp , & ! Dy
         140.68_dp , & ! Ho
         140.68_dp , & ! Er
         140.68_dp , & ! Tm
         140.68_dp , & ! Yb
         140.68_dp , & ! Lu
         105.112_dp , & ! Hf
         81.24_dp , & ! Ta
         81.24_dp , & ! W
         81.24_dp , & ! Re
         81.24_dp , & ! Os
         81.24_dp , & ! Ir
         81.24_dp , & ! Pt
         81.24_dp , & ! Au
         57.364_dp , & ! Hg
         57.254_dp , & ! Tl
         63.162_dp , & ! Pb
         63.540_dp , & ! Bi
         55.283_dp , & ! Po
         57.171_dp , & ! At
         56.64_dp/) ! Rn

!!$    c6(1:max_elem) = (/0.14_dp,0.08_dp,&
!!$         1.61_dp,1.61_dp,3.13_dp,1.75_dp,1.23_dp,0.70_dp,0.75_dp,0.63_dp,&
!!$         5.71_dp,5.71_dp,10.79_dp,9.23_dp,7.84_dp,5.57_dp,5.07_dp,4.61_dp,&
!!$         10.8_dp,10.8_dp,10.8_dp,10.8_dp,10.8_dp,&
!!$         10.8_dp,10.8_dp,10.8_dp,10.8_dp,10.8_dp,10.8_dp,10.8_dp,16.99_dp,&
!!$         17.10_dp,16.37_dp,12.64_dp,12.47_dp,12.01_dp,24.67_dp,24.67_dp,&
!!$         24.67_dp,24.67_dp,24.67_dp,24.67_dp,24.67_dp,24.67_dp,24.67_dp,&
!!$         24.67_dp,24.67_dp,24.67_dp,37.32_dp,38.71_dp,38.44_dp,31.74_dp,&
!!$         31.50_dp,29.99_dp,315.275_dp,226.994_dp,176.252_dp,&
!!$         140.68_dp,140.68_dp,140.68_dp,140.68_dp,140.68_dp,140.68_dp,140.68_dp,&
!!$         140.68_dp,140.68_dp,140.68_dp,140.68_dp,140.68_dp,140.68_dp,140.68_dp,&
!!$         105.112_dp,&
!!$         81.24_dp,81.24_dp,81.24_dp,81.24_dp,81.24_dp,81.24_dp,81.24_dp,&
!!$         57.364_dp,57.254_dp,63.162_dp,63.540_dp,55.283_dp,57.171_dp,56.64_dp /)

    Do is1 = 1, nspecies
       Do ia1 = 1,natoms(is1)
          Do is2 = 1, nspecies
             Do ia2 = 1,natoms(is2)
                C6ab_idxas(idxas(ia1,is1), idxas(ia2,is2))=Sqrt(c6(-nint(spzn(is1)))*c6(-nint(spzn(is2))))
                R0ab_idxas(idxas(ia1,is1), idxas(ia2,is2))=r0(-nint(spzn(is1))) + r0(-nint(spzn(is2)))
             End Do
          End Do
       End Do
    End Do

    !convert to au
    C6ab_idxas = C6ab_idxas * 1d6/J_to_au/(au_to_ang**6)/N_A
    R0ab_idxas = R0ab_idxas/au_to_ang
  End Subroutine loadoldpar
End Module DFT_D2_module

