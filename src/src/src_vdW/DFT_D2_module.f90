Module DFT_D2_module
  Implicit None
  Real(8), Parameter :: au_to_ang = 0.52917726d0
  Real(8), Parameter :: J_to_au = 4.35974417d-18
  Real(8), Parameter :: N_A = 6.02214129d23 !1/mol
  !  Real(8), Parameter :: damping_const = 20d0
  Integer, Parameter :: max_elem = 86
  !  Real(8), Parameter :: cutoff = 95 !au
  !  Real(8), Parameter :: s6 = 0.75 !for pbe: s6 = 0.75
  !  Real(8), Parameter :: sr6 = 1.1 !for pbe: sr6 = 1.1

Contains

  Subroutine loadoldpar(C6ab_idxas,R0ab_idxas)
    Use mod_atoms, Only: nspecies, natoms, spzn, idxas, natmtot
    Implicit None
    Real(8) :: r0(max_elem),c6(max_elem)
    Integer :: is1, ia1, is2, ia2
    Real(8) :: C6ab_idxas(natmtot, natmtot), R0ab_idxas(natmtot, natmtot)

    ! the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799 (tab 1)
    ! refer to the following values multiplied by 1.1 (sr6 in this code)

    r0(1:max_elem) = (/0.91d0 , & ! H
         0.92d0 , & ! He
         0.75d0 , & ! Li
         1.28d0 , & ! Be
         1.35d0 , & ! B
         1.32d0 , & ! C
         1.27d0 , & ! N
         1.22d0 , & ! O
         1.17d0 , & ! F
         1.13d0 , & ! Ne
         1.04d0 , & ! Na
         1.24d0 , & ! Mg
         1.49d0 , & ! Al
         1.56d0 , & ! Si
         1.55d0 , & ! P
         1.53d0 , & ! S
         1.49d0 , & ! Cl
         1.45d0 , & ! Ar
         1.35d0 , & ! K
         1.34d0 , & ! Ca
         1.42d0 , & ! Sc
         1.42d0 , & ! Ti
         1.42d0 , & ! V
         1.42d0 , & ! Cr
         1.42d0 , & ! Mn
         1.42d0 , & ! Fe
         1.42d0 , & ! Co
         1.42d0 , & ! Ni
         1.42d0 , & ! Cu
         1.42d0 , & ! Zn
         1.50d0 , & ! Ga
         1.57d0 , & ! Ge
         1.60d0 , & ! As
         1.61d0 , & ! Se
         1.59d0 , & ! Br
         1.57d0 , & ! Kr
         1.48d0 , & ! Rb
         1.46d0 , & ! Sr
         1.49d0 , & ! Y
         1.49d0 , & ! Zr
         1.49d0 , & ! Nb
         1.49d0 , & ! Mo
         1.49d0 , & ! Tc
         1.49d0 , & ! Ru
         1.49d0 , & ! Rh
         1.49d0 , & ! Pd
         1.49d0 , & ! Ag
         1.49d0 , & ! Cd
         1.52d0 , & ! In
         1.64d0 , & ! Sn
         1.71d0 , & ! Sb
         1.72d0 , & ! Te
         1.72d0 , & ! I
         1.71d0 , & ! Xe
         1.638d0 , & ! Cs
         1.602d0 , & ! Ba
         1.564d0 , & ! La
         1.594d0 , & ! Ce
         1.594d0 , & ! Pr
         1.594d0 , & ! Nd
         1.594d0 , & ! Pm
         1.594d0 , & ! Sm
         1.594d0 , & ! Eu
         1.594d0 , & ! Gd
         1.594d0 , & ! Tb
         1.594d0 , & ! Dy
         1.594d0 , & ! Ho
         1.594d0 , & ! Er
         1.594d0 , & ! Tm
         1.594d0 , & ! Yb
         1.594d0 , & ! Lu
         1.625d0 , & ! Hf
         1.611d0 , & ! Ta
         1.611d0 , & ! W
         1.611d0 , & ! Re
         1.611d0 , & ! Os
         1.611d0 , & ! Ir
         1.611d0 , & ! Pt
         1.611d0 , & ! Au
         1.598d0 , & ! Hg
         1.805d0 , & ! Tl
         1.767d0 , & ! Pb
         1.725d0 , & ! Bi
         1.823d0 , & ! Po
         1.810d0 , & ! At
         1.749d0/)  ! Rn
    ! (/ 0.91d0,0.92d0,&
    !         0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, &
    !         1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,&
    !         1.35d0,1.34d0,&
    !         1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
    !         1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
    !         1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,&
    !         1.48d0,1.46d0,&
    !         1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
    !         1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
    !         1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0,&
    !         1.638d0,1.602d0,1.564d0,1.594d0,1.594d0,1.594d0,1.594d0,&
    !         1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,&
    !         1.594d0,1.594d0,1.594d0,&
    !         1.625d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,&
    !         1.611d0,&
    !         1.598d0,1.805d0,1.767d0,1.725d0,1.823d0,1.810d0,1.749d0/)


    c6(1:max_elem) = (/0.14d0 , & ! H
         0.08d0 , & ! He
         1.61d0 , & ! Li
         1.61d0 , & ! Be
         3.13d0 , & ! B
         1.75d0 , & ! C
         1.23d0 , & ! N
         0.70d0 , & ! O
         0.75d0 , & ! F
         0.63d0 , & ! Ne
         5.71d0 , & ! Na
         5.71d0 , & ! Mg
         10.79d0 , & ! Al
         9.23d0 , & ! Si
         7.84d0 , & ! P
         5.57d0 , & ! S
         5.07d0 , & ! Cl
         4.61d0 , & ! Ar
         10.8d0 , & ! K
         10.8d0 , & ! Ca
         10.8d0 , & ! Sc
         10.8d0 , & ! Ti
         10.8d0 , & ! V
         10.8d0 , & ! Cr
         10.8d0 , & ! Mn
         10.8d0 , & ! Fe
         10.8d0 , & ! Co
         10.8d0 , & ! Ni
         10.8d0 , & ! Cu
         10.8d0 , & ! Zn
         16.99d0 , & ! Ga
         17.10d0 , & ! Ge
         16.37d0 , & ! As
         12.64d0 , & ! Se
         12.47d0 , & ! Br
         12.01d0 , & ! Kr
         24.67d0 , & ! Rb
         24.67d0 , & ! Sr
         24.67d0 , & ! Y
         24.67d0 , & ! Zr
         24.67d0 , & ! Nb
         24.67d0 , & ! Mo
         24.67d0 , & ! Tc
         24.67d0 , & ! Ru
         24.67d0 , & ! Rh
         24.67d0 , & ! Pd
         24.67d0 , & ! Ag
         24.67d0 , & ! Cd
         37.32d0 , & ! In
         38.71d0 , & ! Sn
         38.44d0 , & ! Sb
         31.74d0 , & ! Te
         31.50d0 , & ! I
         29.99d0 , & ! Xe
         315.275d0 , & ! Cs
         226.994d0 , & ! Ba
         176.252d0 , & ! La
         140.68d0 , & ! Ce
         140.68d0 , & ! Pr
         140.68d0 , & ! Nd
         140.68d0 , & ! Pm
         140.68d0 , & ! Sm
         140.68d0 , & ! Eu
         140.68d0 , & ! Gd
         140.68d0 , & ! Tb
         140.68d0 , & ! Dy
         140.68d0 , & ! Ho
         140.68d0 , & ! Er
         140.68d0 , & ! Tm
         140.68d0 , & ! Yb
         140.68d0 , & ! Lu
         105.112d0 , & ! Hf
         81.24d0 , & ! Ta
         81.24d0 , & ! W
         81.24d0 , & ! Re
         81.24d0 , & ! Os
         81.24d0 , & ! Ir
         81.24d0 , & ! Pt
         81.24d0 , & ! Au
         57.364d0 , & ! Hg
         57.254d0 , & ! Tl
         63.162d0 , & ! Pb
         63.540d0 , & ! Bi
         55.283d0 , & ! Po
         57.171d0 , & ! At
         56.64d0/) ! Rn

!!$    c6(1:max_elem) = (/0.14d0,0.08d0,&
!!$         1.61d0,1.61d0,3.13d0,1.75d0,1.23d0,0.70d0,0.75d0,0.63d0,&
!!$         5.71d0,5.71d0,10.79d0,9.23d0,7.84d0,5.57d0,5.07d0,4.61d0,&
!!$         10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,&
!!$         10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,16.99d0,&
!!$         17.10d0,16.37d0,12.64d0,12.47d0,12.01d0,24.67d0,24.67d0,&
!!$         24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,&
!!$         24.67d0,24.67d0,24.67d0,37.32d0,38.71d0,38.44d0,31.74d0,&
!!$         31.50d0,29.99d0,315.275d0,226.994d0,176.252d0,&
!!$         140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,&
!!$         140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,&
!!$         105.112d0,&
!!$         81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,&
!!$         57.364d0,57.254d0,63.162d0,63.540d0,55.283d0,57.171d0,56.64d0 /)

    Do is1 = 1, nspecies
       Do ia1 = 1,natoms(is1)
          Do is2 = 1, nspecies
             Do ia2 = 1,natoms(is2)
                C6ab_idxas(idxas(ia1,is1), idxas(ia2,is2))=Sqrt(c6(-spzn(is1))*c6(-spzn(is2)))
                R0ab_idxas(idxas(ia1,is1), idxas(ia2,is2))=r0(-spzn(is1)) + r0(-spzn(is2))
             End Do
          End Do
       End Do
    End Do

    !convert to au
    C6ab_idxas = C6ab_idxas * 1d6/J_to_au/(au_to_ang**6)/N_A
    R0ab_idxas = R0ab_idxas/au_to_ang
  End Subroutine loadoldpar
End Module DFT_D2_module
