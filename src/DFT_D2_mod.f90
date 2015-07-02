Module DFT_D2_mod
  Implicit None
  Real(8), Parameter :: au_to_ang = 0.52917726d0
  Real(8), Parameter :: J_to_au = 4.35974417d-18
  Real(8), Parameter :: N_A = 6.02214129d23 !1/mol
!  Real(8), Parameter :: damping_const = 20d0
  Integer, Parameter :: max_elem = 86
!  Real(8), Parameter :: cutoff = 95 !au
!  Real(8), Parameter :: s6 = 0.75 !for pbe: s6 = 0.75
!  Real(8), Parameter :: rs6 = 1.1 !for pbe: rs6 = 1.1

Contains

  Subroutine loadoldpar(C6ab_idxas,R0ab_idxas)
    Use mod_atoms, Only: nspecies, natoms, spzn, idxas, natmtot
    Implicit None
    Real(8) :: r0(max_elem),c6(max_elem)
    Integer :: is1, ia1, is2, ia2
    Real(8) :: C6ab_idxas(natmtot, natmtot), R0ab_idxas(natmtot, natmtot)

    ! the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799 (tab 1)
    ! refer to the following values multiplied by 1.1 (rs6 in this code)

    r0(1:max_elem) = (/ 0.91d0,0.92d0,&
         0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, &
         1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,&
         1.35d0,1.34d0,&
         1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
         1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
         1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,&
         1.48d0,1.46d0,&
         1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
         1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
         1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0,&
         1.638d0,1.602d0,1.564d0,1.594d0,1.594d0,1.594d0,1.594d0,&
         1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,&
         1.594d0,1.594d0,1.594d0,&
         1.625d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,&
         1.611d0,&
         1.598d0,1.805d0,1.767d0,1.725d0,1.823d0,1.810d0,1.749d0/)

    c6(1:max_elem) = (/0.14d0,0.08d0,&
         1.61d0,1.61d0,3.13d0,1.75d0,1.23d0,0.70d0,0.75d0,0.63d0,&
         5.71d0,5.71d0,10.79d0,9.23d0,7.84d0,5.57d0,5.07d0,4.61d0,&
         10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,&
         10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,16.99d0,&
         17.10d0,16.37d0,12.64d0,12.47d0,12.01d0,24.67d0,24.67d0,&
         24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,&
         24.67d0,24.67d0,24.67d0,37.32d0,38.71d0,38.44d0,31.74d0,&
         31.50d0,29.99d0,315.275d0,226.994d0,176.252d0,&
         140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,&
         140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,&
         105.112d0,&
         81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,&
         57.364d0,57.254d0,63.162d0,63.540d0,55.283d0,57.171d0,56.64d0 /)

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
End Module DFT_D2_mod
