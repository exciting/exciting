Module vdw_general_routines

Contains

  Function cross(a,b)
    Real(8) :: cross(3)
    Real(8), Intent(in) :: a(3), b(3)

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  End Function cross

  Subroutine getlatticerepetition(latrep, cutoff)
    Use modinput
    Implicit None
    Integer, Intent(out) :: latrep(3)
    Real(8), Intent(in) :: cutoff
    Real(8) :: vec(3,3)
    Integer :: icount

    !orthogonal system
    vec(:,1) = cross(input%structure%crystal%basevect(:,2),input%structure%crystal%basevect(:,3))
    vec(:,2) = cross(input%structure%crystal%basevect(:,1),input%structure%crystal%basevect(:,3))
    vec(:,3) = cross(input%structure%crystal%basevect(:,1),input%structure%crystal%basevect(:,2))
    Do icount = 1,3
       vec(:,icount) = vec(:,icount)/Sqrt(Dot_product(vec(:,icount),vec(:,icount)))!normalize
       latrep(icount) = Int(abs(cutoff/(Dot_product(input%structure%crystal%basevect(:,icount),vec(:,icount))))) + 1
    End Do
  End Subroutine getlatticerepetition

  Real(8) Function vdw_energy_pairwiseC6(s6, sr6, damping_const, cutoff, C6ab, R0ab)
    Use mod_atoms, Only: natmtot, nspecies, natoms, idxas, atposc
    Use modinput
    Implicit None
    Real(8), Intent(in) :: s6, sr6, damping_const, cutoff, C6ab(natmtot, natmtot), R0ab(natmtot, natmtot)
    Real(8) :: xyz(3, natmtot)
    Integer :: latrep(3)
    Integer :: iat, jat, tau_a, tau_b, tau_c, ia, is
    Real(8) :: tau(3), tau_coeff(3), dx, dy, dz, r, r6, damp6
    Do is = 1, nspecies
       Do ia = 1,natoms(is)
          xyz(:, idxas(ia, is)) = atposc(:, ia, is)
       End Do
    End Do

    Call getlatticerepetition(latrep, cutoff)

    vdw_energy_pairwiseC6 = 0d0
    Do iat = 1,natmtot-1
       Do jat = iat+1,natmtot
          Do tau_a = -latrep(1),latrep(1)
             Do tau_b = -latrep(2),latrep(2)
                Do tau_c = -latrep(3),latrep(3)
                   tau_coeff = Dble((/tau_a, tau_b, tau_c/))
                   tau = Matmul(input%structure%crystal%basevect,tau_coeff)
                   dx=xyz(1,iat)-xyz(1,jat)+tau(1)
                   dy=xyz(2,iat)-xyz(2,jat)+tau(2)
                   dz=xyz(3,iat)-xyz(3,jat)+tau(3)
                   r=Sqrt(dx*dx+dy*dy+dz*dz)
                   If(r .Gt. cutoff) Cycle
                   damp6=1d0/(1d0+Exp(-damping_const*(r/(sr6*R0ab(iat, jat))-1d0)))
                   r6=r**6
                   vdw_energy_pairwiseC6 =vdw_energy_pairwiseC6+C6ab(iat, jat)*damp6/r6
                End Do
             End Do
          End Do
       End Do
    End Do

    Do iat = 1,natmtot
       jat = iat
       Do tau_a = -latrep(1),latrep(1)
          Do tau_b = -latrep(2),latrep(2)
             Do tau_c = -latrep(3),latrep(3)
                If (tau_a.Eq.0 .And. tau_b.Eq.0 .And. tau_c.Eq.0) Cycle
                tau_coeff = Dble((/tau_a, tau_b, tau_c/))
                tau = Matmul(input%structure%crystal%basevect,tau_coeff)
                dx = tau(1)
                dy = tau(2)
                dz = tau(3)
                r = Sqrt(dx*dx+dy*dy+dz*dz)
                If(r .Gt. cutoff) Cycle
                damp6=1d0/(1d0+Exp(-damping_const*(r/(sr6*R0ab(iat, jat))-1d0)))
                r6 = r**6
                vdw_energy_pairwiseC6 = vdw_energy_pairwiseC6+C6ab(iat, jat)*damp6/r6*0.50d0
             End Do ! tau_c
          End Do ! tau_b
       End Do ! tau_a
    End Do ! iat
    vdw_energy_pairwiseC6 = -s6*vdw_energy_pairwiseC6
  End Function vdw_energy_pairwiseC6

  Function vdw_force_pairwiseC6(s6, sr6, damping_const, cutoff, C6ab, R0ab)
    Use mod_atoms, Only: natmtot, nspecies, natoms, idxas, atposc
    Use modinput
    Implicit None
    Real(8), Intent(in) :: s6, sr6, damping_const, cutoff, C6ab(natmtot, natmtot), R0ab(natmtot, natmtot)
    Real(8) :: vdw_force_pairwiseC6(3,natmtot)
    Real(8) :: xyz(3, natmtot)
    Integer :: latrep(3)
    Integer :: iat, jat, tau_a, tau_b, tau_c, ia, is
    Real(8) :: tau(3), tau_coeff(3), dx, dy, dz, r, r8, r0ab_, help1_exp, help2
    Do is = 1, nspecies
       Do ia = 1,natoms(is)
          xyz(:, idxas(ia, is)) = atposc(:, ia, is)
       End Do
    End Do

    Call getlatticerepetition(latrep, cutoff)

    vdw_force_pairwiseC6 = 0
    Do iat = 1,natmtot
       Do jat = 1,natmtot
          If (iat .Eq. jat) Cycle
          Do tau_a = -latrep(1),latrep(1)
             Do tau_b = -latrep(2),latrep(2)
                Do tau_c = -latrep(3),latrep(3)
                   tau_coeff = Dble((/tau_a, tau_b, tau_c/))
                   tau = Matmul(input%structure%crystal%basevect,tau_coeff)
                   dx=xyz(1,jat)-xyz(1,iat)+tau(1)
                   dy=xyz(2,jat)-xyz(2,iat)+tau(2)
                   dz=xyz(3,jat)-xyz(3,iat)+tau(3)
                   r=Sqrt(dx*dx+dy*dy+dz*dz)
                   If(r.Gt.cutoff) Cycle
                   r8=r**8
                   r0ab_=sr6*R0ab(iat, jat)
                   help1_exp=damping_const*(r/r0ab_-1.)
                   If(help1_exp .Gt. 100) Then
                      help2=-6*C6ab(iat, jat)/r8
                   Else
                      help1_exp=Exp(help1_exp)
                      help2=C6ab(iat, jat)*help1_exp*(damping_const*r/(1+help1_exp)**2-6*r0ab_/(1+help1_exp))/r8/r0ab_
                   End If
                   vdw_force_pairwiseC6(1,iat) =vdw_force_pairwiseC6(1,iat)+help2*dx
                   vdw_force_pairwiseC6(2,iat) =vdw_force_pairwiseC6(2,iat)+help2*dy
                   vdw_force_pairwiseC6(3,iat) =vdw_force_pairwiseC6(3,iat)+help2*dz
                End Do
             End Do
          End Do
       End Do
    End Do
    vdw_force_pairwiseC6 = -s6*vdw_force_pairwiseC6
  End Function vdw_force_pairwiseC6

  Subroutine set_default_vdW_parameters
    Use modinput
    Use inputdom
    Implicit None
    
    If ( .Not. (associated(input%groundstate%DFTD2parameters))) Then
       ! set the default values if element not present
       input%groundstate%DFTD2parameters => getstructDFTD2parameters (emptynode)
    End If
    If ( .Not. (associated(input%groundstate%TSvdWparameters))) Then
       ! set the default values if element not present
       input%groundstate%TSvdWparameters => getstructTSvdWparameters (emptynode)
    End If
  End Subroutine set_default_vdW_parameters
  
End Module vdw_general_routines
