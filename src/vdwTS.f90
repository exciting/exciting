Subroutine vdwTS
  !Tkatchenko-Scheffler van der Waals correction 
  Use vdwTS_module, Only:  integrand_numerator, integrand_denominator,&
       current_atom, current_species, sph_int, &
       atoms_in_ambit,&
       list_of_positions_hirshfeld, list_of_species_hirshfeld, num_of_atoms_in_sphere_hirshfeld, &
       get_free_atom_vdw_param, R0_free
  Use mod_atoms, Only: sprmax, atposc, nspecies, natoms, spzn, idxas, natmtot
  Use modinput
  Use modspdeflist
  Implicit None
  Integer :: nsph, nr
  Real(8) :: I_numerator, I_denominator
  Real(8) :: V_ratio
  Real(8) :: max_sprmax
  Real(8) :: C6_free(nspecies), alpha_free_is(nspecies)
  Integer :: is
  Real(8), Allocatable :: C6_eff(:), alpha_free_idxas(:), R0_eff(:)
  Real(8), Allocatable :: xyz(:,:)
  Real(8) :: e_vdwTS
  max_sprmax = maxval(sprmax)
  Allocate(R0_free(nspecies))
  Do is = 1, nspecies
     Call get_free_atom_vdw_param(-spzn(is), C6_free(is), alpha_free_is(is), R0_free(is))
  End Do

  Call init0
  ! read density from file
  Call readstate

  Allocate(C6_eff(natmtot), alpha_free_idxas(natmtot), R0_eff(natmtot))
  Allocate(xyz(3, natmtot))
  nsph=590 !possible numbers are: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3740, 3890, 4334, 4802, 5294, 5810
  nr=120
!!!! quick calc:
!  nsph=146
!  nr=80
!!!!
  Do current_species = 1, nspecies
     Do current_atom = 1,natoms(current_species)
        Call atoms_in_ambit(max_sprmax + sprmax(current_species), atposc(:, current_atom, current_species), list_of_positions_hirshfeld, list_of_species_hirshfeld, num_of_atoms_in_sphere_hirshfeld)
        I_numerator = sph_int(atposc(:, current_atom, current_species), R0_free(current_species), nsph,nr,integrand_numerator)
        I_denominator = sph_int( (/ 0d0, 0d0, 0d0 /), R0_free(current_species), 1, 80, integrand_denominator)
        V_ratio = I_numerator/I_denominator
        C6_eff(idxas(current_atom, current_species)) = V_ratio**2 * C6_free(current_species)
        alpha_free_idxas(idxas(current_atom, current_species)) =  alpha_free_is(current_species)
        R0_eff(idxas(current_atom, current_species)) = V_ratio**(1d0/3) * R0_free(current_species)
        xyz(:, idxas(current_atom, current_species)) = atposc(:, current_atom, current_species)
     End Do ! current_atom
  End Do ! current_species

  Call calc_e_vdwTS()
  Write(*,*)'e_vdwTS = ', e_vdwTS
Contains

  Subroutine calc_e_vdwTS
    Use DFT_D2_subroutines, Only: getlatticerepetition
    Implicit None
    Real(8), Parameter :: s6=1d0, rs6=0.94d0, damping_const=20d0, cutoff=95d0
    Integer :: latrep(3)
    integer :: iat, jat, tau_a, tau_b, tau_c
    real(8) :: C6ab, tau(3), tau_coeff(3), dx, dy, dz, r, r6, damp6
    Call getlatticerepetition(latrep, cutoff)
    e_vdwTS = 0d0
    Do iat = 1,natmtot-1
       Do jat = iat+1,natmtot
          C6ab = 2d0 / ( alpha_free_idxas(jat)/(alpha_free_idxas(iat)*C6_eff(jat) + alpha_free_idxas(iat)/(alpha_free_idxas(jat)*C6_eff(iat) )))
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
                   damp6=1d0/(1d0+Exp(-damping_const*(r/(rs6*(R0_eff(iat)+R0_eff(jat)))-1d0)))
                   r6=r**6
                   e_vdwTS =e_vdwTS+C6ab*damp6/r6
                End Do
             End Do
          End Do
       End Do
    End Do

    Do iat = 1,natmtot
       jat = iat
       C6ab = C6_eff(iat)
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
                damp6=1d0/(1d0+Exp(-damping_const*(r/(rs6*(R0_eff(iat)+R0_eff(jat)))-1d0)))
                r6 = r**6
                e_vdwTS = e_vdwTS+C6ab*damp6/r6*0.50d0
             End Do ! tau_c
          End Do ! tau_b
       End Do ! tau_a
    End Do ! iat
    e_vdwTS = -s6*e_vdwTS

  End Subroutine calc_e_vdwTS
End Subroutine vdwTS
