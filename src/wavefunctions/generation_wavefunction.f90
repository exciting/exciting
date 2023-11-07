!> Routines to generate first and second variation wavefunction
Module generation_wavefunction
  Use precision, only: dp 
  Use asserts, only: assert
  Use constants, only: zzero
  Use modinput, only: input, issvlo
  Use modmpi, only: mpiglobal
  Use errors_warnings, only: terminate_if_true
  implicit none
  private
  public :: generate_basisfunction_secondvariation_MT

contains
  !> Generates the second variation basis in the muffin-tin region. The basis can be
  !> generated in two ways:
  !> <ol>
  !> <li> <b>Second variation with local orbital (issvlo()==true).</b> <br>
  !> The total number of second variation basis functions (num_of_basis_funs_sv) is the 
  !> sum of the first variational states (nstfv) and the total number of local
  !> orbitals (nlotot). 
  !> For the first nstfv basis functions the second variation basis is the first variational
  !> wavefunction in terms of a spherical harmonic expansion, only of the APW contribution. 
  !> For atom \(\alpha\) and a particular \(k\)-point \(\mathbf{p}\), the \(r\)-dependent 
  !> \((l,m)\)-coefficient of the \(i\)th basis function where \(1 < i \leq \mathrm{nstfv}\) are given by
  !> \[
  !>   \Phi^{i\mathbf{p}}_{\alpha lm}\left(r\right)=\sum_{\mathbf{G}}b^{i\mathbf{p}}_{\mathbf{G}}
  !>   \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}\left(\mathbf{G+p}\right)u^{\alpha}_{jl}\left(r\right),
  !> \]
  !> where \(b^{i\mathbf{p}}\) is the \(i\)th eigenvector returned from routine
  !> seceqn; \(A^{\alpha}_{jlm}\left(\mathbf{G+p}\right)\) is the matching coefficient;
  !> \(M^{\alpha}_l\) is the order of the APW; \(u^{\alpha}_{jl}\) is the APW radial
  !> function. For the \(i\)th basis where \(\mathrm{nstfv}+1 < i \leq \mathrm{nstfv+nlotot}\), 
  !> the second variation basis is instead defined as the \(j\)th local orbital radial function
  !> \(v_{j}\left(r\right)\), where \(1 < j \leq \mathrm{nlotot}\).
  !> </li>
  !> <li> <b> Traditional second variation. </b><br>
  !> Reference: Li, Chun, <i>et al.</i>, <i>Phys. Rev. B</i> <b>42</b>, 5433, 1990
  !> The total number of second variation basis function (num_of_basis_funs_sv) is equal to the
  !> number of the first variational states (nstfv). The basis function is therefore 
  !> the first variational wavefunction constructed as an expansion of spherical harmonic 
  !> expansion of the APW and local orbitals. For atom \(\alpha\) and a particular 
  !> \(k\)-point \(\mathbf{p}\), the \(r\)-dependent \((l,m)\)-coefficients of the 
  !> wavefunction for the \(i\)th state are given by
  !> \[
  !>   \Phi^{i\mathbf{p}}_{\alpha lm}\left(r\right)=\sum_{\mathbf{G}}b^{i\mathbf{p}}_{\mathbf{G}}
  !>   \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}\left(\mathbf{G+p}\right)u^{\alpha}_{jl}\left(r\right)
  !>   +\sum_{j=1}^{N^{\alpha}}b^{i\mathbf{p}}_{(\alpha,j,m)}v^{\alpha}_{j}\left(r\right) \delta_{\alpha,\alpha_j}
  !>   \delta_{l,l_j} \delta_{m,m_j}, 
  !> \]
  !> where \(N^{\alpha}\) is the number of local-orbitals; \(v^{\alpha}_j\) is the \(j\)th 
  !> local-orbital radial function; and \(\left(\alpha,j,m\right)\) is a compound
  !> index for the location of the local-orbital in the eigenvector. 
  !> </li>
  !> </ol>       
  subroutine generate_basisfunction_secondvariation_MT(lmax, lmmax,&
       & ia, is, ngp, apwalm, evecfv, wfmtfv)
    Use mod_Gkvector, only: ngkmax
    Use mod_atoms, only: natmtot
    Use mod_APW_LO, only: apwordmax, nlotot
    Use mod_muffin_tin, only: nrcmtmax, lmmaxapw
    Use mod_eigenvalue_occupancy, only: nstfv
    Use svlo, only: get_num_of_basis_funs_sv
    Use mod_eigensystem, only: nmatmax
    implicit none
    !> maximum angular momentum for potentials and densities (lmax = input%groundstate%lmaxvr)
    integer, intent(in) :: lmax
    !> \(lmmax = (lmaxvr + 1)^2\)
    integer, intent(in) :: lmmax
    !> Index atom
    integer, intent(in) :: ia
    !> Index species
    integer, intent(in) :: is
    !> Number of G+k vector for the used k-point
    integer, intent(in) :: ngp
    !> APW matching coefficients  
    complex(dp), intent(in) :: apwalm(ngkmax, apwordmax, lmmaxapw, natmtot)
    !> First variation eigenvector
    complex(dp), intent(in) :: evecfv (nmatmax,nstfv)
    !> First variation wavefunction in the MT region
    complex(dp), intent(out) :: wfmtfv(lmmax, nrcmtmax, get_num_of_basis_funs_sv())

    Integer :: num_of_basis_funs_sv
    num_of_basis_funs_sv = get_num_of_basis_funs_sv()
    
    wfmtfv(:,:,:)=zzero 

    If (issvlo()) Then
       call terminate_if_true(mpiglobal,nlotot==0, "SVLO does not work with nlotot = 0")
       call assert(num_of_basis_funs_sv == nstfv + nlotot, "For svlo basis, &
            & expect total number of basis functions to be equal to & 
            & first variational states plus total number of local &
            & orbitals")
       call generate_basisfunction_secondvariation_MT_APW(lmax, &
            & lmmax, ia, is, ngp, apwalm, evecfv(:,1:nstfv), wfmtfv(:,:,1:nstfv))
       call generate_basisfunction_secondvariation_MT_lo(lmmax, &
            & ia, is, wfmtfv(:,:,nstfv+1:num_of_basis_funs_sv))
    Else
       call assert(num_of_basis_funs_sv == nstfv , "For standard sv &
            & basis, expect total number of basis functions to be &
            & equal to first variational states ")
       call generate_wavefunction_firstvariation_MT_APW_and_lo(lmax,&
            & lmmax, ia, is, ngp, apwalm ,evecfv, wfmtfv) 
    End If
  end subroutine generate_basisfunction_secondvariation_MT

  !> Generation of the first variation wavefunction, in standard second variation 
  !> are the second variation basis function  
  subroutine generate_wavefunction_firstvariation_MT_APW_and_lo(lmax,&
       & lmmax, ia, is, ngp, apwalm,  evecfv, wfmtfv) 
    Use mod_Gkvector, only: ngkmax
    Use mod_atoms, only: natmtot
    Use mod_APW_LO, only: apwordmax
    Use mod_muffin_tin, only: nrcmtmax, lmmaxapw
    Use mod_eigenvalue_occupancy, only: nstfv
    Use mod_eigensystem, only: nmatmax
    implicit none
    !> maximum angular momentum for potentials and densities (lmax = input%groundstate%lmaxvr)
    integer, intent(in) :: lmax
    !> \(lmmax = (lmax + 1)^2\)
    integer, intent(in) :: lmmax
    !> Index atom
    integer, intent(in) :: ia
    !> Index species
    integer, intent(in) :: is
    !> Number of G+k vector for the used k-point
    integer, intent(in) :: ngp
    !> APW matching coefficients  
    complex(dp), intent(in) :: apwalm(ngkmax, apwordmax, lmmaxapw, natmtot)
    !> First variation eigenvector
    complex(dp), intent(in) :: evecfv (nmatmax,nstfv)
    !> First variation wavefunction in the MT region
    complex(dp), intent(out) :: wfmtfv(lmmax, nrcmtmax, nstfv)
    !> Internal variables
    integer :: ist

    ! Loop over the number of first variation wavefunction 
    Do ist = 1, nstfv
       Call wavefmt (input%groundstate%lradstep, &
            & lmax, is, ia, ngp, apwalm, &
            & evecfv(:, ist), lmmax, wfmtfv(:, :, ist))
    End Do
  end subroutine generate_wavefunction_firstvariation_MT_APW_and_lo

  !> Generation of second variation basis function from 1 to nstfv. The basis functions
  !> are construct in terms of APW without local orbitals 
  subroutine generate_basisfunction_secondvariation_MT_APW(lmax, &
       & lmmax, ia, is, ngp, apwalm, evecfv, wfmtfv) 
    Use mod_Gkvector, only: ngkmax
    Use mod_atoms, only: natmtot
    Use mod_APW_LO, only: apwordmax
    Use mod_muffin_tin, only: nrcmtmax, lmmaxapw
    Use mod_eigenvalue_occupancy, only: nstfv
    Use mod_eigensystem, only: nmatmax
    implicit none
    !> maximum angular momentum for potentials and densities (lmax = input%groundstate%lmaxvr)
    integer, intent(in) :: lmax
    !> \(lmmax = (lmax + 1)^2\)
    integer, intent(in) :: lmmax
    !> Index atom
    integer, intent(in) :: ia
    !> Index species
    integer, intent(in) :: is
    !> Number of G+k vector for the used k-point
    integer, intent(in) :: ngp
    !> APW matching coefficients  
    complex(dp), intent(in) :: apwalm(ngkmax, apwordmax, lmmaxapw, natmtot)
    !> First variation eigenvector
    complex(dp), intent(in) :: evecfv (nmatmax,nstfv)
    !> First variation wavefunction in the MT region
    complex(dp), intent(out) :: wfmtfv(lmmax, nrcmtmax, nstfv)
    !> Internal variables
    integer :: ist

    ! Loop over the first variation states 
    Do ist = 1, nstfv
       Call wavefmt_apw (input%groundstate%lradstep, &
            & lmax, is, ia, ngp, apwalm, &
            & evecfv(:, ist), lmmax, wfmtfv(:, :, ist))
    End Do

  end subroutine generate_basisfunction_secondvariation_MT_APW

  !> Generation of second variation basis function from nstsv+1 to nstfv+nlotot,
  !> which are the local orbitals 
  subroutine generate_basisfunction_secondvariation_MT_lo(lmmax, &
       & ia, is, wfmtfv) 
    Use mod_Gkvector, only: ngkmax
    Use mod_atoms, only: idxas
    Use mod_APW_LO, only: nlorb, lorbl, lofr, nlotot
    Use mod_muffin_tin, only: nrcmtmax, idxlm, nrmt
    Use mod_eigensystem, only: idxlo
    implicit none
    !> \(lmmax = (lmax + 1)^2\) with lmax, maximum angular momentum for potentials and densities
    integer, intent(in) :: lmmax
    !> Index atom
    integer, intent(in) :: ia
    !> Index species
    integer, intent(in) :: is
    !> First variation wavefunction in the MT region
    complex(dp), intent(out) :: wfmtfv(lmmax, nrcmtmax, nlotot)
    ! Internal variables
    !> Index to atoms and species 
    integer :: ias
    !> Local orbital angular momentum 
    integer :: l
    !> Index to (l,m) pairs      
    integer :: lm 
    !> Position of the local-orbitals in the Hamiltonian and overlap matrices
    integer :: ist      
    !> Number of muffin-tin radial points for the species is
    integer :: nrmt_is
    !> Step size going from the fine to the coarse radial mesh
    integer :: lradstep
    !> Counter over the coarse radial grid points
    integer :: nr
    !> Loop indixes
    integer :: ir, ilo, m 

    ias = idxas(ia,is)
    nrmt_is = nrmt(is) 
    lradstep = input%groundstate%lradstep

    ! Loop over the local orbital of the species is
    Do ilo = 1, nlorb(is) 
       ! Local orbital angular momentum 
       l = lorbl(ilo,is)
       If (l .Le.input%groundstate%lmaxvr ) Then
          Do m = -l,l
             lm = idxlm(l,m)
             ist = idxlo (lm, ilo, ias)
             nr = 0
             ! Loop over the coarse grid.  
             Do ir = 1, nrmt_is, lradstep 
                nr = nr+1 
                wfmtfv(lm,nr,ist) = lofr(ir,1,ilo,ias)
             End Do
          End Do
       End If
    End Do

  end subroutine generate_basisfunction_secondvariation_MT_lo

end module generation_wavefunction
