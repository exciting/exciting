!> This module contains procedures for the calculation of the 
!> Hamiltonian and overlap matrix response upon a phonon-like
!> perturbation and the solution of the corresponding Sternheimer
!> equation.
module phonons_eigensystem
  use dfpt_variables
  use dfpt_eigensystem
  use phonons_variables

  use precision, only: dp

  implicit none
  private

  public :: ph_eig_gen_dSHmat
  public :: ph_eig_sternheimer
  public :: gen_dapwalm, gen_dcfun_ig

  contains

    !> This subroutine computes the response of the overlap and Hamiltonian matrix 
    !> for the given \({\bf k}\) point.
    !>
    !> The total response of the overlap matrix is given by
    !> \[ \delta^{\bf q}_{I \mu} S_{mn}({\bf k}) = \sum_{\kappa,\alpha} p^{I \mu}_{\kappa\alpha}({\bf q})\,
    !>    \delta^{\bf q}_{\kappa\alpha} S_{mn}({\bf k}) \;, \]
    !> with
    !> \[ \delta^{\bf q}_{\kappa\alpha} S_{mn}({\bf k})
    !>    = \int\limits_{{\rm MT} \kappa} {\rm d}{\bf r} \left[ 
    !>      \left( \breve{\bar{\delta}}\phantom{}^{\bf -q}_{\kappa\alpha} \psi_{m{\bf k+q}}({\bf r}) \right)^\ast \,
    !>      \psi_{n{\bf k}}({\bf r}) + 
    !>      \psi_{m{\bf k+q}}^\ast({\bf r}) \,
    !>      \breve{\bar{\delta}}\phantom{}^{\bf q}_{\kappa\alpha} \psi_{n{\bf k}}({\bf r}) \right]
    !>    - \oint\limits_{\partial {\rm MT}\kappa} \left[ \psi_{m{\bf k+q}}^\ast({\bf r})\,
    !>      \psi_{n{\bf k}}({\bf r}) \right]_{\rm IR} \hat{e}_\alpha\, {\rm d}S \;, \]
    !> where \(\breve{\bar{\delta}}\phantom{}^{\bf q}_{\kappa\alpha} \psi_{n{\bf k}}\)
    !> is the contribution to the wavefunction response coming from the soft part of 
    !> the basis function response (i.e., without the gradient).
    !>
    !> Similarly, the full Hamiltonian matrix response is described in terms of the Canonical
    !> response
    !> \[ \delta^{\bf q}_{\kappa\alpha} H_{mn}({\bf k})
    !>    = \int\limits_{{\rm MT} \kappa} {\rm d}{\bf r} \left[ 
    !>      \left( \breve{\bar{\delta}}\phantom{}^{\bf -q}_{\kappa\alpha} \psi_{m{\bf k+q}}({\bf r}) \right)^\ast \,
    !>      \hat{\bf H}\, \psi_{n{\bf k}}({\bf r}) + 
    !>      \psi_{m{\bf k+q}}^\ast({\bf r})\, \hat{\bf H}\,
    !>      \breve{\bar{\delta}}\phantom{}^{\bf q}_{\kappa\alpha} \psi_{n{\bf k}}({\bf r}) \right]
    !>    - \oint\limits_{\partial {\rm MT}\kappa} \left[ \psi_{m{\bf k+q}}^\ast({\bf r})\, \hat{\bf H}\, 
    !>      \psi_{n{\bf k}}({\bf r}) \right]_{\rm IR} \hat{e}_\alpha\, {\rm d}S \\
    !>    + \int\limits_\Omega {\rm d}{\bf r}\, \psi_{m{\bf k+q}}^\ast({\bf r})\, 
    !>      \breve{\delta}\phantom{}^{\bf q}_{\kappa\alpha} V_{\rm eff}({\bf r})\,
    !>      \psi_{n{\bf k}}({\bf r}) \;, \]
    !> where \(\breve{\delta}\phantom{}^{\bf q}_{\kappa\alpha} V_{\rm eff}\) is the soft part
    !> of the effective potential response (i.e., without the gradient). Accordingly, the third
    !> term (also called \(\delta^{\bf q}_{\kappa\alpha} {\bf H}^0\)) is the only part of the
    !> Hamiltonian response that changes during the self-consistency cycle. Hence, the 
    !> first two terms form the so called constant part of the Hamiltonoan matrix response.
    !>
    !> If the displacement pattern `pat` is given, then this subroutine computes the
    !> overlap matrix response \(\delta^{\bf q}_{I \mu} {\bf S}\) and the
    !> constant part of the Hamiltonian matrix response \(\delta^{\bf q}_{I \mu} {\bf H}\).
    !> These are given by two contributions. One coming from muffin-tin integrals 
    !> (see [[gen_dSH0_mt(subroutine)]]) and one coming from surface integrals
    !> (see [[gen_dSH0_ir(subroutine)]]).
    !>
    !> If the effective potential response is given in form of the radial muffin-tin
    !> integrals times Gaunt coefficients `dHmat_mt_basis` and the product with the 
    !> characteristic function in reciprocal space `dpot_cfun_ig`, then this subroutine
    !> computes the corresponding contribution to the Hamiltonian matrix response, 
    !> \(\delta^{\bf q}_{I \mu} {\bf H}^0\) (see [[dfpt_eig_gen_dHmat(subroutine)]]).
    !> 
    !> @note The result is added to the input matrix! @endnote
    subroutine ph_eig_gen_dSHmat( ik, Gkset, Gkqset, fst, lst, eveck, eveckq, apwalmk, apwalmkq, dSmat, dHmat, &
        pat, dHmat_mt_basis, dpot_cfun_ig, dkin_cfun_ig )
      use mod_kpointset, only: Gk_set
      use mod_atoms, only: natmtot
      use mod_APW_LO, only: nlotot
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k}\) and \({\bf G+k+q}\) vectors
      type(Gk_set), intent(in) :: Gkset, Gkqset
      !> first and last state for which the matrix elements are calculated
      integer, intent(in) :: fst, lst
      !> eigenvectors at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: eveck(:,:), eveckq(:,:)
      !> (L)APW matching coefficients \(A^\kappa_{{\bf G+p},lm,\xi}\) at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)
      !> overlap response and constant part of Hamiltonian response
      complex(dp), intent(inout) :: dSmat(:,:), dHmat(:,:)
      !> displacement pattern \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), optional, intent(in) :: pat(3, natmtot)
      !> radial integrals of effective potential response times Gaunt coefficients
      complex(dp), optional, intent(in) :: dHmat_mt_basis(:,:,:)
      !> interstitial effective potential response times characteristic function in reciprocal space
      complex(dp), optional, intent(in) :: dpot_cfun_ig(:)
      !> interstitial (scalar relativistic) kinetic energy response times characteristic function in reciprocal space
      complex(dp), optional, intent(in) :: dkin_cfun_ig(:)

      integer :: nmatkq

      nmatkq = Gkqset%ngk(1, ik) + nlotot

      if( present( pat ) ) then
        ! contribution from muffin-tin integrals
        call gen_dSH0_mt( ik, Gkset, Gkqset, 1, nmatkq, fst, lst, &
               eveck, eveckq, apwalmk, apwalmkq, pat, dSmat, dHmat )
        ! contribution from interstitial integrals
        call gen_dSH0_ir( ik, Gkset, Gkqset, 1, nmatkq, fst, lst, &
               eveck, eveckq, pat, dSmat, dHmat )
      end if

      if( present( dHmat_mt_basis ) .and. present( dpot_cfun_ig ) .and. present( dkin_cfun_ig ) ) then
        ! contribution from potential response
        call dfpt_eig_gen_dHmat( ik, Gkqset, Gkset, 1, nmatkq, fst, lst, &
               eveckq, eveck, apwalmkq, apwalmk, dHmat_mt_basis, dpot_cfun_ig, dkin_cfun_ig, dHmat, &
               Gset=ph_Gqset )
      end if
    end subroutine ph_eig_gen_dSHmat

    !> This subroutine solves the Sternheimer equation for a phonon-like perturbation
    !> for the eigenvalue and eigenvector response at a given \({\bf k}\) point.
    !> 
    !> In general, Sternheimer's equation for the state \(\psi_{n{\bf k}}\) reads
    !> \[ \left[ \hat{\bf H} - \epsilon_{n{\bf k}} \right] \delta \psi_{n{\bf k}}({\bf r})
    !>    = -\left[ \delta \hat{\bf H} - \delta \epsilon_{n{\bf k}} \right] \psi_{n{\bf k}}({\bf r}) \;. \]
    !> In the case of a phonon-like perturbation the wavefunction response consists
    !> of two parts. One coming from the response of the expansion coefficients and the other
    !> coming from the response of the basis functions
    !> \[ \delta^{\bf q}_{I \mu} \psi_{n{\bf k}}({\bf r})
    !>    = \sum_\nu \delta^{\bf q}_{I \mu} C_{\nu n}({\bf k}) \, \phi_{\nu{\bf k}}({\bf r})
    !>      + \sum_\nu C_{\nu n}({\bf k}) \, \delta^{\bf q}_{I \mu} \phi_{\nu{\bf k}}({\bf r}) \;. \]
    !> We want to solve Sternheimer's equation for the response of the expansion coefficients 
    !> (eigenvectors) and the Kohn-Sham energies (eigenvalues). The latter is zero unless
    !> \({\bf q} = {\bf 0}\). To this extend, we expand the eigenvector response in the eigenvectors
    !> of the unperturbed system at wavevector \({\bf k+q}\)
    !> \[ \delta^{\bf q}_{I \mu} C_{\nu n}({\bf k}) 
    !>    = \sum_m \delta^{\bf q}_{I \mu} X_{mn}({\bf k})\, C_{\nu m}({\bf k+q}) \;. \]
    !> This yields
    !> \[ \delta^{\bf q}_{I \mu} X_{mn}({\bf k}) = - \begin{cases}
    !>    \frac{\delta^{\bf q}_{I \mu} H_{mn}({\bf k}) - \epsilon_{n{\bf k}} \delta^{\bf q}_{I \mu} S_{mn}({\bf k})}
    !>    {\epsilon_{m{\bf k+q}} - \epsilon_{n{\bf k}}} & \text{if } \epsilon_{m{\bf k+q}} \neq \epsilon_{n{\bf k}} \\
    !>    \frac{1}{2} \delta^{\bf q}_{I \mu} S_{mn}({\bf k}) & \text{if } \epsilon_{m{\bf k+q}} = \epsilon_{n{\bf k}} 
    !>    \end{cases} \;, \]
    !> where \(\delta^{\bf q}_{I \mu} {\bf S}({\bf k})\) and \(\delta^{\bf q}_{I \mu} {\bf H}({\bf k})\)
    !> are the response of the overlap and Hamiltonian matrix, respectively, as obtained from
    !> [[ph_eig_gen_dSHmat(subroutine)]].
    !>
    !> If `projector=.true.`, then the eigenvector response corresponding to the projection of the
    !> wavefunction response onto the manifold of unoccupied states, 
    !> \(\hat{\bf P}_u \delta^{\bf q}_{I \mu} \psi_{n{\bf k}}({\bf r})\), is computed. This is 
    !> sufficient to compute the density response and ensures (given that \(n\) corresponds to
    !> an occupied state) that the energie differences \(\epsilon_{m{\bf k+q}} - \epsilon_{n{\bf k}}\)
    !> are always greater than zero.
    subroutine ph_eig_sternheimer( ik, Gkqset, fst, lst, evalk, occk, evalkq, occkq, eveckq, dSmat, dHmat, gamma, deval, devec, &
        projector, eps_deg )
      use constants, only: zzero, zone
      use mod_kpointset, only: Gk_set
      use mod_APW_LO, only: nlotot
      use mod_eigenvalue_occupancy, only: occmax
      use modinput
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k+q}\) vectors
      type(Gk_set), intent(in) :: Gkqset
      !> first and last state for which the Sternheimer equation is solved
      integer, intent(in) :: fst, lst
      !> eigenvalues at \({\bf k}\) and \({\bf k+q}\)
      real(dp), intent(in) :: evalk(:), evalkq(:)
      !> occupation numbers at \({\bf k}\) and \({\bf k+q}\)
      real(dp), intent(in) :: occk(:), occkq(:)
      !> eigenvectors at \({\bf k+q}\)
      complex(dp), intent(in) :: eveckq(:,:)
      !> overlap and Hamiltonian matrix response
      complex(dp), intent(in) :: dSmat(:,:), dHmat(:,:)
      !> true for Gamma point phonons
      logical, intent(in) :: gamma
      !> eigenvalue response at \({\bf k}\)
      real(dp), intent(out) :: deval(:)
      !> eigenvector response at \({\bf k}\)
      complex(dp), intent(out) :: devec(:,:)
      !> calculate the eigenvector response corresponding to the projection
      !> of the wavefunction response onto the manifold of unoccupied states
      !> (default: `.true.`)
      logical, optional, intent(in) :: projector
      !> tolerance for degenerate eigenvalues (default: `1e-5`)
      real(dp), optional, intent(in) :: eps_deg

      integer :: nmatkq, ist, jst, nst
      real(dp) :: eps, dev, sig, t1, t2
      logical :: proj

      complex(dp), allocatable :: dX(:,:)

      proj = .true.
      if( present( projector ) ) proj = projector
      eps = 1e-5_dp
      if( present( eps_deg ) ) eps = eps_deg

      nmatkq = Gkqset%ngk(1, ik) + nlotot
      nst = lst - fst + 1

      allocate( dX(nmatkq, fst:lst) )

      deval = 0.0_dp
      do ist = fst, lst
        dX(:, ist) = dHmat(1:nmatkq, ist) - evalk(ist) * dSmat(1:nmatkq, ist)
        if( gamma ) deval(ist) = dble( dX(ist, ist) )
        if( occk(ist) < input%groundstate%epsocc ) cycle
        do jst = 1, nmatkq
          dev = evalk(ist) - evalkq(jst)
          sig = sign( 1.0_dp, dev )
          dev = abs( dev )
          if( proj ) then
            if( dev < eps ) then
              t1 = 0.0_dp; t2 = -0.5_dp
            else
              t1 = sig * max( 0.0_dp, occmax - occkq(jst) ) / (dev + eps)
              t2 = -0.5_dp * occkq(jst)
            end if
            t1 = t1 / occmax; t2 = t2 / occmax
          else
            if( gamma .and. ist == jst ) then
              t1 = 0.0_dp; t2 = 0.0_dp
            else if( dev < eps ) then
              t1 = 0.0_dp; t2 = -0.5_dp
            else
              t1 = sig / dev; t2 = 0.0_dp
            end if
          end if
          dX(jst, ist) = t1 * dX(jst, ist) + t2 * dSmat(jst, ist)
        end do
      end do

      call zgemm( 'n', 'n', nmatkq, nst, nmatkq, zone, &
             eveckq, size( eveckq, dim=1 ), &
             dX, nmatkq, zzero, &
             devec, size( devec, dim=1 ) )

      deallocate( dX )
    end subroutine ph_eig_sternheimer

    !> This subroutine computes the contribution to the overlap matrix response
    !> and the constant part of the Hamiltonian matrix response coming from the 
    !> muffin-tin integrals.
    !>
    !> I.e., it computes
    !> \[ \delta^{\bf q}_{I \mu} S_{mn}^{\rm MT}({\bf k})
    !>    = \sum\limits_{\kappa, i} p^{I \mu}_{\kappa\alpha}({\bf q})\, \int\limits_{{\rm MT} \kappa} {\rm d}{\bf r} \left[ 
    !>      \left( \breve{\bar{\delta}}\phantom{}^{\bf -q}_{\kappa\alpha} \psi_{m{\bf k+q}}({\bf r}) \right)^\ast \,
    !>      \psi_{n{\bf k}}({\bf r}) + 
    !>      \psi_{m{\bf k+q}}^\ast({\bf r}) \,
    !>      \breve{\bar{\delta}}\phantom{}^{\bf q}_{\kappa\alpha} \psi_{n{\bf k}}({\bf r}) \right] \]
    !> and
    !> \[ \delta^{\bf q}_{I \mu} H_{mn}^{\rm MT}({\bf k})
    !>    = \sum\limits_{\kappa,\alpha} p^{I \mu}_{\kappa\alpha}({\bf q})\, \int\limits_{{\rm MT} \kappa} {\rm d}{\bf r} \left[ 
    !>      \left( \breve{\bar{\delta}}\phantom{}^{\bf -q}_{\kappa\alpha} \psi_{m{\bf k+q}}({\bf r}) \right)^\ast \,
    !>      \hat{\bf H}\, \psi_{n{\bf k}}({\bf r}) + 
    !>      \psi_{m{\bf k+q}}^\ast({\bf r})\, \hat{\bf H}\,
    !>      \breve{\bar{\delta}}\phantom{}^{\bf q}_{\kappa\alpha} \psi_{n{\bf k}}({\bf r}) \right] \;, \]
    !> where 
    !> \[ \breve{\bar{\delta}}\phantom{}^{\bf q}_{\kappa\alpha} \psi_{n{\bf k}}({\bf r}_\kappa)
    !>    = \sum\limits_{\bf G+k} C_{{\bf G+k}\,n}({\bf k}) \sum_{l,m} \sum_\xi \delta^{\bf q}_{\kappa\alpha} A^{\kappa}_{{\bf G+k},lm,\xi}\, 
    !>      u^\kappa_{l,\xi}(r_\kappa)\, Y_{lm}(\hat{\bf r}_\kappa) \;.\]
    !> See [[gen_dapwalm(subroutine)]] for the calculation of the response of the matching coefficients
    !> \(A^\kappa_{{\bf G+k},lm,\xi}\).
    !>
    !> @note The result is added to the input matrix! @endnote
    subroutine gen_dSH0_mt( ik, Gkset, Gkqset, fst1, lst1, fst2, lst2, eveck, eveckq, apwalmk, apwalmkq, pat, dSmat, dHmat, &
        diagonal )
      use constants, only: zone
      use matrix_elements
      use mod_kpointset, only: Gk_set
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k}\) and \({\bf G+k+q}\) vectors
      type(Gk_set), intent(in) :: Gkset, Gkqset
      !> first and last state on the left for which the matrix elements are calculated
      integer, intent(in) :: fst1, lst1
      !> first and last state on the right for which the matrix elements are calculated
      integer, intent(in) :: fst2, lst2
      !> eigenvectors at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: eveck(:,:), eveckq(:,:)
      !> (L)APW matching coefficients \(A^\kappa_{{\bf G+p},lm,\xi}\) at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)
      !> displacement pattern \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> overlap and Hamiltonian response
      complex(dp), intent(inout) :: dSmat(:,:), dHmat(:,:)
      !> conmpute only `fst2` to `lst2` diagonal elements 
      !> and store them in the first dimension (default: `.false.`)
      logical, optional, intent(in) :: diagonal

      integer :: ngk, ngkq, nst1, nst2
      integer :: is, ia, ias
      logical :: unpert, diag

      complex(dp), allocatable :: dapwalmk(:,:,:), dapwalmkq(:,:,:)

      diag = .false.
      if( present( diagonal ) ) diag = diagonal
      unpert = (sum( abs( pat ) ) < 1e-12_dp)

      ngk = Gkset%ngk(1, ik)
      ngkq = Gkqset%ngk(1, ik)
      nst1 = lst1 - fst1 + 1
      nst2 = lst2 - fst2 + 1

      allocate( dapwalmk, source=apwalmk(:, :, :, 1) )
      allocate( dapwalmkq, source=apwalmkq(:, :, :, 1) )

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          if( .not. unpert .and. sum( abs( pat(:, ias ) ) ) < 1e-12_dp ) cycle
          ! get perturbed matching coefficients at k
          call gen_dapwalm( ngk, Gkset%vgkc(:, :, 1, ik), pat(:, ias), apwalmk(:, :, :, ias), dapwalmk )
          ! get perturbed matching coefficients at k+q
          ! Note: conjugated pattern for pertubation with -q
          call gen_dapwalm( ngkq, Gkqset%vgkc(:, :, 1, ik), conjg( pat(:, ias) ), apwalmkq(:, :, :, ias), dapwalmkq )
          ! overlap response
          call me_mt_mat( is, ias, ngkq, ngk, apwalmkq(:, :, :, ias), dapwalmk, &
                 eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
                 zone, Smat_mt_basis(:, :, ias), zone, dSmat, &
                 diagonal_only=diag, left_local_orbitals=.true., right_local_orbitals=unpert )
          call me_mt_mat( is, ias, ngkq, ngk, dapwalmkq, apwalmk(:, :, :, ias), &
                 eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
                 zone, Smat_mt_basis(:, :, ias), zone, dSmat, &
                 diagonal_only=diag, left_local_orbitals=unpert, right_local_orbitals=.true. )
          ! Hamiltonian
          call me_mt_mat( is, ias, ngkq, ngk, apwalmkq(:, :, :, ias), dapwalmk, &
                 eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
                 zone, Hmat_mt_basis(:, :, ias), zone, dHmat, &
                 diagonal_only=diag, left_local_orbitals=.true., right_local_orbitals=unpert )
          call me_mt_mat( is, ias, ngkq, ngk, dapwalmkq, apwalmk(:, :, :, ias), &
                 eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
                 zone, Hmat_mt_basis(:, :, ias), zone, dHmat, &
                 diagonal_only=diag, left_local_orbitals=unpert, right_local_orbitals=.true. )
        end do
      end do

      deallocate( dapwalmk, dapwalmkq )
    end subroutine gen_dSH0_mt

    !> This subroutine computes the contribution to the overlap matrix response
    !> and the constant part of the Hamiltonian matrix response coming from the 
    !> variation of the characteristic function.
    !>
    !> I.e., it computes
    !> \[ \delta^{\bf q}_{I \mu} S_{mn}^{\rm IR}({\bf k})
    !>    = \sum\limits_{\kappa,\alpha} p^{I \mu}_{\kappa\alpha}({\bf q}) 
    !>      \int\limits_{\Omega} \left[ \psi_{m{\bf k+q}}^\ast({\bf r})\,
    !>      \psi_{n{\bf k}}({\bf r}) \right]_{\rm IR}\, \delta^{\bf q}_{\kappa\alpha} \Theta({\bf r})\, {\rm d}^3r \]
    !> and
    !> \[ \delta^{\bf q}_{I \mu} H_{mn}^{\rm IR}({\bf k})
    !>    = \sum\limits_{\kappa,\alpha} p^{I \mu}_{\kappa\alpha}({\bf q}) 
    !>      \int\limits_{\Omega} \left[ \psi_{m{\bf k+q}}^\ast({\bf r})\,
    !>      \hat{\bf H}\, \psi_{n{\bf k}}({\bf r}) \right]_{\rm IR}\, \delta^{\bf q}_{\kappa\alpha} \Theta({\bf r})\, {\rm d}^3r  \;, \]
    !> where \([\dots]_{\rm IR}\) means that the interstitial represenation of the integrand
    !> has to be used.
    !>
    !> See [[gen_dcfun_ig(subroutine)]] for the Foruier transform of the smooth characteristic 
    !> function response \(\delta^{\bf q}_{\kappa\alpha} \Theta({\bf r})\).
    !>
    !> @note The result is added to the input matrix! @endnote
    subroutine gen_dSH0_ir( ik, Gkset, Gkqset, fst1, lst1, fst2, lst2, eveck, eveckq, pat, dSmat, dHmat, &
        diagonal )
      use matrix_elements
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_opig
      use constants, only: zzero, zone
      use physical_constants, only: alpha
      use mod_kpointset, only: Gk_set
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_potential_and_density, only: pot_ir => veffir
      use modinput
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k}\) and \({\bf G+k+q}\) vectors
      type(Gk_set), intent(in) :: Gkset, Gkqset
      !> first and last state on the left for which the matrix elements are calculated
      integer, intent(in) :: fst1, lst1
      !> first and last state on the right for which the matrix elements are calculated
      integer, intent(in) :: fst2, lst2
      !> eigenvectors at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: eveck(:,:), eveckq(:,:)
      !> displacement pattern \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> overlap and Hamiltonian response
      complex(dp), intent(inout) :: dSmat(:,:), dHmat(:,:)
      !> conmpute only `fst2` to `lst2` diagonal elements 
      !> and store them in the first dimension (default: `.false.`)
      logical, optional, intent(in) :: diagonal

      integer :: ngk, ngkq, nst1, nst2, i
      integer :: is, ia, ias
      logical :: unpert, diag

      real(dp), allocatable :: kin_ir(:)
      complex(dp), allocatable :: dcfun_ig(:), dcfun_ir(:), pot_dcfun_ig(:), kin_dcfun_ig(:), zfft(:)

      diag = .false.
      if( present( diagonal ) ) diag = diagonal
      unpert = (sum( abs( pat ) ) < 1e-12_dp)

      ngk = Gkset%ngk(1, ik)
      ngkq = Gkqset%ngk(1, ik)
      nst1 = lst1 - fst1 + 1
      nst2 = lst2 - fst2 + 1

      allocate( dcfun_ig(ph_Gqset%ngvec), dcfun_ir(dfpt_2Gset%ngrtot), kin_ir(dfpt_Gset%ngrtot) )
      allocate( pot_dcfun_ig(ph_Gqset%ngvec), kin_dcfun_ig(ph_Gqset%ngvec) )
      allocate( zfft(dfpt_2Gset%ngrtot) )

      ! compute interstitial kinetic energy
      if( input%groundstate%ValenceRelativity == 'none' ) then
        kin_ir = 0.5_dp
      else
        kin_ir = 0.5_dp / (1.0_dp - 0.5_dp * alpha**2 * pot_ir )
      end if
      ! get characteristic function response on G+q points
      dcfun_ir = zzero
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          if( .not. unpert .and. sum( abs( pat(:, ias) ) ) < 1e-12_dp ) cycle
          call gen_dcfun_ig( ph_2Gqset%ngvec, ph_2Gqset%gc, ph_2Gqset%vgc, is, ia, pat(:, ias), dcfun_ir )
        end do
      end do
      call ph_2Gqset%change_set( ph_Gqset, dcfun_ir, dcfun_ig, 'pull', ng=ph_Gqset%ngvec )
      ! get characteristic function response on dense FFT grid
      call ph_2Gqset%change_set( dfpt_2Gset, dcfun_ir, zfft, 'push' )
      call dfpt_2Gset%ig2igfft( zfft, dcfun_ir )
      ! transform characteristic function response to dense real-space grid
      call zfftifc( 3, dfpt_2Gset%ngrid, 1, dcfun_ir )
      ! multiply with potential and transform back to reciprocal space
      call me_lapwlo_ir_opig( zone, pot_ir, dfpt_Gset, dcfun_ir, zzero, pot_dcfun_ig, ph_Gqset )
      ! multiply with kinetic energy and transform back to reciprocal space
      call me_lapwlo_ir_opig( zone, kin_ir, dfpt_Gset, dcfun_ir, zzero, kin_dcfun_ig, ph_Gqset )

      ! ** overlap
      call me_ir_mat( Gkqset, ik, Gkset, ik, &
             eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
             zone, dcfun_ig, zone, dSmat, &
             Gset_op=ph_Gqset, diagonal_only=diag )

      ! ** Hamiltonian
      ! potential
      call me_ir_mat( Gkqset, ik, Gkset, ik, &
             eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
             zone, pot_dcfun_ig, zone, dHmat, &
             Gset_op=ph_Gqset, diagonal_only=diag )
      ! kinetic energy
      do i = 1, 3
        call me_ir_mat( Gkqset, ik, Gkset, ik, &
               eveckq(:, fst1:lst1), eveck(:, fst2:lst2), &
               zone, kin_dcfun_ig, zone, dHmat, &
               left_gradient=i, right_gradient=i, &
               Gset_op=ph_Gqset, diagonal_only=diag )
      end do

      deallocate( dcfun_ig, dcfun_ir, kin_ir, pot_dcfun_ig, kin_dcfun_ig, zfft )
    end subroutine gen_dSH0_ir

    !> This subroutine computes the perturbed (L)APW matching coefficients
    !> for a single atom.
    !>
    !> They are given by
    !> \[ \delta^{\bf q}_{I \mu} A^\kappa_{{\bf G+k},lm,\xi} 
    !>    = \sum_\alpha p^{I \mu}_{\kappa\alpha}({\bf q})\, {\rm i} (G_\alpha + k_\alpha)\, A^\kappa_{{\bf G+k},lm,\xi} \;.\]
    subroutine gen_dapwalm( ngp, vgpc, pat, apwalm, dapwalm )
      !> number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> \({\bf G+p}\) vectors in Cartesian coordinates
      real(dp), intent(in) :: vgpc(3,*)
      !> dispalcement pattern \(p^{I \mu}_{\kappa\alpha}({\bf q})\) of the atom \(\kappa\)
      complex(dp), intent(in) :: pat(3)
      !> unperturbed matching coefficients \(A^\kappa_{{\bf G+k},lm,\xi}\)
      complex(dp), intent(in) :: apwalm(:,:,:)
      !> perturbed matching coefficients \(\delta^{\bf q}_{I \mu} A^\kappa_{{\bf G+k},lm,\xi}\)
      complex(dp), intent(out) :: dapwalm(:,:,:)

      integer :: i, j

      real(dp), allocatable :: t(:)
      complex(dp), allocatable :: z(:)

      allocate( t(2*ngp), z(ngp) )
      call dgemv( 't', 3, ngp, -1._dp, vgpc, 3, aimag(pat), 1, 0._dp, t(1), 1 )
      call dgemv( 't', 3, ngp,  1._dp, vgpc, 3, dble(pat),  1, 0._dp, t(ngp+1), 1 )
      z = cmplx( t(:ngp), t(ngp+1:), dp )

!$omp parallel default(shared) private(i,j)
!$omp do collapse(2)
      do j = 1, size( dapwalm, dim=3 )
        do i = 1, size( dapwalm, dim=2 )
          dapwalm(1:ngp, i, j) = z * apwalm(1:ngp, i, j)
        end do
      end do
!$omp end do
!$omp end parallel

      deallocate( t, z )
    end subroutine gen_dapwalm

    !> Compute the Fourier coefficients of the smooth characteristic function response
    !> for a phonon-like perturbation of a given atom.
    !>
    !> The smooth characteristic function response is given by
    !> \[ \delta^{\bf q}_{\kappa\alpha} \Theta({\bf r}) = 
    !>    \sum\limits_{\bf G} \delta^{\bf q}_{\kappa\alpha}\hat{\Theta}({\bf G}+{\bf q})\, 
    !>    {\rm e}^{{\rm i} ({\bf G}+{\bf q}) \cdot {\bf r}} \;, \]
    !> with
    !> \[ \delta^{\bf q}_{\kappa\alpha}\hat{\Theta}({\bf G}+{\bf q}) = 
    !>    \frac{4\pi\, {\rm i}\, R_\kappa^3}{\Omega} \frac{j_1(|{\bf G}+{\bf q}|R_\kappa)}{|{\bf G}+{\bf q}|R_\kappa}
    !>    {\rm e}^{-{\rm i} ({\bf G}+{\bf q}) \cdot {\bf r}}\, (G_\alpha + q_\alpha) \;. \]
    !> This routine will add
    !> \[ \delta^{\bf q}_{I \mu, \kappa} \hat{\Theta}({\bf G}+{\bf q}) = \sum_\alpha p^{I \mu}_{\kappa\alpha}({\bf q}) \,
    !>    \delta^{\bf q}_{\kappa\alpha}\hat{\Theta}({\bf G}+{\bf q}) \]
    !> to the input array `dcfun_ig`.
    subroutine gen_dcfun_ig( ngq, gqc, vgqc, is, ia, pat, dcfun_ig )
      use constants, only: fourpi
      use mod_lattice, only: omega
      use mod_atoms, only: atposc
      use mod_muffin_tin, only: rmt
      use modinput
      !> number of \({\bf G}+{\bf q}\) vectors
      integer, intent(in) :: ngq
      !> length of \({\bf G}+{\bf q}\) vectors
      real(dp), intent(in) :: gqc(*)
      !> \({\bf G}+{\bf q}\) vectors in Cartesian coordinates
      real(dp), intent(in) :: vgqc(3, *)
      !> species index
      integer, intent(in) :: is
      !> atom index
      integer, intent(in) :: ia
      !> displacement pattern \({\bf p}^{I \mu}_\kappa ({\bf q})\) of the atom
      complex(dp), intent(in) :: pat(3)
      !> Fourier coefficients \(\delta^{\bf q}_{\kappa\alpha}\hat{\Theta}({\bf G}+{\bf q})\)
      complex(dp), intent(inout) :: dcfun_ig(*)

      integer :: i, igq
      real(dp) :: t1, t2
      complex(dp) :: z1

      integer, allocatable :: finite_g_indices(:)

      t1 = fourpi / omega

      finite_g_indices = pack( [(igq,igq=1,ngq)], [(gqc(igq) > input%structure%epslat, igq=1, ngq)] )
      
!$omp parallel default( shared ) private( i, igq, t2, z1 )
!$omp do
      do i = 1, size( finite_g_indices )
        igq = finite_g_indices(i)
        t2 = gqc(igq) * rmt(is)
        z1 = t1 * dot_product( vgqc(:, igq), pat ) * (sin(t2) - t2 * cos(t2)) / (gqc(igq)**3)
        t2 = dot_product( vgqc(:, igq), atposc(:, ia, is) )
        dcfun_ig(igq) = dcfun_ig(igq) + z1 * cmplx( sin(t2), cos(t2), dp )
      end do
!$omp end do
!$omp end parallel
    end subroutine gen_dcfun_ig

end module phonons_eigensystem
