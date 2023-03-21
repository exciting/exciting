!> This module contains procedures for the calculation of the response
!> of the force acting on the nuclei in the crystal upon a phonon-like
!> perturbation \(\delta^{\bf q}_{I \mu}\).
!>
!> The response of the force \({\bf F}_\beta\) acting on nuclei \(\beta\)
!> upon a canonical phonon-like perturbation \(\delta^{\bf q}_{\alpha i}\)
!> reads
!> \[ \delta^{\bf q}_{\alpha i} F_{\beta j}
!>    = \sum_{\bf R} {\rm e}^{{\rm i}{\bf q}\cdot{\bf R}}\,
!>      \frac{\partial F_{\beta j}}{\partial \tau_{\alpha {\bf R} i}}
!>    = \sum_{\bf R} {\rm e}^{{\rm i}{\bf q}\cdot{\bf R}}\,
!>      \frac{\partial^2 E}{\partial \tau_{\alpha {\bf R} i}\, \partial \tau_{\beta {\bf 0} j}}
!>    = \sum_{\bf R} {\rm e}^{{\rm i}{\bf q}\cdot{\bf R}}\,
!>      \Phi_{\alpha i, \beta j}({\bf R},{\bf 0})
!>    = \sqrt{M_\alpha M_\beta}\, D_{\alpha i, \beta j}({\bf q}) \;, \]
!> which is exactly the Fourier transform of the interatomic force constants
!> \({\bf \Phi}({\bf R},{\bf R}')\) and hence directly give
!> the entries of the dynamical matrix \({\bf D}({\bf q})\).
!>
!> The total force can be separated into three terms: The Hellmann-Feynman
!> force \({\bf F}^{\rm HF}_\beta\), the Pulay force \({\bf F}^{\rm Pulay}_\beta\)
!> and the surface force \({\bf F}^{\rm SF}_\beta\).
!> The Pulay force has two origins. First, that the wavefunctions are only
!> variational eigenfunctions of the Hamiltonian and second, that we use
!> an atom position dependent basis and partitioning of the unit cell.
!> The surface force accounts for the discontinuity of the density and the potential
!> on the muffin-tin sphere surfaces.
module phonons_force
  use dfpt_variables
  use phonons_variables

  use precision, only: dp

  implicit none
  private

  public :: ph_frc_dpulay_k, ph_frc_dsurf_k
  public :: ph_frc_dpulay_int, ph_frc_dhf!, ph_frc_df0
  public :: ph_frc_symmetrize

  contains

    !> This subroutine calculates the contribution to the sum over states part
    !> of the Pulay force response coming from the given \({\bf k}\) point.
    !>
    !> The Pulay force is given by
    !> \[ F^{\rm Pulay}_{\beta j}
    !>    = \sum\limits_{n,{\bf k}}^{\rm valence} w_{\bf k}\, f_{n{\bf k}} \left[
    !>      \langle \psi_{n{\bf k}} | \hat{\bf H} - \epsilon_{n{\bf k}} | 
    !>      \breve{\bar{\psi}}\phantom{}^{\beta j}_{n{\bf k}} \rangle_{{\rm MT}\beta}
    !>    + \langle \breve{\bar{\psi}}\phantom{}^{\beta j}_{n{\bf k}} |
    !>      \hat{\bf H} - \epsilon_{n{\bf k}} | \psi_{n{\bf k}} \rangle_{{\rm MT}\beta} \right]
    !>    - \int\limits_{{\rm MT}\beta} \nabla_j n({\bf r})\, V_{\rm eff}({\bf r})\, {\rm d}{\bf r} \;, \]
    !> where \(\breve{\bar{\psi}}\phantom{}^{\beta j}_{n{\bf k}}\) is the contribution
    !> to \(\frac{{\rm d} \psi_{n{\bf k}}}{{\rm d} \tau_{\beta j}}\) that comes from 
    !> the change in the matching coefficients, i.e., the soft part of the basis function
    !> change.
    !>
    !> This subroutine computes the response of the contribution from a single \({\bf k}\) 
    !> point to the first term upon a phonon-like perturbation \(\delta^{\bf q}_{I \mu}\).
    !>
    !> @note The result is added to the input array! @endnote
    subroutine ph_frc_dpulay_k( ik, kset, Gkset, Gkqset, fst, lst, evalk, devalk, occk, docck, eveck, deveck, apwalmk, apwalmkq, pat, gamma, dforce, &
        dHmat_mt_basis, order )
      use dfpt_eigensystem, only: Smat_mt_basis, Hmat_mt_basis
      use phonons_eigensystem, only: gen_dapwalm
      use matrix_elements
      use constants, only: zzero, zone
      use mod_kpointset, only: k_set, Gk_set
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf k}\) vectors
      type(k_set), intent(in) :: kset
      !> set of \({\bf G+k}\) and \({\bf G+k+q}\) vectors
      type(Gk_set), intent(in) :: Gkset, Gkqset
      !> first and last state to consider in the sum over states
      integer, intent(in) :: fst, lst
      !> eigenvalues at \({\bf k}\)
      real(dp), intent(in) :: evalk(:)
      !> eigenvalue response at \({\bf k}\)
      real(dp), intent(in) :: devalk(:)
      !> occupation numbers at \({\bf k}\)
      real(dp), intent(in) :: occk(:)
      !> occupation number response at \({\bf k}\)
      real(dp), intent(in) :: docck(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> eigenvector response at \({\bf k}\)
      complex(dp), intent(in) :: deveck(:,:)
      !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)
      !> displacement pattern \(p^{I \mu}_{\alpha i}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> true for Gamma point phonons
      logical, intent(in) :: gamma
      !> force response
      complex(dp), intent(inout) :: dforce(3, natmtot)
      !> radial integrals of effective potential response times Gaunt coefficients
      complex(dp), optional, intent(in) :: dHmat_mt_basis(:,:,:)
      !> perturbation order of basis functions, either `0` or `1` (default: `0`)
      integer, optional, intent(in) :: order

      integer :: ord, ngk, ngkq, nst
      integer :: ip, is, ia, ias, ist
      real(dp) :: t1
      complex(dp) :: pat0(3,3)

      complex(dp), allocatable :: dSmat(:,:), dHmat(:,:)
      complex(dp), allocatable :: dapwalmk(:,:,:), ddapwalmk(:,:,:), dapwalmkq(:,:,:)

      ord = 0
      if( present( order ) ) ord = order

      ngk = Gkset%ngk(1, ik)
      ngkq = Gkqset%ngk(1, ik)
      nst = lst - fst + 1

      allocate( dapwalmk, source=apwalmk(:, :, :, 1) )
      allocate( ddapwalmk, source=apwalmk(:, :, :, 1) )
      allocate( dapwalmkq, source=apwalmkq(:, :, :, 1) )
      allocate( dSmat(fst:lst, 1), dHmat(fst:lst, 1) )
      pat0 = zzero
      do ip = 1, 3
        pat0(ip, ip) = zone
      end do

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ip = 1, 3

            dSmat = zzero; dHmat = zzero

            ! ** contribution from wave function response
            ! get canonically perturbed matching coefficients at k
            call gen_dapwalm( ngk, Gkset%vgkc(:, :, 1, ik), pat0(:, ip), apwalmk(:, :, :, ias), dapwalmk )
            ! get canonically perturbed matching coefficients at k+q
            call gen_dapwalm( ngkq, Gkqset%vgkc(:, :, 1, ik), pat0(:, ip), apwalmkq(:, :, :, ias), dapwalmkq )
            if( ord == 1 ) then
              ! get irrep perturbed matching coefficients at k
              call gen_dapwalm( ngk, Gkset%vgkc(:, :, 1, ik), pat(:, ias), apwalmk(:, :, :, ias), ddapwalmk )
            else
              ddapwalmk = apwalmk(:, :, :, ias)
            end if
            ! overlap
            call me_mt_mat( is, ias, ngk, ngkq, ddapwalmk, dapwalmkq, &
                   eveck(:, fst:lst), deveck(:, fst:lst), &
                   zone, Smat_mt_basis(:, :, ias), zone, dSmat, &
                   diagonal_only=.true., left_local_orbitals=(ord==0), right_local_orbitals=.false. )
            ! Hamiltonian
            call me_mt_mat( is, ias, ngk, ngkq, ddapwalmk, dapwalmkq, &
                   eveck(:, fst:lst), deveck(:, fst:lst), &
                   zone, Hmat_mt_basis(:, :, ias), zone, dHmat, &
                   diagonal_only=.true., left_local_orbitals=(ord==0), right_local_orbitals=.false. )

            if( ord == 1 ) then
              ! get doubly perturbed matching coefficients at k
              call gen_dapwalm( ngk, Gkset%vgkc(:, :, 1, ik), pat(:, ias), dapwalmk, ddapwalmk )
            else
              ddapwalmk = dapwalmk
            end if
            ! overlap
            call me_mt_mat( is, ias, ngk, ngkq, ddapwalmk, apwalmkq(:, :, :, ias), &
                   eveck(:, fst:lst), deveck(:, fst:lst), &
                   zone, Smat_mt_basis(:, :, ias), zone, dSmat, &
                   diagonal_only=.true., left_local_orbitals=.false., right_local_orbitals=.true. )
            ! Hamiltonian
            call me_mt_mat( is, ias, ngk, ngkq, ddapwalmk, apwalmkq(:, :, :, ias), &
                   eveck(:, fst:lst), deveck(:, fst:lst), &
                   zone, Hmat_mt_basis(:, :, ias), zone, dHmat, &
                   diagonal_only=.true., left_local_orbitals=.false., right_local_orbitals=.true. )
            dSmat = 2 * dSmat; dHmat = 2 * dHmat

            ! ** contribution from potential response
            if( present( dHmat_mt_basis ) ) then
              call me_mt_mat( is, ias, ngk, ngk, dapwalmk, apwalmk(:, :, :, ias), &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, dHmat_mt_basis(:, :, ias), zone, dHmat, &
                     diagonal_only=.true., left_local_orbitals=.false., right_local_orbitals=.true. )
              call me_mt_mat( is, ias, ngk, ngk, apwalmk(:, :, :, ias), dapwalmk, &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, dHmat_mt_basis(:, :, ias), zone, dHmat, &
                     diagonal_only=.true., left_local_orbitals=.true., right_local_orbitals=.false. )
            end if

            ! sum over states
            do ist = fst, lst
              t1 = kset%wkpt(ik) * occk(ist)
              dforce(ip, ias) = dforce(ip, ias) + t1 * (dHmat(ist,1) - evalk(ist) * dSmat(ist,1))
            end do

            ! ** contribution from eigenvalue and occupation response
            if( gamma .and. ord == 0 ) then
              dSmat = zzero; dHmat = zzero
              ! overlap
              call me_mt_mat( is, ias, ngk, ngk, dapwalmk, apwalmk(:, :, :, ias), &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, Smat_mt_basis(:, :, ias), zone, dSmat, &
                     diagonal_only=.true., left_local_orbitals=.false., right_local_orbitals=.true. )
              call me_mt_mat( is, ias, ngk, ngk, apwalmk(:, :, :, ias), dapwalmk, &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, Smat_mt_basis(:, :, ias), zone, dSmat, &
                     diagonal_only=.true., left_local_orbitals=.true., right_local_orbitals=.false. )
              ! Hamiltonian
              call me_mt_mat( is, ias, ngk, ngk, dapwalmk, apwalmk(:, :, :, ias), &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, Hmat_mt_basis(:, :, ias), zone, dHmat, &
                     diagonal_only=.true., left_local_orbitals=.false., right_local_orbitals=.true. )
              call me_mt_mat( is, ias, ngk, ngk, apwalmk(:, :, :, ias), dapwalmk, &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, Hmat_mt_basis(:, :, ias), zone, dHmat, &
                     diagonal_only=.true., left_local_orbitals=.true., right_local_orbitals=.false. )
              ! sum over states
              do ist = fst, lst
                t1 = kset%wkpt(ik) * occk(ist)
                dforce(ip, ias) = dforce(ip, ias) - t1 * devalk(ist) * dSmat(ist, 1)
                t1 = kset%wkpt(ik) * docck(ist)
                dforce(ip, ias) = dforce(ip, ias) + t1 * (dHmat(ist, 1) - evalk(ist) * dSmat(ist, 1))
              end do
            end if

          end do
        end do
      end do

      deallocate( dapwalmk, ddapwalmk, dapwalmkq, dSmat, dHmat )
    end subroutine ph_frc_dpulay_k

    !> This subroutine calculates the contribution to the sum over states part
    !> of the Surface force response coming from the given \({\bf k}\) point.
    !>
    !> The Surface force is given by
    !> \[ F^{\rm SF}_{\beta j}
    !>    = \sum\limits_{n,{\bf k}}^{\rm valence} w_{\bf k}\, f_{n{\bf k}}\,
    !>      \oint\limits_{\partial{\rm MT}\beta} \left[
    !>      \psi_{n{\bf k}}^\ast({\bf r}) \left[ \hat{\bf T} - \epsilon_{n{\bf k}} \right]
    !>      \psi_{n{\bf k}}({\bf r}) \right]_{\rm IR} \hat{e}_j\, {\rm d}S
    !>    - \oint\limits_{\partial{\rm MT}\beta} \left[ n({\bf r})
    !>      \left( V_{\rm C}({\bf r}) + \epsilon_{\rm xc}({\bf r}) \right) \right]_{\rm SF}
    !>      \hat{e}_j\, {\rm d}S \;. \]
    !>
    !> This subroutine computes the response of the contribution from a single \({\bf k}\) 
    !> point to the first term upon a phonon-like perturbation \(\delta^{\bf q}_{I \mu}\).
    !>
    !> @note The result is added to the input array! @endnote
    subroutine ph_frc_dsurf_k( ik, kset, Gkset, Gkqset, fst, lst, evalk, devalk, occk, docck, eveck, deveck, pat, gamma, dforce, &
        dpot_ir, order )
      use phonons_eigensystem, only: gen_dcfun_ig
      use matrix_elements
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_opig
      use constants, only: zzero, zone, zi
      use physical_constants, only: alpha
      use mod_kpointset, only: k_set, Gk_set
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_potential_and_density, only: pot_ir => veffir
      use modinput
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf k}\) vectors
      type(k_set), intent(in) :: kset
      !> set of \({\bf G+k}\) and \({\bf G+k+q}\) vectors
      type(Gk_set), intent(in) :: Gkset, Gkqset
      !> first and last state to consider in the sum over states
      integer, intent(in) :: fst, lst
      !> eigenvalues at \({\bf k}\)
      real(dp), intent(in) :: evalk(:)
      !> eigenvalue response at \({\bf k}\)
      real(dp), intent(in) :: devalk(:)
      !> occupation numbers at \({\bf k}\)
      real(dp), intent(in) :: occk(:)
      !> occupation number response at \({\bf k}\)
      real(dp), intent(in) :: docck(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> eigenvector response at \({\bf k}\)
      complex(dp), intent(in) :: deveck(:,:)
      !> displacement pattern \(p^{I \mu}_{\alpha i}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> true for Gamma point phonons
      logical, intent(in) :: gamma
      !> force response
      complex(dp), intent(inout) :: dforce(3, natmtot)
      !> interstitial effective potential response on real space FFT grid
      complex(dp), optional, intent(in) :: dpot_ir(:)
      !> perturbation order of basis functions, either `0` or `1` (default: `0`)
      integer, optional, intent(in) :: order

      integer :: ord, ngk, ngkq, nst
      integer :: ip, is, ia, ias, ist
      integer :: igq, i
      real(dp) :: t1
      complex(dp) :: pat0(3,3)

      complex(dp), allocatable :: dSmat(:,:), dHmat(:,:)
      real(dp), allocatable :: kin_ir(:)
      complex(dp), allocatable :: dkin_ir(:), dcfun0_ir(:), dcfun0_ig(:), kin_dcfun0_ig(:), dkin_dcfun0_ig(:), zfft(:)

      ord = 0
      if( present( order ) ) ord = order

      ngk = Gkset%ngk(1, ik)
      ngkq = Gkqset%ngk(1, ik)
      nst = lst - fst + 1

      allocate( kin_ir(dfpt_Gset%ngrtot) )
      allocate( dcfun0_ir(dfpt_2Gset%ngrtot), dcfun0_ig(ph_Gqset%ngvec) )
      allocate( kin_dcfun0_ig(ph_Gqset%ngvec) )
      allocate( zfft(dfpt_2Gset%ngrtot) )
      allocate( dSmat(fst:lst, 1), dHmat(fst:lst, 1) )

      pat0 = zzero
      do ip = 1, 3
        pat0(ip, ip) = zone
      end do

      ! compute interstitial kinetic energy
      if( input%groundstate%ValenceRelativity == 'none' ) then
        kin_ir = 0.5_dp
      else
        kin_ir = 0.5_dp / (1.0_dp - 0.5_dp * alpha**2 * pot_ir )
      end if
      ! compute interstitial kinetic energy response
      if( present( dpot_ir ) ) then
        allocate( dkin_ir(dfpt_Gset%ngrtot) )
        allocate( dkin_dcfun0_ig(ph_Gqset%ngvec) )
        if( input%groundstate%ValenceRelativity == 'none' ) then
          dkin_ir = zzero
        else
          dkin_ir = 0.5_dp * (0.5_dp * alpha**2 * dpot_ir) / (1.0_dp - 0.5_dp * alpha**2 * pot_ir )**2
        end if
      end if

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ip = 1, 3
            
            ! ** contribution from wave function response
            dSmat = zzero; dHmat = zzero
            ! get canonically perturbed characteristic function on G-vectors
            dcfun0_ir = zzero
            call gen_dcfun_ig( ph_2Gqset%ngvec, ph_2Gqset%gc, ph_2Gqset%vgc, is, ia, pat0(:, ip), dcfun0_ir )
            if( ord == 1 ) then
              ! get doubly perturbed characteristic function on G+q-vectors
              do igq = 1, ph_2Gqset%ngvec
                dcfun0_ir(igq) = -zi * dot_product( ph_2Gqset%vgc(:, igq), pat(:, ias) ) * dcfun0_ir(igq)
              end do
            end if
            call ph_2Gqset%change_set( ph_Gqset, dcfun0_ir, dcfun0_ig, 'pull', ng=ph_Gqset%ngvec )
            ! get canonically perturbed characteristic function on dense real-space grid
            call ph_2Gqset%change_set( dfpt_2Gset, dcfun0_ir, zfft, 'push' )
            call dfpt_2Gset%ig2igfft( zfft, dcfun0_ir )
            call zfftifc( 3, dfpt_2Gset%ngrid, 1, dcfun0_ir )
            ! multiply with kinetic energy and transform back to reciprocal space
            call me_lapwlo_ir_opig( zone, kin_ir, dfpt_Gset, dcfun0_ir, zzero, kin_dcfun0_ig, ph_Gqset )

            ! overlap
            call me_ir_mat( Gkqset, ik, Gkset, ik, &
                   deveck(:, fst:lst), eveck(:, fst:lst), &
                   zone, dcfun0_ig, zone, dSmat, &
                   Gset_op=ph_Gqset, diagonal_only=.true. )
            ! Hamiltonian
            do i = 1, 3
              call me_ir_mat( Gkqset, ik, Gkset, ik, &
                     deveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, kin_dcfun0_ig, zone, dHmat, &
                     left_gradient=i, right_gradient=i, &
                     Gset_op=ph_Gqset, diagonal_only=.true. )
            end do
            dSmat = 2 * dSmat; dHmat = 2 * dHmat

            ! ** contribution from potential response
            if( present( dpot_ir ) .and. ord == 0 ) then
              ! multiply with kinetic energy response and transform back to reciprocal space
              call me_lapwlo_ir_opig( zone, dkin_ir, ph_Gqset, dcfun0_ir, zzero, dkin_dcfun0_ig, ph_Gqset )
              ! Hamiltonian
              do i = 1, 3
                call me_ir_mat( Gkset, ik, Gkset, ik, &
                       eveck(:, fst:lst), eveck(:, fst:lst), &
                       zone, dkin_dcfun0_ig, zone, dHmat, &
                       left_gradient=i, right_gradient=i, &
                       Gset_op=ph_Gqset, diagonal_only=.true. )
              end do
            end if

            ! sum over states
            do ist = fst, lst
              t1 = kset%wkpt(ik) * occk(ist)
              if( ord == 1 ) t1 = t1 / 2
              dforce(ip, ias) = dforce(ip, ias) + t1 * conjg( dHmat(ist, 1) - evalk(ist) * dSmat(ist, 1) )
            end do

            ! ** contribution from eigenvalue and occupation response
            if( gamma .and. ord == 0 ) then
              dSmat = zzero; dHmat = zzero
              ! overlap
              call me_ir_mat( Gkset, ik, Gkset, ik, &
                     eveck(:, fst:lst), eveck(:, fst:lst), &
                     zone, dcfun0_ig, zone, dSmat, &
                     Gset_op=ph_Gqset, diagonal_only=.true. )
              ! Hamiltonian
              do i = 1, 3
                call me_ir_mat( Gkset, ik, Gkset, ik, &
                       eveck(:, fst:lst), eveck(:, fst:lst), &
                       zone, kin_dcfun0_ig, zone, dHmat, &
                       left_gradient=i, right_gradient=i, &
                       Gset_op=ph_Gqset, diagonal_only=.true. )
              end do
              ! sum over states
              do ist = fst, lst
                t1 = kset%wkpt(ik) * occk(ist)
                dforce(ip, ias) = dforce(ip, ias) - t1 * devalk(ist) * dSmat(ist, 1)
                t1 = kset%wkpt(ik) * docck(ist)
                dforce(ip, ias) = dforce(ip, ias) + t1 * (dHmat(ist, 1) - evalk(ist) * dSmat(ist,1))
              end do
            end if

          end do
        end do
      end do

      deallocate( kin_ir, dcfun0_ir, dcfun0_ig, kin_dcfun0_ig, zfft )
      if( present( dpot_ir ) ) &
        deallocate( dkin_ir, dkin_dcfun0_ig )

    end subroutine ph_frc_dsurf_k

    !> This subroutine calculates the contribution of the integral term
    !> to the Pulay force response.
    !>
    !> The Pulay force is given by
    !> \[ F^{\rm Pulay}_{\beta j}
    !>    = \sum\limits_{n,{\bf k}}^{\rm valence} w_{\bf k}\, f_{n{\bf k}} \left[
    !>      \langle \psi_{n{\bf k}} | \hat{\bf H} - \epsilon_{n{\bf k}} | 
    !>      \breve{\bar{\psi}}\phantom{}^{\beta j}_{n{\bf k}} \rangle_{{\rm MT}\beta}
    !>    + \langle \breve{\bar{\psi}}\phantom{}^{\beta j}_{n{\bf k}} |
    !>      \hat{\bf H} - \epsilon_{n{\bf k}} | \psi_{n{\bf k}} \rangle_{{\rm MT}\beta} \right]
    !>    - \int\limits_{{\rm MT}\beta} \nabla_j n({\bf r})\, V_{\rm eff}({\bf r})\, {\rm d}{\bf r} \;. \]
    !>
    !> This subroutine computes the response of the second integral term upon a 
    !> phonon-like perturbation \(\delta^{\bf q}_{I \mu}\).
    !>
    !> @note The result is added to the input array! @endnote
    subroutine ph_frc_dpulay_int( drho_mt, dpot_mt, dforce )
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      use mod_muffin_tin, only: lmmaxvr, nrmtmax, nrmt
      use mod_potential_and_density, only: rho_mt => rhomt, pot_mt => veffmt
      use modinput
      !> soft muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(in) :: drho_mt(:,:,:)
      !> soft muffin-tin effective potential response as complex spherical harmonics expansion
      complex(dp), intent(in) :: dpot_mt(:,:,:)
      !> force response
      complex(dp), intent(inout) :: dforce(3, natmtot)

      integer :: is, ia, ias, ip, nr, ir

      complex(dp), allocatable :: zfmt(:,:,:)
      
      complex(dp), external :: zfmtinp

      allocate( zfmt(lmmaxvr, nrmtmax, 0:3) )
      do is = 1, nspecies
        nr = nrmt(is)
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! transform effective potential into complex function
          do ir = 1, nr
            call rtozflm( input%groundstate%lmaxvr, pot_mt(:, ir, ias), zfmt(:, ir, 0) )
          end do
          ! get gradient of density response
          call gradzfmt( input%groundstate%lmaxvr, nr, spr(1:nr, is), lmmaxvr, nrmtmax, drho_mt(:, 1:nr, ias), zfmt(:, :, 1:3) )
          ! add integral to force response
          do ip = 1, 3
            dforce(ip, ias) = dforce(ip, ias) &
              - zfmtinp( .true., input%groundstate%lmaxvr, nr, spr(1:nr, is), lmmaxvr, zfmt(:, :, 0), zfmt(:, :, ip) )
          end do
          ! transform density into complex function
          do ir = 1, nr
            call rtozflm( input%groundstate%lmaxvr, rho_mt(:, ir, ias), zfmt(:, ir, 0) )
          end do
          ! get gradient of density
          call gradzfmt( input%groundstate%lmaxvr, nr, spr(1:nr, is), lmmaxvr, nrmtmax, zfmt(:, 1:nr, 0), zfmt(:, :, 1:3) )
          ! add integral to force response
          do ip = 1, 3
            dforce(ip, ias) = dforce(ip, ias) &
              - zfmtinp( .true., input%groundstate%lmaxvr, nr, spr(1:nr, is), lmmaxvr, zfmt(:, :, ip), dpot_mt(:, :, ias) )
          end do
        end do
      end do
      deallocate( zfmt )
    end subroutine ph_frc_dpulay_int

    !> This subroutine calculates the response of the Hellmann-Feynman force.
    !>
    !> The Hellmann-Feynman force acting of atom \(\beta\) is defined as the negative 
    !> gradient of the Coulomb potential induced by all charges except the nucleus
    !> charge of atom \(\beta\) evaluated at the position of the nucleus
    !> \[ F^{\rm HF}_{\beta j}
    !>    = -Z_{\beta} \sum_{\bf R} \lim\limits_{{\bf r} \rightarrow {\bf \tau}_{\beta {\bf R}}}
    !>      \nabla_j \left[ \sum_{\alpha \neq \beta} \frac{Z_{\alpha}}{|{\bf r}-{\bf \tau}_\alpha|} 
    !>      -\int \frac{n({\bf r}')}{|{\bf r}-{\bf r}'|} {\rm d}{\bf r}' \right] \;. \]
    !>
    !> This subroutine computes the response of the Hellmann-Feynman force upon a 
    !> phonon-like perturbation \(\delta^{\bf q}_{I \mu}\).
    !>
    !> @note The result is added to the input array! @endnote
    subroutine ph_frc_dhf( drho_mt, dpot_coul_mt, dforce )
      use constants, only: fourpi, sqrt_two
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr, spzn
      use mod_muffin_tin, only: nrmtmax, nrmt, rmt
      !> soft muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(in) :: drho_mt(:,:,:)
      !> soft muffin-tin Coulomb potential response as complex spherical harmonics expansion
      complex(dp), intent(in) :: dpot_coul_mt(:,:,:)
      !> force response
      complex(dp), intent(inout) :: dforce(3, natmtot)

      integer :: is, ia, ias, nr, ir, m, lm
      real(dp) :: t1, t2, t3

      real(dp), allocatable :: gr(:), cf(:,:)
      complex(dp), allocatable :: zfmt(:)

      t3 = sqrt( fourpi / 3.0_dp )
      allocate( gr(nrmtmax), cf(3, nrmtmax) )
      allocate( zfmt(nrmtmax) )
      do is = 1, nspecies
        nr = nrmt(is)
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do m = -1, 1
            lm = 3 + m
            do ir = 1, nr
              zfmt(ir) = drho_mt(lm, ir, ias) * (1.0_dp - (spr(ir, is) / rmt(is))**3)
            end do
            call fderiv( -1, nr, spr(1:nr, is), dble( zfmt ), gr, cf )
            t1 = spzn(is) * (t3 * gr(nr) + dble( dpot_coul_mt(lm, nr, ias) ) / (rmt(is) * t3))
            call fderiv( -1, nr, spr(1:nr, is), aimag( zfmt ), gr, cf)
            t2 = spzn(is) * (t3 * gr(nr) + aimag( dpot_coul_mt(lm, nr, ias) ) / (rmt(is) * t3))
            if( m == -1 ) then
              dforce(1, ias) = dforce(1, ias) + cmplx( t1,  t2, dp ) / sqrt_two
              dforce(2, ias) = dforce(2, ias) + cmplx( t2, -t1, dp ) / sqrt_two
            else if( m == 1) then
              dforce(1, ias) = dforce(1, ias) - cmplx( t1,  t2, dp ) / sqrt_two
              dforce(2, ias) = dforce(2, ias) + cmplx( t2, -t1, dp ) / sqrt_two
            else
              dforce(3, ias) = dforce(3, ias) + cmplx( t1, t2, dp )
            end if
          end do
        end do
      end do
      deallocate( gr, cf, zfmt )
    end subroutine ph_frc_dhf

    !> This subroutine symmetrizes the phonon-like perturbed atomic force
    !> in the basis of an irrep. See also [[phonons_symmetry(module)]].
    !>
    !> If `acoustic_sum_rule=.true.`, then the acoustic sum rule is imposed, i.e.,
    !> \[ \sum \delta^{\bf q}_{I \mu} {\bf F}_{\beta} = 0 \;. \]
    !> Note, that this only holds for \({\bf q} = {\bf 0}\).
    subroutine ph_frc_symmetrize( dforce, vql, dirrep, nsym, isym, ivsym, symmat, &
        acoustic_sum_rule )
      use constants, only: zzero, twopi
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, natmmax
      use mod_symmetry, only: lsplsymc, symlat, symlatc, ieqatom
      use modinput
      !> force response
      complex(dp), intent(inout) :: dforce(3, natmtot, *)
      !> Bloch wavevector \({\bf q}\) in lattice coordinates
      real(dp), intent(in) :: vql(3)
      !> dimension \(d_I\) of the irrep
      integer, intent(in) :: dirrep
      !> number of symmetry operations in small group of \({\bf q}\)
      integer, intent(in) :: nsym
      !> indices of symmetries in global arrays
      integer, intent(in) :: isym(:)
      !> lattice vectors that map \({\bf \rm S}{\bf q}\) back to 1st BZ
      integer, intent(in) :: ivsym(3,*)
      !> matrix representation of symmetries in the basis of the irrep \({\bf \texttt{S}}^I\)
      complex(dp), intent(in) :: symmat(:,:,:)
      !> impose acoustic sum rule (default: `.false.`)
      logical, optional, intent(in) :: acoustic_sum_rule

      integer :: is, ia, ias, ja, jas, d, dd, s, lspl, ip
      real(dp) :: sl(3,3), sc(3,3), a(3), b(3), c(3), phi
      complex(dp) :: z1
      logical :: sumrule

      complex(dp), allocatable :: dfc(:,:,:)

      sumrule = .false.
      if( present( acoustic_sum_rule ) ) sumrule = acoustic_sum_rule

      allocate( dfc(3, natmmax, dirrep) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! make a copy of the input force response and nullify it
          do d = 1, dirrep
            dfc(:, ia, d) = dforce(:, ias, d)
            dforce(:, ias, d) = zzero
          end do
        end do
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! loop over symmetries
          do s = 1, nsym
            lspl = lsplsymc(isym(s))                ! index of spatial rotation
            sl = dble( symlat(:, :, lspl) )         ! lattice spatial rotation
            sc = symlatc(:, :, lspl)                ! Cartesian spatial rotation
            ja = ieqatom(ia, is, isym(s))           ! equivalent atom
            jas = idxas(ja, is)
            a = input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord
            call r3mv( sl, a, b ) 
            b = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord - b! + vtlsymc(:,isym(s))
            phi = -twopi * dot_product( vql, b )
            do d = 1, dirrep
              ! apply rotation to force response
              a = dble(  dfc(:, ia, d) )
              b = aimag( dfc(:, ia, d) )
              call r3mtv( sc, a, c ); a = c
              call r3mtv( sc, b, c ); b = c
              ! add contribution to input array
              do dd = 1, dirrep
                dforce(:, jas, dd) = dforce(:, jas, dd) + symmat(d, dd, s) * cmplx( cos( phi ), sin( phi ), dp ) * cmplx( a, b, dp )
              end do
            end do
          end do
        end do
      end do
      
      ! normalize
      do d = 1, dirrep
        dforce(:, :, d) = dforce(:, :, d) / nsym
      end do

      ! impose acoustic sum rule
      if( sumrule ) then
        do d = 1, dirrep
          do ip = 1, 3
            z1 = sum( dforce(ip, :, d) ) / natmtot
            dforce(ip, :, d) = dforce(ip, :, d) - z1
          end do
        end do
      end if

     deallocate( dfc)
    end subroutine ph_frc_symmetrize

end module phonons_force
