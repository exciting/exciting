!> This module contains procedures for the calculation of the 
!> Hamiltonian response upon a perturbing electric field
!> and the solution of the corresponding Sternheimer equation.
module efield_eigensystem
  use dfpt_variables
  use dfpt_eigensystem

  use precision, only: dp

  implicit none
  private

  ! radial integrals times Gaunt coefficients of momentum operator
  complex(dp), allocatable :: Pmat_mt_basis(:,:,:,:)

  public :: ef_eig_init, ef_eig_free
  public :: ef_eig_gen_dHmat
  public :: ef_eig_sternheimer
  public :: ef_eig_rotate_devec

  contains

    !> This subroutine initializes variables for the calculation of matrix elements
    !> that remain constant during the entire calculation.
    !>
    !> This includes: 
    !> 
    !> * calculation of radial muffin-tin integrals times Gaunt coefficients 
    !>   for momentum matrix elements
    subroutine ef_eig_init
      ! compute MT radial integrals times Gaunt coefficients
      call gen_momentum_mt_basis( Pmat_mt_basis, lmax_apw=dfpt_lmaxapw )
    end subroutine ef_eig_init

    !> This subroutine frees memory from the module variables.
    subroutine ef_eig_free
      if( allocated( Pmat_mt_basis ) ) deallocate( Pmat_mt_basis )
    end subroutine ef_eig_free

    !> This subroutine computes the response of the Hamiltonian matrix 
    !> for the given \({\bf k}\) point.
    !>
    !> The full Hamiltonian response is given by
    !> \[ \delta_{\mathcal{E}_\alpha} H_{mn}({\bf k})
    !>    = \int\limits_{{\rm MT} \alpha} {\rm d}{\bf r}
    !>      \psi_{m{\bf k}}^\ast({\bf r})\, 
    !>      \left[ \delta_{\mathcal{E}_\alpha} V_{\rm eff}({\bf r}) - r_\alpha \right]
    !>      \psi_{n{\bf k}}({\bf r}) \;, \]
    !> where the matrix elements of the position operator form the so called constant
    !> part of the Hamiltonian response that does not change during the self-consistency
    !> cycle.
    !>
    !> If the polarization direction `ip` is given, then this subroutine computes the
    !> constant part of the Hamiltonian response (see [[gen_dH0(subroutine)]]).
    !>
    !> If the effective potential response is given in form of the radial muffin-tin
    !> integrals times Gaunt coefficients `dHmat_mt_basis` and the product with the 
    !> characteristic function in reciprocal space `dpot_cfun_ig`, then this subroutine
    !> computes the corresponding contribution to the Hamiltonian matrix response, 
    !> \(\delta^{\bf q}_{I \mu} {\bf H}^0\) (see [[dfpt_eig_gen_dHmat(subroutine)]]).
    !> 
    !> @note The result is added to the input matrix! @endnote
    subroutine ef_eig_gen_dHmat( ik, Gkset, fst, lst, evalk, eveck, apwalmk, dHmat, &
        ip, dHmat_mt_basis, dpot_cfun_ig, dkin_cfun_ig )
      use mod_kpointset, only: Gk_set
      use mod_APW_LO, only: nlotot
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k}\) vectors
      type(Gk_set), intent(in) :: Gkset
      !> first and last state for which the matrix elements are calculated
      integer, intent(in) :: fst, lst
      !> eigenvalues at \({\bf k}\)
      real(dp), intent(in) :: evalk(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) at \({\bf k}\)
      complex(dp), intent(in) :: apwalmk(:,:,:,:)
      !> Hamiltonian response
      complex(dp), intent(inout) :: dHmat(:,:)
      !> polarization direction
      integer, optional, intent(in) :: ip
      !> radial integrals of effective potential response times Gaunt coefficients
      complex(dp), optional, intent(in) :: dHmat_mt_basis(:,:,:)
      !> interstitial effective potential response times characteristic function in reciprocal space
      complex(dp), optional, intent(in) :: dpot_cfun_ig(:)
      !> interstitial (scalar relativistic) kinetic energy response times characteristic function in reciprocal space
      complex(dp), optional, intent(in) :: dkin_cfun_ig(:)

      integer :: nmatk

      nmatk = Gkset%ngk(1, ik) + nlotot

      if( present( ip ) ) then
        call gen_dH0( ik, Gkset, 1, nmatk, fst, lst, evalk, eveck, apwalmk, ip, dHmat )
      end if

      if( present( dHmat_mt_basis ) .and. present( dpot_cfun_ig ) .and. present( dkin_cfun_ig ) ) then
        ! contribution from potential response
        call dfpt_eig_gen_dHmat( ik, Gkset, Gkset, 1, nmatk, fst, lst, &
          eveck, eveck, apwalmk, apwalmk, dHmat_mt_basis, dpot_cfun_ig, dkin_cfun_ig, dHmat )
      end if
    end subroutine ef_eig_gen_dHmat

    !> This subroutine solves the Sternheimer equation for a perturbing electric field
    !> for the eigenvalue and eigenvector response at a given \({\bf k}\) point.
    !> 
    !> In general, Sternheimer's equation for the state \(\psi_{n{\bf k}}\) reads
    !> \[ \left[ \hat{\bf H} - \epsilon_{n{\bf k}} \right] \delta \psi_{n{\bf k}}({\bf r})
    !>    = -\left[ \delta \hat{\bf H} - \delta \epsilon_{n{\bf k}} \right] \psi_{n{\bf k}}({\bf r}) \;. \]
    !> In the case of a perturbing electric field the wavefunction response reads
    !> \[ \delta_{\mathcal{E}_\alpha} \psi_{n{\bf k}}({\bf r})
    !>    = \sum_\nu \delta_{\mathcal{E}_\alpha \mu} C_{\nu n}({\bf k}) \, \phi_{\nu{\bf k}}({\bf r})
    !>   \;, \]
    !> We want to solve Sternheimer's equation for the response of the expansion coefficients 
    !> (eigenvectors) and the Kohn-Sham energies (eigenvalues). 
    !> To this extend, we expand the eigenvector response in the eigenvectors
    !> of the unperturbed system at wavevector \({\bf k}\)
    !> \[ \delta_{\mathcal{E}_\alpha} C_{\nu n}({\bf k}) 
    !>    = \sum_m \delta_{\mathcal{E}_\alpha} X_{mn}({\bf k})\, C_{\nu m}({\bf k}) \;. \]
    !> This yields
    !> \[ \delta_{\mathcal{E}_\alpha} X_{mn}({\bf k}) = - \begin{cases}
    !>    \frac{\delta_{\mathcal{E}_\alpha} H_{mn}({\bf k})} 
    !>    {\epsilon_{m{\bf k}} - \epsilon_{n{\bf k}}} & \text{if } \epsilon_{m{\bf k}} \neq \epsilon_{n{\bf k}} \\
    !>    0 & \text{if } \epsilon_{m{\bf k}} = \epsilon_{n{\bf k}} 
    !>    \end{cases} \;, \]
    !> where \(\delta_{\mathcal{E}_\alpha} {\bf H}({\bf k})\) is the response of the Hamiltonian matrix as obtained from
    !> [[ef_eig_gen_dHmat(subroutine)]].
    !>
    !> If `projector=.true.`, then the eigenvector response corresponding to the projection of the
    !> wavefunction response onto the manifold of unoccupied states, 
    !> \(\hat{\bf P}_u \delta_{\mathcal{E}_\alpha} \psi_{n{\bf k}}({\bf r})\), is computed. This is 
    !> sufficient to compute the density response and ensures (given that \(n\) corresponds to
    !> an occupied state) that the energy differences \(\epsilon_{m{\bf k+q}} - \epsilon_{n{\bf k}}\)
    !> are always greater than zero.
    subroutine ef_eig_sternheimer( ik, Gkset, fst, lst, evalk, occk, eveck, dHmat, deval, devec, &
        projector, eps_deg )
      use constants, only: zzero, zone
      use mod_kpointset, only: Gk_set
      use mod_APW_LO, only: nlotot
      use mod_eigenvalue_occupancy, only: occmax
      use modinput
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k}\) vectors
      type(Gk_set), intent(in) :: Gkset
      !> first and last state for which the Sternheimer equation is solved
      integer, intent(in) :: fst, lst
      !> eigenvalues at \({\bf k}\)
      real(dp), intent(in) :: evalk(:)
      !> occupation numbers at \({\bf k}\)
      real(dp), intent(in) :: occk(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> Hamiltonian matrix response
      complex(dp), intent(in) :: dHmat(:,:)
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

      integer :: nmatk, ist, jst, nst
      real(dp) :: eps, dev, sig, t1
      logical :: proj

      complex(dp), allocatable :: dX(:,:)

      proj = .true.
      if( present( projector ) ) proj = projector
      eps = 1e-5_dp
      if( present( eps_deg ) ) eps = eps_deg

      nmatk = Gkset%ngk(1, ik) + nlotot
      nst = lst - fst + 1

      allocate( dX(nmatk, fst:lst) )

      deval = 0.0_dp
      do ist = fst, lst
        dX(:, ist) = dHmat(1:nmatk, ist)
        deval(ist) = dble( dX(ist, ist) )
        if( occk(ist) < input%groundstate%epsocc ) cycle
        do jst = 1, nmatk
          dev = evalk(ist) - evalk(jst)
          sig = sign( 1.0_dp, dev )
          dev = abs( dev )
          if( proj ) then
            if( dev < eps ) then
              t1 = 0.0_dp
            else
              t1 = sig * max( 0.0_dp, occmax - occk(jst) ) / (dev + eps)
            end if
            t1 = t1 / occmax
          else
            if( dev < eps ) then
              t1 = 0.0_dp
            else
              t1 = sig / dev
            end if
          end if
          dX(jst, ist) = t1 * dX(jst, ist)
        end do
      end do

      call zgemm( 'n', 'n', nmatk, nst, nmatk, zone, &
             eveck, size( eveck, dim=1 ), &
             dX, nmatk, zzero, &
             devec, size( devec, dim=1 ) )

      deallocate( dX )
    end subroutine ef_eig_sternheimer

    !> This subroutine computes the radial muffin-tin integrals times Gaunt coefficients
    !> for the momentum matrix.
    subroutine gen_momentum_mt_basis( Pmat_mt_basis, &
        lmax_apw )
      use constants, only: y00, zzero, zone
      use matrix_elements
      use mod_muffin_tin, only: nrmtmax
      use mod_atoms, only: nspecies, natoms, idxas
      !> momentum radial muffin-tin integrals times Gaunt coefficients
      complex(dp), allocatable, intent(inout) :: Pmat_mt_basis(:,:,:,:)
      !> maximum angular momentum \(l\) for APWs (default: from input file)
      integer, optional, intent(in) :: lmax_apw

      integer :: lmaxapw, is, ia, ias, ip
      real(dp), allocatable :: rfun(:,:)

      lmaxapw = dfpt_lmaxapw
      if( present( lmax_apw ) ) lmaxapw = lmax_apw

      ! allocate integrals
      call me_mt_alloc( Pmat_mt_basis, 3 )

      allocate( rfun(1, nrmtmax) )

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)

          rfun = 1.0_dp / y00

          do ip = 1, 3
            call me_mt_prepare( is, ias, 0, zone, rfun, zzero, Pmat_mt_basis(:, :, ias, ip), &
              right_gradient=ip )
          end do

        end do
      end do

      deallocate( rfun )
    end subroutine gen_momentum_mt_basis

    !> This subroutine computes the constant part of the Hamiltonian matrix response.
    !>
    !> I.e., it computes
    !> \[ -\langle \psi_{m{\bf k}} | r_\alpha | \psi_{n{\bf k}} \rangle 
    !>     = - \frac{ \langle \psi_{m{\bf k}} | [\hat{\bf H}, r_\alpha] | \psi_{n{\bf k}} \rangle }{\epsilon_{m{\bf k}} - \epsilon_{n{\bf k}}}
    !>     = {\rm i} \frac{ \langle \psi_{m{\bf k}} | p_\alpha | \psi_{n{\bf k}} \rangle }{\epsilon_{m{\bf k}} - \epsilon_{n{\bf k}}}
    !>     = \frac{ \langle \psi_{m{\bf k}} | \nabla_\alpha | \psi_{n{\bf k}} \rangle }{\epsilon_{m{\bf k}} - \epsilon_{n{\bf k}}}
    !>     \;. \]
    subroutine gen_dH0( ik, Gkset, fst1, lst1, fst2, lst2, evalk, eveck, apwalmk, ip, dHmat, &
        diagonal )
      use constants, only: zzero, zone
      use matrix_elements
      use mod_kpointset, only: Gk_set
      use mod_atoms, only: nspecies, natoms, idxas
      !> index of the \({\bf k}\) point
      integer, intent(in) :: ik
      !> set of \({\bf G+k}\) vectors
      type(Gk_set), intent(in) :: Gkset
      !> first and last state on the left for which the matrix elements are calculated
      integer, intent(in) :: fst1, lst1
      !> first and last state on the right for which the matrix elements are calculated
      integer, intent(in) :: fst2, lst2
      !> eigenvalues at \({\bf k}\)
      real(dp), intent(in) :: evalk(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) at \({\bf k}\)
      complex(dp), intent(in) :: apwalmk(:,:,:,:)
      !> polarization direction
      integer, intent(in) :: ip
      !> Hamiltonian response
      complex(dp), intent(out) :: dHmat(fst1:, fst2:)
      !> conmpute only `fst2` to `lst2` diagonal elements 
      !> and store them in the first dimension (default: `.false.`)
      logical, optional, intent(in) :: diagonal

      integer :: ngk, ist, jst
      integer :: is, ia, ias
      real(dp) :: dev
      logical :: diag

      diag = .false.
      if( present( diagonal ) ) diag = diagonal

      dHmat = zzero

      ngk = Gkset%ngk(1, ik)

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, ngk, ngk, apwalmk(:, :, :, ias), apwalmk(:, :, :, ias), &
                 eveck(:, fst1:lst1), eveck(:, fst2:lst2), &
                 zone, Pmat_mt_basis(:, :, ias, ip), zone, dHmat, &
                 diagonal_only=diag )
        end do
      end do
      call me_ir_mat( Gkset, ik, Gkset, ik, &
             eveck(:, fst1:lst1), eveck(:, fst2:lst2), &
             zone, cfun_ig, zone, dHmat, &
             right_gradient=ip, diagonal_only=diag )

      ! make matrix anti-Hermitian
      if( diag ) then
        dHmat = cmplx( 0, aimag( dHmat ), dp )
      else
        do jst = max( fst1, fst2 ), min( lst1, lst2 )
          do ist = jst, min( lst1, lst2 )
            dHmat(ist, jst) = (dHmat(ist, jst) - conjg( dHmat(jst, ist) )) / 2
            dHmat(jst, ist) = - conjg( dHmat(ist, jst) )
          end do
        end do
      end if

      ! divide by energy differences
      if( diag ) then
        dHmat = zzero
      else
        do jst = fst2, lst2
          do ist = fst1, lst1
            dev = evalk(ist) - evalk(jst)
            if( abs( dev ) < 1e-5_dp ) then
              dHmat(ist, jst) = zzero
            else
              dHmat(ist, jst) = dHmat(ist, jst) / dev
            end if
          end do
        end do
      end if
    end subroutine gen_dH0

    !> Rotates the eigenvector response corresponding to wavefunctions \(\delta_{\mathcal{E}_\alpha} \psi_{n{\bf p}}\)
    !> to the ones corresponing to \(\delta_{\mathcal{E}_\alpha} \psi_{n{\bf p}'}\) using the symmetry operation
    !> `isym` which rotates \({\bf p}\) into \({\bf p}'\), i.e., 
    !> \({\bf p}' = {\bf S}^\top \cdot {\bf p}\). 
    subroutine ef_eig_rotate_devec( isym, vpl, vprl, ngp, vgpl, vgprl, devec, ld1, ld2, nst )
      use constants, only: zzero, zone
      use mod_symmetry, only: symlatc, lsplsymc
      !> index of symmetry operation
      integer, intent(in) :: isym
      !> k-point \({\bf p}\) in lattice coordinates
      real(dp), intent(in) :: vpl(3)
      !> rotated k-point \({\bf p}' = {\bf S}^\top \cdot {\bf p}\) in lattice coordinates
      real(dp), intent(in) :: vprl(3)
      !> number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> \({\bf G+p}\) vectors at \({\bf p}\)
      real(dp), intent(in) :: vgpl(3, *)
      !> \(\bf G+p}'\) vectors at \({bf p}'\)
      real(dp), intent(in) :: vgprl(3, *)
      !> leading dimensions of `devec`
      integer, intent(in) :: ld1, ld2
      !> on input: eigenvector response at \({\bf p}\);
      !> on output: eigenvector response at \({\bf p}'\)
      complex(dp), intent(inout) :: devec(ld1, ld2, *)
      !> number of states
      integer, intent(in) :: nst

      integer :: ip
      complex(dp) :: zsymmat(3, 3)

      complex(dp), allocatable :: devec_tmp(:,:,:)

      allocate( devec_tmp, source=devec(:, 1:nst, 1:3) )
      zsymmat = cmplx( symlatc(:, :, lsplsymc(isym)), 0, dp )

      do ip = 1, 3
        call rotate_evecfv( isym, vpl, vprl, ngp, vgpl, vgprl, devec_tmp(:, :, ip), ld1, nst )
      end do
      call zgemm( 'n', 'c', ld1*nst, 3, 3, zone, &
             devec_tmp, ld1*nst, &
             zsymmat, 3, zzero, &
             devec, ld1*ld2 )

      deallocate( devec_tmp )
    end subroutine ef_eig_rotate_devec

end module efield_eigensystem
