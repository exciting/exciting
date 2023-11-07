!> This module contains procedures for the calculation of the response
!> of the macroscopic polarzation upon a perturbing electric field
!> \(\delta_{\mathcal{E}_\alpha}\).
module efield_polarization
  use dfpt_variables
  use dfpt_polarization, only: PWmat_mt_basis

  use precision, only: dp

  implicit none
  private

  public :: ef_pol_dpol_k, &
            ef_pol_symmetrize

  contains
    
  !> This subroutine calculates the contribution to the polarization 
  !> response coming from the given \({\bf k}\)-point.
  !>
  !> The calculation is based on the Berry phase approach to the macroscopic polarization
  !> outlined in 
  !> [King-Smith, Vanderbilt, *Phys. Rev. B* **47**, 1651-1654 (1993).](https://doi.org/10.1103/PhysRevB.47.1651)
  !>
  !> The berry phase along a 1D thread through the BZ is given by
  !> \[ \varphi({\bf k}_\parallel) = 2\, \Im \left[ \ln \prod_{j=0}^{J-1} 
  !> \det \langle u_{{\bf k}_j} | u_{{\bf k}_{j+1}} \rangle \right] \;. \]
  !> Using Jacobi's formula for the derivative of the determinant, we find for the variation of the
  !> Berry phase
  !> \[ \delta\varphi({\bf k}_\parallel) = 2\, \Im 
  !>    \frac{ \sum_{j=0}^{J-1} \delta\det \langle u_{{\bf k}_j} | u_{{\bf k}_{j+1}} \rangle 
  !>    \prod_{k\neq j}^{J-1} \det \langle u_{{\bf k}_j} | u_{{\bf k}_{j+1}} \rangle }
  !>    { \prod_{j=0}^{J-1} \det \langle u_{{\bf k}_j} | u_{{\bf k}_{j+1}} \rangle } 
  !>  = 2\, \Im \sum_{j=0}^{J-1} \operatorname{tr} \left[ \langle u_{{\bf k}_j} | u_{{\bf k}_{j+1}} \rangle^{-1}\,
  !>    \delta \langle u_{{\bf k}_j} | u_{{\bf k}_{j+1}} \rangle \right] \;. \]
  !> The overlap between cell-periodic parts of the wavefunctions are given by the plane wave matrix elements
  !> \[ M_{mn}({\bf k}_j) = \langle u_{m{\bf k}_j} | u_{n{\bf k}_{j+1}} \rangle 
  !>    = \langle \psi_{m{\bf k}_j} | \rm{e}^{-{\rm i} ({\bf k}_{j+1} - {\bf k}_j) \cdot {\bf r}} | \psi_{n{\bf k}_{j+1}} \rangle \;. \]
  subroutine ef_pol_dpol_k( iknr, ikbnr, kset, Gkset, fst, lst, eveck, eveckb, deveck, deveckb, dpol )
    use dfpt_eigensystem, only: cfun_ig

    use constants, only: zzero, zone, zi
    use matrix_elements
    use matrix_elements_lapw_lo, only: me_lapwlo_ir_mat
    use mod_kpointset, only: k_set, Gk_set
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_Gkvector, only: ngkmax_ptr
    use mod_muffin_tin, only: lmmaxapw
    use mod_APW_LO, only: apwordmax
    use mod_eigenvalue_occupancy, only: occmax
    !> index of the \({\bf k}\) point in full BZ
    integer, intent(in) :: iknr
    !> index of the \({\bf k+b}\) point in full BZ
    integer, intent(in) :: ikbnr
    !> set of \({\bf k}\) vectors
    type(k_set), intent(in) :: kset
    !> set of \({\bf G+k}\) vectors
    type(Gk_set), intent(in) :: Gkset
    !> first and last state to consider in the sum over states
    integer, intent(in) :: fst, lst
    !> eigenvectors at \({\bf k}\)
    complex(dp), intent(in) :: eveck(:,:)
    !> eigenvectors at \({\bf k+b}\)
    complex(dp), intent(in) :: eveckb(:,:)
    !> eigenvector response at \({\bf k}\)
    complex(dp), intent(in) :: deveck(:,:)
    !> eigenvector response at \({\bf k+b}\)
    complex(dp), intent(in) :: deveckb(:,:)
    !> polarization response
    complex(dp), intent(inout) :: dpol

    integer :: nst, ngk, ngkb, shift(3), info, ist
    integer :: is, ia, ias
    integer, target :: ngkmax
    real(dp) :: vbl(3), vbc(3)
    complex(dp) :: z1

    integer, allocatable :: ipiv(:)
    complex(dp), allocatable :: apwalmk(:,:,:,:), apwalmkb(:,:,:,:), M(:,:), dM(:,:)

    nst = lst - fst + 1
    ngk = Gkset%ngknr(1, iknr)
    ngkb = Gkset%ngknr(1, ikbnr)
    ngkmax = Gkset%ngkmax
    ngkmax_ptr => ngkmax

    allocate( apwalmk(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
    allocate( apwalmkb(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
    allocate( M(nst, nst), source=zzero )
    allocate( dM(nst, nst), source=zzero )
    allocate( ipiv(nst) )

    ! get b vector and G shift
    vbl = kset%vklnr(:, ikbnr) - kset%vklnr(:, iknr)
    call r3frac( 1e-6_dp, vbl, shift )
    call r3mv( kset%bvec, vbl, vbc )
    shift = - shift
    ! get matching coefficients at k
    call match( ngk, Gkset%gknrc(:, 1, iknr), Gkset%tpgknrc(:, :, 1, iknr), Gkset%sfacgknr(:, :, 1, iknr), &
      apwalmk )
    ! get matching coefficients at k+b
    call match( ngkb, Gkset%gknrc(:, 1, ikbnr), Gkset%tpgknrc(:, :, 1, ikbnr), Gkset%sfacgknr(:, :, 1, ikbnr), &
      apwalmkb )

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        ! M = <u_k|u_k+b>
        call me_mt_mat( is, ias, ngk, ngkb, apwalmk(:, :, :, ias), apwalmkb(:, :, :, ias), &
          eveck(:, fst:lst), eveckb(:, fst:lst), &
          zone, PWmat_mt_basis(:, :, ias), zone, M )
        ! dM = d<u_k|u_k+b>
        call me_mt_mat( is, ias, ngk, ngkb, apwalmk(:, :, :, ias), apwalmkb(:, :, :, ias), &
          deveck(:, fst:lst), eveckb(:, fst:lst), &
          zone, PWmat_mt_basis(:, :, ias), zone, dM )
        call me_mt_mat( is, ias, ngk, ngkb, apwalmk(:, :, :, ias), apwalmkb(:, :, :, ias), &
          eveck(:, fst:lst), deveckb(:, fst:lst), &
          zone, PWmat_mt_basis(:, :, ias), zone, dM )
      end do
    end do
    ! M = <u_k|u_k+b>
    call me_lapwlo_ir_mat( Gkset, iknr, Gkset, ikbnr, dfpt_Gset, &
      zone, cfun_ig, zone, M, &
      left_evec=eveck(:, fst:lst), right_evec=eveckb(:, fst:lst), &
      G_shift=shift, non_reduced_p1=.true., non_reduced_p2=.true. )
    ! dM = d<u_k|u_k+b>
    call me_lapwlo_ir_mat( Gkset, iknr, Gkset, ikbnr, dfpt_Gset, &
      zone, cfun_ig, zone, dM, &
      left_evec=deveck(:, fst:lst), right_evec=eveckb(:, fst:lst), &
      G_shift=shift, non_reduced_p1=.true., non_reduced_p2=.true. )
    call me_lapwlo_ir_mat( Gkset, iknr, Gkset, ikbnr, dfpt_Gset, &
      zone, cfun_ig, zone, dM, &
      left_evec=eveck(:, fst:lst), right_evec=deveckb(:, fst:lst), &
      G_shift=shift, non_reduced_p1=.true., non_reduced_p2=.true. )

    ! solve linear system M.S = dM to find S = M^-1.dM
    call zgesv( nst, nst, M, nst, ipiv, dM, nst, info )

    ! sum up trace
    z1 = zi * kset%wkptnr(iknr) * occmax
    do ist = 1, nst
      dpol = dpol - z1 * dM(ist, ist)
    end do

    deallocate( apwalmk, apwalmkb, M, dM, ipiv )
  end subroutine ef_pol_dpol_k

  !> Symmetrize polarization response.
  subroutine ef_pol_symmetrize( dpol, nsym, isym )
    use constants, only: zzero
    use mod_symmetry, only: symlatc, lsplsymc
    !> polarization response
    complex(dp), intent(inout) :: dpol(3,3)
    !> number of symmetry operations
    integer, intent(in) :: nsym
    !> indices of symmetries in global arrays
    integer, intent(in) :: isym(:)

    integer :: s, lspl, ip, ipp
    real(dp) :: sc(3,3), a(3), b(3), c(3)
    complex(dp) :: dpc(3,3)

    ! make a copy of the input polarization response and nullify it
    dpc = dpol
    dpol = zzero

    ! loop over symmetries
    do s = 1, nsym
      lspl = lsplsymc(isym(s))
      sc = symlatc(:, :, lspl)
      do ip = 1, 3
        ! apply rotation to polarization response
        a = dble(  dpc(:, ip) )
        b = aimag( dpc(:, ip) )
        call r3mtv( sc, a, c ); a = c
        call r3mtv( sc, b, c ); b = c
        ! add contribution to input array
        do ipp = 1, 3
          dpol(:, ipp) = dpol(:, ipp) + sc(ip, ipp) * cmplx( a, b, dp )
        end do
      end do
    end do

    ! normalize
    dpol = dpol / nsym
  end subroutine ef_pol_symmetrize

end module efield_polarization
