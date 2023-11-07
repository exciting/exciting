!> This module contains procedures for the calculation of the response
!> of the macroscopic polarzation upon a phonon-like
!> perturbation \(\delta^{\bf q}_{I \mu}\).
module phonons_polarization
  use dfpt_variables
  use dfpt_polarization
  use phonons_variables

  use precision, only: dp

  implicit none
  private

  !> variation of interstitial characteristic function for all irrep members
  complex(dp), allocatable :: dcfun_ig(:,:)

  public :: ph_pol_init, ph_pol_free, &
            ph_pol_dpol_k, ph_pol_dpol_ion, &
            ph_pol_symmetrize

  contains

    !> Initialize plane wave matrix elements for a given irrep.
    subroutine ph_pol_init( irr )
      use phonons_symmetry, only: irrep
      use phonons_eigensystem, only: gen_dcfun_ig

      use constants, only: zzero
      use mod_atoms, only: nspecies, natoms, idxas
      !> irrep
      type(irrep), intent(in) :: irr

      integer :: id, is, ia, ias

      allocate( dcfun_ig(ph_Gqset%ngvec, irr%dim), source=zzero )

      do id = 1, irr%dim
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            call gen_dcfun_ig( ph_Gqset%ngvec, ph_Gqset%gc, ph_Gqset%vgc, is, ia, irr%pat(:, ias, id), dcfun_ig(:, id) )
          end do
        end do
      end do

    end subroutine ph_pol_init

    !> Free memory from unneeded variables.
    subroutine ph_pol_free
      call dfpt_pol_free
      if( allocated( dcfun_ig ) ) deallocate( dcfun_ig )
    end subroutine ph_pol_free

    !> This subroutine calculates the contribution to the polarization 
    !> response coming from the given \({\bf k}\)-point.
    !>
    !> The calulcation is based on the Berry phase approach to the macroscopic polarizatiob
    !> outlined in 
    !> [King-Smith, Vanderbilt, *Phys. Rev. B* **47**, 1651-1654 (1993).](https://doi.org/10.1103/PhysRevB.47.1651)
    !>
    !> Ther berry phase along a 1D thread through the BZ is given by
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
    subroutine ph_pol_dpol_k( iknr, ikbnr, kset, Gkset, fst, lst, eveck, eveckb, deveck, deveckb, pat, id, dpol )
      use dfpt_eigensystem, only: cfun_ig
      use phonons_eigensystem, only: gen_dapwalm

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
      !> displacement pattern \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> index of irrep member
      integer, intent(in) :: id
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

      allocate( apwalmk(ngkmax_ptr, apwordmax, lmmaxapw, 0:natmtot) )
      allocate( apwalmkb(ngkmax_ptr, apwordmax, lmmaxapw, 0:natmtot) )
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
        apwalmk(:, :, :, 1:natmtot) )
      ! get matching coefficients at k+b
      call match( ngkb, Gkset%gknrc(:, 1, ikbnr), Gkset%tpgknrc(:, :, 1, ikbnr), Gkset%sfacgknr(:, :, 1, ikbnr), &
        apwalmkb(:, :, :, 1:natmtot) )

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! get matching coefficient response at k
          call gen_dapwalm( Gkset%ngknr(1, iknr), Gkset%vgknrc(:, :, 1, iknr), conjg( pat(:, ias) ), apwalmk(:, :, :, ias), &
            apwalmk(:, :, :, 0) )
          ! get matching coefficient response at k+b
          call gen_dapwalm( Gkset%ngknr(1, ikbnr), Gkset%vgknrc(:, :, 1, ikbnr), pat(:, ias), apwalmkb(:, :, :, ias), &
            apwalmkb(:, :, :, 0) )
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
          call me_mt_mat( is, ias, ngk, ngkb, apwalmk(:, :, :, 0), apwalmkb(:, :, :, ias), &
            eveck(:, fst:lst), eveckb(:, fst:lst), &
            zone, PWmat_mt_basis(:, :, ias), zone, dM, &
            left_local_orbitals=.false. )
          call me_mt_mat( is, ias, ngk, ngkb, apwalmk(:, :, :, ias), apwalmkb(:, :, :, 0), &
            eveck(:, fst:lst), eveckb(:, fst:lst), &
            zone, PWmat_mt_basis(:, :, ias), zone, dM, &
            right_local_orbitals=.false. )
          z1 = - zi * dot_product( vbc, pat(:, ias) )
          call me_mt_mat( is, ias, ngk, ngkb, apwalmk(:, :, :, ias), apwalmkb(:, :, :, ias), &
            eveck(:, fst:lst), eveckb(:, fst:lst), &
            z1, PWmat_mt_basis(:, :, ias), zone, dM )
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
      call me_lapwlo_ir_mat( Gkset, iknr, Gkset, ikbnr, ph_Gqset, &
        zone, dcfun_ig(:, id), zone, dM, &
        left_evec=eveck(:, fst:lst), right_evec=eveckb(:, fst:lst), &
        G_shift=shift, non_reduced_p1=.true., non_reduced_p2=.true. )

      ! solve linear system M.S = dM to find S = M^-1.dM
      call zgesv( nst, nst, M, nst, ipiv, dM, nst, info )

      ! sum up trace
      z1 = zi * kset%wkptnr(iknr) * occmax
      do ist = 1, nst
        dpol = dpol - z1 * dM(ist, ist)
      end do

      deallocate( apwalmk, apwalmkb, M, dM, ipiv )
    end subroutine ph_pol_dpol_k

    !> Add ionic controbution to polarization response.
    subroutine ph_pol_dpol_ion( dirrep, pat, dpol )
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spze, spnst, spocc, spcore
      !> irrep dimension
      integer, intent(in) :: dirrep
      !> displacement pattern \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot, *)
      !> polarization response
      complex(dp), intent(inout) :: dpol(3, *)

      integer :: is, ia, ias, id, ist
      real(dp) :: zion
      
      do is = 1, nspecies
        zion = spze(is)
        do ist = 1, spnst(is)
          if( spcore(ist, is) ) zion = zion - spocc(ist, is)
        end do
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do id = 1, dirrep
            dpol(:, id) = dpol(:, id) + zion * pat(:, ias, id)
          end do
        end do
      end do
    end subroutine ph_pol_dpol_ion

    !> Symmetrize polarization response in basis of an irrep.
    subroutine ph_pol_symmetrize( dpol, dirrep, nsym, isym, ivsym, symmat )
      use constants, only: zzero
      use mod_symmetry, only: symlatc, lsplsymc
      !> polarization response
      complex(dp), intent(inout) :: dpol(3, *)
      !> dimension of the irrep
      integer, intent(in) :: dirrep
      !> number of symmetry operations in small group of q
      integer, intent(in) :: nsym
      !> indices of symmetries in global arrays
      integer, intent(in) :: isym(:)
      !> lattice vectors that map Sq back to 1st BZ
      integer, intent(in) :: ivsym(3, *)
      !> matrix representation of symmetries in the basis of the irrep
      complex(dp), intent(in) :: symmat(:,:,:)
      
      integer :: s, lspl, d, dd
      real(dp) :: sc(3,3), a(3), b(3), c(3)

      complex(dp), allocatable :: dpc(:,:)

      ! make a copy of the input polarization response and nullify it
      allocate( dpc, source=dpol(:, 1:dirrep) )
      dpol(:, 1:dirrep) = zzero

      ! loop over symmetries
      do s = 1, nsym
        lspl = lsplsymc(isym(s))
        sc = symlatc(:, :, lspl)
        do d = 1, dirrep
          ! apply rotation to polarization response
          a = dble(  dpc(:, d) )
          b = aimag( dpc(:, d) )
          call r3mtv( sc, a, c ); a = c
          call r3mtv( sc, b, c ); b = c
          ! add contribution to input array
          do dd = 1, dirrep
            dpol(:, dd) = dpol(:, dd) + symmat(d, dd, s) * cmplx( a, b, dp )
          end do
        end do
      end do

      ! normalize
      do d = 1, dirrep
        dpol(:, d) = dpol(:, d) / nsym
      end do

      deallocate( dpc )
    end subroutine ph_pol_symmetrize

end module phonons_polarization
