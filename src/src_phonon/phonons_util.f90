module phonons_util
  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: terminate_if_false

  implicit none

  contains

    !> Impose the acoustic sum rule on a set of dynamical matrices \({\bf D}({\bf p})\).
    !>
    !> On exit, the dynamical matrices are scaled such that the three smallest
    !> eigenvalues \({\bf p}=\Gamma\) are zero. \(\Gamma\) must be part of point set.
    subroutine ph_util_sumrule_dyn( vpl, dynp )
      use mod_symmetry, only: symlat, lsplsymc, nsymcrys, find_equivalent_wavevectors
      use math_utils, only: is_square
      use m_linalg, only: zhediag
      !> list of reciprocal space points \({\bf p}\) in lattice coordinates
      real(dp), intent(in) :: vpl(:, :)
      !> dynamical matrices \({\bf D}({\bf p})\)
      complex(dp), intent(inout) :: dynp(:,:,:)

      integer :: np, ip, ip0, n, i, j, k
      
      integer, allocatable :: ip_idx(:), isym_idx(:)
      real(dp), allocatable :: eval(:)
      complex(dp), allocatable :: dyn0(:,:), evec(:,:)

      call assert( size( vpl, dim=1 ) == 3, 'Wavevectors are not 3-dimensional.' )
      np = size( vpl, dim=2 )
      call assert( is_square( dynp(:,:,1) ), 'Dynamical matrices are not square.' )
      call assert( size( dynp, dim=3 ) == np, 'Different number of dynamical matrices and wavevectors.' )
      n = size( dynp, dim=1 )

      allocate( eval(n), dyn0(n,n), evec(n,n) )

      ! find Gamma point in list
      call find_equivalent_wavevectors( 3, [0.0_dp, 0.0_dp, 0.0_dp], vpl, np, symlat(:, :, lsplsymc(1:nsymcrys)), nsymcrys, &
        ip_idx, isym_idx, first_only=.true. )
      call terminate_if_false( size( ip_idx ) > 0, '(ph_util_sumrule_dyn) &
        Gamma point not found in list of wavevectors.' )
      ip0 = ip_idx(1)
      dyn0 = 0.5d0 * (dynp(:, :, ip0) + conjg( transpose( dynp(:, :, ip0) ) ))
      call zhediag( dyn0, eval, evec=evec )

      do ip = 1, np
        do i = 1, n
          do j = 1, n
            do k = 1, 3
              dynp(i, j, ip) = dynp(i, j, ip) - eval(k) * evec(i, k) * conjg( evec(j, k) )
            end do
          end do
        end do
      end do

      deallocate( eval, dyn0, evec, ip_idx, isym_idx )
    end subroutine ph_util_sumrule_dyn

    !> Impose acoustic sum rule on a set of Born effective charge tensos \({\bf Z}^\ast_\kappa\).
    !>
    !> On exit, the charge tensors sum up to zero.
    subroutine ph_util_sumrule_borncharge( zstar, correction )
      use mod_atoms, only: natmtot
      !> Born effective charge tensors \({\bf Z}^\ast_\kappa\)
      real(dp), intent(inout) :: zstar(3, 3, natmtot)
      !> correction that was applied to all atoms to impose sum rule
      real(dp), intent(out) :: correction(3, 3)

      integer :: ias

      correction = sum( zstar, dim=3 ) / natmtot
      do ias = 1, natmtot
        zstar(:, :, ias) = zstar(:, :, ias) - correction
      end do
    end subroutine ph_util_sumrule_borncharge

    !> Symmetrize a dynamical matrix \({\bf D}({\bf p}\) with a given
    !> set of crystal symmetries.
    !> See also [[ph_util_symapp_dyn(subroutine)]].
    subroutine ph_util_symmetrize_dyn( vpl, dyn, nsym, isym )
      use constants, only: zzero
      use math_utils, only: is_square
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> dynamical matrix \({\bf D}({\bf p})\)
      complex(dp), intent(inout) :: dyn(:,:)
      !> number of crystal symmetries
      integer, intent(in) :: nsym
      !> indices of crystal symmetries in global symmetry arrays
      integer, intent(in) :: isym(:)

      integer :: n, jsym
      complex(dp), allocatable :: tmp(:,:)

      call assert( is_square(dyn), 'Dynamical matrix is not square.' )
      n = size( dyn, dim=1 )

      ! make Hermitian
      dyn = (dyn + conjg( transpose( dyn ) )) / 2

      ! symmetrize
      allocate( tmp(n, n), source=zzero )
      do jsym = 1, nsym
        call ph_util_symapp_dyn( isym(jsym), vpl, dyn, tmp, dir=1, matrix='T' )
      end do
      tmp = tmp / nsym

      dyn = zzero
      do jsym = 1, nsym
        call ph_util_symapp_dyn( isym(jsym), vpl, tmp, dyn, dir=-1, matrix='T' )
      end do
      dyn = dyn / nsym

      deallocate( tmp )
    end subroutine ph_util_symmetrize_dyn

    !> Apply a symmetry operation to a dynamical matrix \({\bf D}({\bf p})\).
    !>
    !> This is done by applying the sandwich product
    !> \[ {\bf S}({\bf p})\, {\bf D}({\bf p})\, {\bf S}^\dagger({\bf p}) \;, \]
    !> where \({\bf S}({\bf p})\) is either the symmetry matrix \({\bf \Gamma}({\bf p})\) or
    !> \({\bf T}({\bf p})\) according to *Maradudin, Vosko, Rev. Mod. Phys. **40**, 1 (1968)*
    !> if `dir >= 0` or the respective Hermitian conjugate if `dir < 0`.
    !> The result is added to the input argument `dyns`.
    !> See also [[ph_util_symmetry_G(subroutine)]] and [[ph_util_symmetry_T(subroutine)]].
    subroutine ph_util_symapp_dyn( isym, vpl, dyn, dyns, &
        dir, matrix )
      use constants, only: zzero, zone
      use mod_atoms, only: natmtot
      !> index of symmetry operation in global symmetry arrays
      integer, intent(in) :: isym
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> dynamical matrix \({\bf D}({\bf p})\)
      complex(dp), intent(in) :: dyn(:,:)
      !> symmetrized dynamical matrix
      complex(dp), intent(inout) :: dyns(:,:)
      !> direction of symmetry application (default: `1`)
      integer, optional, intent(in) :: dir
      !> type of symmetry matrix 
      !>(`G` for \({\bf \Gamma}\), `T` for \({\bf T}\); default: `G`)
      character, optional, intent(in) :: matrix

      integer :: direction
      character :: m, c1, c2
      complex(dp), allocatable :: aux(:,:), S(:,:)

      direction = 1
      if( present( dir ) ) direction = dir
      if( present( matrix ) ) m = matrix
      if( m /= 'G' .or. m /= 'T') m = 'G'

      allocate( aux(3*natmtot, 3*natmtot) )

      if( m == 'T') then
        S = ph_util_symmetry_T( isym, vpl )
      else
        S = ph_util_symmetry_G( isym, vpl )
      end if
      if( direction >= 0 ) then
        c1 = 'n'
        c2 = 'c'
      else
        c1 = 'c'
        c2 = 'n'
      end if
      call zgemm( c1, 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             S, 3*natmtot, &
             dyn, size( dyn, dim=1 ), zzero, &
             aux, 3*natmtot )
      call zgemm( 'n', c2, 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             aux, 3*natmtot, &
             S, 3*natmtot, zone, &
             dyns, size( dyns, dim=1 ) )

      deallocate( S, aux)
    end subroutine ph_util_symapp_dyn

    !> Get the symmetry matrix \({\bf \Gamma}({\bf p})\) for a given
    !> symmetry opertation
    !>
    !> according to *Maradudin, Vosko, Rev. Mod. Phys. **40**, 1 (1968)*.
    !> Note: In exciting, the translation is applied before the rotation!
    function ph_util_symmetry_G( isym, vpl ) result( G )
      use constants, only: zzero, zone, twopi
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use math_utils, only: all_zero
      use mod_symmetry, only: lsplsymc, isymlat, symlat, symlatc, ieqatom, vtlsymc
      use modinput
      !> index of symmetry operation in global symmetry arrays
      integer, intent(in) :: isym
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> symemtry matrix
      complex(dp), allocatable :: G(:,:)

      integer :: lspl, ilspl, is, ia, ias, ja, jas, i, j
      real(dp) :: sl(3, 3), sc(3, 3), v1(3), v2(3), phase
      complex(dp) :: aux(3*natmtot, 3*natmtot)
      character(256) :: errmsg

      if( allocated( G ) ) deallocate( G )
      allocate( G(3*natmtot, 3*natmtot), source=zzero )

      lspl = lsplsymc(isym)
      ilspl = isymlat(lspl)
      sl = dble( symlat(:, :, ilspl) )
      sc = symlatc(:, :, lspl)
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          i = (ias - 1) * 3 + 1
          ja = ieqatom(ia, is, isym)
          jas = idxas(ja, is)
          j = (jas - 1) * 3 + 1
          v1 = input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord + vtlsymc(:, isym)
          call r3mv( sl, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord, v2 )
          phase = twopi * dot_product( vpl, v2 - v1 )
          G(i:i+2, j:j+2) = cmplx( cos( phase ), sin( phase ), dp ) * sc
        end do
      end do
      ! check if Gamma is unitary
      call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             G, 3*natmtot, &
             G, 3*natmtot, zzero, &
             aux, 3*natmtot )
      do ias = 1, 3*natmtot
        aux(ias, ias) = aux(ias, ias) - zone
      end do
      write( errmsg, '("(ph_util_symmetry_G) Non-unitary Gamma for symmetry ",i2," and q = ",3f13.6)' ) isym, vpl
      call terminate_if_false( all_zero( aux, tol=input%structure%epslat ), trim (errmsg ) )
    end function ph_util_symmetry_G

    !> Get the symmetry matrix \({\bf T}({\bf p})\) for a given
    !> symmetry opertation
    !>
    !> according to *Maradudin, Vosko, Rev. Mod. Phys. **40**, 1 (1968)*.
    !> Note: In exciting, the translation is applied before the rotation!
    function ph_util_symmetry_T( isym, vpl ) result( T )
      use constants, only: zzero, zone, twopi
      use math_utils, only: all_zero
      use mod_symmetry, only: lsplsymc, symlat, vtlsymc
      use mod_atoms, only: natmtot
      use modinput
      !> index of symmetry operation in global symmetry arrays
      integer, intent(in) :: isym
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> symemtry matrix
      complex(dp), allocatable :: T(:,:)

      integer :: lspl, ias
      real(dp) :: sl(3, 3), v1(3), phase
      complex(dp) :: aux(3*natmtot, 3*natmtot)
      character(256) :: errmsg

      T = ph_util_symmetry_G( isym, vpl )
      lspl = lsplsymc(isym)
      sl = dble( symlat(:, :, lspl) )
      call r3mv( sl, vtlsymc(:, isym), v1 )
      phase = twopi * dot_product( vpl, v1 )
      T = cmplx( cos( phase ), sin( phase ), dp ) * T
      ! check if T is unitary
      call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             T, 3*natmtot, &
             T, 3*natmtot, zzero, &
             aux, 3*natmtot )
      do ias = 1, 3*natmtot
        aux(ias, ias) = aux(ias, ias) - zone
      end do
      write( errmsg, '("(ph_util_symmetry_T): Non-unitary T for symmetry ",i2," and q = ",3f13.6)' ) isym, vpl
      call terminate_if_false( all_zero( aux, tol=input%structure%epslat ), trim( errmsg ))
    end function ph_util_symmetry_T

    !> Diagonalize dynamical matrix and obtain phonon frequencies
    !> and eigenvectors.
    !>
    !> Given a dynamical matrix \(D_{\kappa\alpha, \kappa'\beta}\), this calculates
    !> the eigenvalues \(\omega_\nu^2\) and eigenvectors \(e_{\kappa\alpha,\nu}\)
    !> such that 
    !> \[ {\bf M}\, {\bf D}\, {\bf M}
    !>    = {\bf e}\, \operatorname{diag}(\omega_\nu^2)\, {\bf e}^\dagger \;, \]
    !> where \({\bf M} = {\bf M}^0\) if `basis_transform` is not specified
    !> \({\bf M} = {\bf B}^\dagger \, {\bf M}^0\, {\bf B}\) otherwise. Here,
    !> \(M^0_{\kappa\alpha,\kappa'\beta} = 1/\sqrt{M_\kappa}\, \delta_{\kappa \kappa'}\, \delta_{\alpha\beta}\)
    !> and \(B_{\kappa\alpha, n}\) describes a basis transform from the Cartesian
    !> \((\kappa,\alpha)\) basis into another basis \(n\).
    !>
    !> If `normalized_evec = .false.`, \({\bf u} = {\bf M}\, {\bf e}\) is returned
    !> as the eigenvectors.
    subroutine ph_util_diag_dynmat( dyn, w, evec, &
        basis_transform, normalized_evec )
      use constants, only: zone, zzero
      use mod_atoms, only: natmtot, nspecies, natoms, spmass
      use m_linalg, only: zhediag, zhegauge
      !> dynamical matrix \({\bf D}\)
      complex(dp), intent(in) :: dyn(:,:)
      !> phonon freuqncies \(\omega_\nu\)
      real(dp), intent(out) :: w(:)
      !> phonon eigenvectors \({\bf e}_\nu\)
      complex(dp), intent(out) :: evec(:,:)
      !> unitary basis transformation (default: none)
      complex(dp), optional, intent(in) :: basis_transform(:,:)
      !> return normalized eigenvectors (default: `.true.`)
      logical, optional, intent(in) :: normalized_evec

      integer :: i, is, ia, ip
      logical :: norm

      complex(dp), allocatable :: mass_matrix(:, :), tmp(:, :)

      norm = .true.
      if( present( normalized_evec ) ) norm = normalized_evec

      allocate( mass_matrix(3*natmtot, 3*natmtot), tmp(3*natmtot, 3*natmtot) )

      ! set up mass matrix
      if( .not. present( basis_transform ) ) then
        mass_matrix = zzero
        i = 0
        do is = 1, nspecies
          do ia = 1, natoms(is)
            do ip = 1, 3
              i = i + 1
              if( spmass(is) > 0.0_dp ) &
                mass_matrix(i, i) = 1.0_dp / sqrt( spmass(is) )
            end do
          end do
        end do
      else
        tmp = zzero
        i = 0
        do is = 1, nspecies
          do ia = 1, natoms(is)
            do ip = 1, 3
              i = i + 1
              if( spmass(is) > 0.0_dp ) &
                tmp(i, :) = basis_transform(i, :) / sqrt( spmass(is) )
            end do
          end do
        end do
        call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
               basis_transform, 3*natmtot, &
               tmp, 3*natmtot, zzero, &
               mass_matrix, 3*natmtot )
      end if
      ! set up dynamical matrix (add mass term)
      tmp = 0.5_dp * (dyn + conjg( transpose( dyn ) ))
      call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             tmp, 3*natmtot, &
             mass_matrix, 3*natmtot, zzero, &
             evec, 3*natmtot )
      call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             mass_matrix, 3*natmtot, &
             evec, 3*natmtot, zzero, &
             tmp, 3*natmtot )
      ! diagonalize
      call zhediag( tmp, w, evec )
      ! assign phonon frequencies
      do i = 1, 3*natmtot
        if( w(i) .ge. 0.0_dp ) then
          w(i) = sqrt( w(i))
        else
          w(i) = -sqrt( abs( w(i) ) )
        end if
      end do
      ! fix unique gauge of degenerate eigenvectors
      call zhegauge( w, evec, 1e-12_dp )

      if( .not. norm ) then
        tmp = evec
        call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
               mass_matrix, 3*natmtot, &
               tmp, 3*natmtot, zzero, &
               evec, 3*natmtot )
      end if

      deallocate( mass_matrix, tmp )
    end subroutine ph_util_diag_dynmat

    !> Setup a phonon Fourier interpolation.
    !>
    !> This tries to read the dynamical matrices \(D({\bf q})\) from file, unless they are
    !> explicitly given via `dynmatq`. If both \({\bf \epsilon}^\infty\) (`dielten`) and 
    !> \({\bf Z}^\ast_\kappa\) (`borncharge`) is given, the non-analytic part is substracted
    !> from \(D({\bf q})\) and needs to be added back again after the interpolation.
    !>
    !> This routine returns a matrix Fourier interpolation object `mfi` and the dynamical 
    !> matrices in real space `dynmatr` (interatomic force constants), \(D({\bf R})\). 
    !> To interpolate back to reciprocal space \(D({\bf p})\), use [[ph_util_interpolate(subroutine)]].
    subroutine ph_util_setup_interpolation( bvec, ngridq, nq, ivq, vql, mfi, dynmatr, &
        sumrule, symmetrize, dielten, borncharge, directory, dynmatq )
      use constants, only: zzero, zone
      use phonons_io_util, only: ph_io_read_dynmat_grid
      use matrix_fourier_interpolation, only: mfi_type
      use mod_atoms, only: natmtot, nspecies, natoms, atposc
      use mod_symmetry, only: symlat, lsplsymc, nsymcrys, find_equivalent_wavevectors
      !> reciprical lattice vectors (column wise)
      real(dp), intent(in) :: bvec(3, 3)
      !> phonon \({\bf q}\)-grid division
      integer, intent(in) :: ngridq(3)
      !> number of (reduced) \({\bf q}\)-points on grid
      integer, intent(in) :: nq
      !> integer positions of \({\bf q}\)-points on grid
      integer, intent(in) :: ivq(3, nq)
      !> lattice coordinates of \({\bf q}\)-points
      real(dp), intent(in) :: vql(3, nq)
      !> matrix Fourier interpolation object
      type(mfi_type), intent(out) :: mfi
      !> dynamical matrices in real space (interatomic force constants)
      complex(dp), allocatable, intent(out) :: dynmatr(:, :, :)
      
      !> impose acoustic sumrule (default: `.true.`)
      logical, optional, intent(in) :: sumrule
      !> symmetrize dynamical matrices (default: `.true.` )
      logical, optional, intent(in) :: symmetrize
      !> high-frequency dielectric tensor \({\bf \epsilon}^\infty\)
      real(dp), optional, intent(in) :: dielten(3, 3)
      !> Born effective charges \({\bf Z}^\ast_\kappa\)
      real(dp), optional, intent(in) :: borncharge(3, 3, natmtot)
      !> directory to search for dynamical matrices in (default: `'./'`)
      character(*), optional, intent(in) :: directory
      !> dynamical matrices on \({\bf q}\)-grid
      complex(dp), optional, intent(in) :: dynmatq(3*natmtot, 3*natmtot, nq)

      integer :: iq, ias, i, j, k
      real(dp) :: vql_off(3)
      logical :: smrl, symm, success

      integer, allocatable  :: iq_idx(:), isym_idx(:)
      real(dp), allocatable :: vql_full(:, :)
      complex(dp), allocatable :: dynq(:, :, :), dynq_full(:, :, :)

      smrl = .true.
      if( present( sumrule ) ) smrl = sumrule
      symm = .true.
      if( present( symmetrize) ) symm = symmetrize

      ! check if dynamical matrices were given, otherwise try to read from file
      if( present( dynmatq ) ) then
        allocate( dynq, source=dynmatq )
      else
        allocate( dynq(3*natmtot, 3*natmtot, nq) )
        if( present( directory ) ) then
          call ph_io_read_dynmat_grid( ivq, ngridq, nspecies, natoms, dynq, success, directory=directory )
        else
          call ph_io_read_dynmat_grid( ivq, ngridq, nspecies, natoms, dynq, success )
        end if
        call terminate_if_false( success, '(ph_util_setup_interpolation): &
          failed to read dynamical matrices from file.' )
      end if

      ! symmetrize dynamical matrices
      if( symm ) then
        do iq = 1, nq
          call find_equivalent_wavevectors( 3, vql(:, iq), vql(:, iq:iq), 1, symlat(:, :, lsplsymc(1:nsymcrys)), nsymcrys, iq_idx, isym_idx )
          call ph_util_symmetrize_dyn( vql(:, iq), dynq(:, :, iq), size( isym_idx ), isym_idx )
        end do
      end if

      ! apply accoustic sum rule
      if( smrl ) &
        call ph_util_sumrule_dyn( vql(:, 1:nq), dynq )

      ! subtract non-analytic part for polar materials
      if( present( dielten ) .and. present( borncharge ) ) &
        call ph_util_dynmat_lr( bvec, nq, vql, dielten, borncharge, -zone, dynq )

      ! set up interpolation object
      vql_off = vql(:, 1) * ngridq
      vql_off = vql_off - floor( vql_off )
      vql_full = reshape( [(((([i, j, k] + vql_off) / ngridq, i=0, ngridq(1)-1), j=0, ngridq(2)-1), k=0, ngridq(3)-1)], [3, product(ngridq)] )
      mfi = mfi_type( ngridq, bvec, vql_full )
      call mfi%set_localization_centers( &
        left_centers=reshape( [((atposc(:, i, j), i=1, natoms(j)), j=1, nspecies)], [3, natmtot] ), &
        right_centers=reshape( [((atposc(:, i, j), i=1, natoms(j)), j=1, nspecies)], [3, natmtot] ), &
        coordinates='cartesian' )

      ! set dynamical matrices of q-grid
      allocate( dynq_full(3*natmtot, 3*natmtot, mfi%np), source=zzero )
      do iq = 1, mfi%np
        call find_equivalent_wavevectors( 3, mfi%vpl(:, iq), vql, nq, symlat(:, :, lsplsymc(1:nsymcrys)), nsymcrys, iq_idx, isym_idx, first_only=.true. )
        call terminate_if_false( size( iq_idx ) > 0, '(ph_util_setup_interpolation) &
          No equivalent q-point found in list.' )
        call ph_util_symapp_dyn( isym_idx(1), vql(:, iq_idx(1)), dynq(:, :, iq_idx(1)), dynq_full(:, :, iq) )
      end do

      ! generate dynamical matrices in real space
      allocate( dynmatr(3*natmtot, 3*natmtot, mfi%nr) )
      call mfi%transform_p2R( [3*natmtot, 3*natmtot], 1, dynq_full, (3*natmtot)**2, 1, dynmatr, (3*natmtot)**2, 1 )

      ! clean up
      deallocate( dynq, dynq_full, vql_full )
      if( allocated( iq_idx ) ) deallocate( iq_idx )
      if( allocated( isym_idx ) ) deallocate( isym_idx )
    end subroutine ph_util_setup_interpolation

    !> Interpolate dynamical matrix on an arbitrary set of \({\bf q}\)-points.
    !>
    !> If both \({\bf \epsilon}^\infty\) (`dielten`) and \({\bf Z}^\ast_\kappa\) (`borncharge`) is given, 
    !> the non-analytic part is added to \(D({\bf q})\).
    !>
    !> Given a matrix Fourier interpolation object `mfi` and the dynamical 
    !> matrices in real space `dynmatr` (interatomic force constants), \(D({\bf R})\), as returned
    !> by [[ph_util_setup_interpolation(subroutine)]], this subroutine interpolates the dynamical matrix 
    !> \(D({\bf q})\) on a given set of \({\bf q}\)-points.
    !> (See also [[transform_R2p(subroutine)]].)
    subroutine ph_util_interpolate( nq, vql, mfi, dynmatr, &
        dynmatq, &
        minimal_distances, symmetrize, dielten, borncharge )
      use constants, only: zone
      use matrix_fourier_interpolation, only: mfi_type
      use mod_atoms, only: natmtot
      use mod_symmetry, only: symlat, lsplsymc, nsymcrys, find_equivalent_wavevectors
      !> number of \({\bf q}\)-points
      integer, intent(in) :: nq
      !> lattice coordinates of \({\bf q}\)-points
      real(dp), intent(in) :: vql(3, nq)
      !> matrix Fourier interpolation object
      type(mfi_type), intent(in) :: mfi
      !> dynamical matrices in real space (interatomic force constants)
      complex(dp), intent(in) :: dynmatr(3*natmtot, 3*natmtot, mfi%nr)
      !> dynamical matrices at \({\bf q}\)-points
      complex(dp), allocatable, intent(out) :: dynmatq(:, :, :)
      
      !> use minimal distance interpolation (default: `.true.` )
      logical, optional, intent(in) :: minimal_distances
      !> symmetrize dynamical matrices (default: `.true.` )
      logical, optional, intent(in) :: symmetrize
      !> high-frequency dielectric tensor \({\bf \epsilon}^\infty\)
      real(dp), optional, intent(in) :: dielten(3, 3)
      !> Born effective charges \({\bf Z}^\ast_\kappa\)
      real(dp), optional, intent(in) :: borncharge(3, 3, natmtot)

      integer :: iq
      logical :: symm, mindist
      
      integer, allocatable :: iq_idx(:), isym_idx(:)

      symm = .true.
      if( present( symmetrize) ) symm = symmetrize
      mindist = .true.
      if( present( minimal_distances ) ) mindist = minimal_distances

      allocate( dynmatq(3*natmtot, 3*natmtot, nq) )
      call mfi%transform_R2p( [3*natmtot, 3*natmtot], 1, dynmatr, (3*natmtot)**2, 1, dynmatq, (3*natmtot)**2, 1, vql, &
        minimal_distances=mindist )

      ! add non-analytic part for polar materials
      if( present( dielten ) .and. present( borncharge ) ) &
        call ph_util_dynmat_lr( mfi%bvec, nq, vql, dielten, borncharge, zone, dynmatq )

      ! symmetrize dynamical matrices
      if( symm ) then
        do iq = 1, nq
          call find_equivalent_wavevectors( 3, vql(:, iq), vql(:, iq:iq), 1, symlat(:, :, lsplsymc(1:nsymcrys)), nsymcrys, iq_idx, isym_idx )
          call ph_util_symmetrize_dyn( vql(:, iq), dynmatq(:, :, iq), size( isym_idx ), isym_idx )
        end do
      end if

      ! clean up
      if( allocated( iq_idx ) ) deallocate( iq_idx )
      if( allocated( isym_idx ) ) deallocate( isym_idx )
    end subroutine ph_util_interpolate

    !> Compute non-analytic contribution to dynamical matrix in polar materials.
    !>
    !> This implements Eqs. (71-74) of 
    !> [Gonze, Lee, *Phys. Rev. B* **55**, 10355-10368 (1997).](https://doi.org/10.1103/PhysRevB.55.10355)
    subroutine ph_util_dynmat_lr( bvec, nq, vql, dielten, borncharge, alpha, dynq )
      use constants, only: twopi, fourpi, zzero, zone
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, atposc
      !> reciprocal lattice vectors (columnwise)
      real(dp), intent(in) :: bvec(3, 3)
      !> number of \({\bf q}\)-points
      integer, intent(in) :: nq
      !> \({\bf q}\)-points in lattice coordinates
      real(dp), intent(in) :: vql(3, nq)
      !> dielectric tensor \({\bf \epsilon}^\infty\)
      real(dp), intent(in) :: dielten(3, 3)
      !> Born-effective charges \({\bf Z}^\ast_\alpha\)
      real(dp), intent(in) :: borncharge(3, 3, natmtot)
      !> prefactor \(\alpha\)
      complex(dp), intent(in) :: alpha
      !> dynamical matrices \(D({\bf q})\)
      complex(dp), intent(inout) :: dynq(3*natmtot, 3*natmtot, nq)

      real(dp), parameter :: tol = 1e-12_dp          ! tolerance for terms to include
      real(dp), parameter :: gmax = -log(tol)        ! exp(-gmax) ~ tol
      real(dp), parameter :: rmax = sqrt(-log(tol))  ! exp(-rmax^2) ~ erfc(rmax) ~ tol

      integer :: ng(3), nr(3), i, ia, ias, ja, jas
      real(dp) :: lambda
      real(dp) :: avec(3, 3), dielten_inv(3, 3), det_dielten_inv, m(3, 3)

      complex(dp), allocatable :: dyn0(:,:,:), tmp(:,:)

      real(8), external :: r3mdet

      call r3minv( bvec, avec )
      avec = twopi * transpose( avec )
      call r3minv( dielten, dielten_inv )
      det_dielten_inv = r3mdet( dielten_inv )

      ! get optimal lambda
      lambda = maxval( [(dot_product( bvec(:, i), matmul( dielten, bvec(:, i) ) ), i=1, 3)] )
      lambda = sqrt( lambda * rmax / (sqrt( gmax ) * fourpi ) )

      ! get limits for recpirocal and real space sums
      m = matmul( transpose( bvec ), matmul( dielten, bvec ) )
      ng = [(ceiling( sqrt( gmax * 4.0_dp * lambda**2 / m(i, i) ) ), i=1, 3 )]
      m = matmul( transpose( avec ), matmul( dielten_inv, avec ) )
      nr = [(ceiling( rmax / sqrt( lambda**2 * m(i, i) ) ), i=1, 3)]

      ! q=0 term
      allocate( dyn0(3, 3, natmtot), source=zzero )
      allocate( tmp(3*natmtot, 3*natmtot), source=zzero )
      ! reciprocal space sum
      call reciprocal_space_sum( [0.0_dp, 0.0_dp, 0.0_dp], ng, tmp )
      ! real space sum
      call real_space_sum( [0.0_dp, 0.0_dp, 0.0_dp], nr, tmp )
      ! limiting contribution
      call limiting_contribution( tmp )      
      ! sum over atoms
      do ia = 1, natmtot
        ias = (ia - 1) * 3 + 1
        do ja = 1, natmtot
          jas = (ja - 1) * 3 + 1
          dyn0(:, :, ia) = dyn0(:, :, ia) + tmp(ias:ias+2, jas:jas+2)
        end do
      end do
      deallocate( tmp )

      ! loop over q
      !$omp parallel default( shared ) private( ia, ias )
      !$omp do
      do i = 1, nq
        ! reciprocal space sum
        call reciprocal_space_sum( vql(:, i), ng, dynq(:, :, i) )
        ! real space sum
        call real_space_sum( vql(:, i), nr, dynq(:, :, i) )
        ! limiting contribution
        call limiting_contribution( dynq(:, :, i) )      
        ! atom-diagonal part
        do ia = 1, natmtot
          ias = (ia - 1) * 3 + 1
          dynq(ias:ias+2, ias:ias+2, i) = dynq(ias:ias+2, ias:ias+2, i) - dyn0(:, :, ia)
        end do
      end do
      !$omp end do
      !$omp end parallel

      contains

        subroutine reciprocal_space_sum( vql, ng, dynq )
          real(dp), intent(in) :: vql(3)
          integer, intent(in) :: ng(3)
          complex(dp), intent(inout) :: dynq(3*natmtot, 3*natmtot)

          integer :: g1, g2, g3, is, ia, ias
          real(dp) :: t1, vt(3*natmtot), vgqc(3), gqc, gqegq, dotp
          complex(dp) :: z1, z2, vz(3*natmtot)

          z1 = alpha * fourpi / abs( r3mdet( avec ) )
          do g3 = -ng(3), ng(3)
            do g2 = -ng(2), ng(2)
              do g1 = -ng(1), ng(1)
                vgqc = matmul( bvec, dble( [g1, g2, g3] ) + vql )
                gqc = norm2( vgqc )
                if( gqc < tol ) cycle
                vgqc = vgqc / gqc
                gqegq = dot_product( vgqc, matmul( dielten, vgqc ) )
                t1 = gqc**2 * gqegq / (4.0_dp * lambda**2)
                if( (gqegq < tol) .or. (t1 > gmax) ) cycle
                z2 = z1 * exp( -t1 ) / gqegq
                call dgemv( 't', 3, 3*natmtot, 1.0_dp, borncharge, 3, vgqc, 1, 0.0_dp, vt, 1 )
                do is = 1, nspecies
                  do ia = 1, natoms(is)
                    ias = (idxas(ia, is) - 1) * 3 + 1
                    dotp = gqc * dot_product( vgqc, atposc(:, ia, is) )
                    vz(ias:ias+2) = vt(ias:ias+2) * cmplx( cos( dotp ), sin( dotp ), dp )
                  end do
                end do
                call zgemm( 'n', 'c', 3*natmtot, 3*natmtot, 1, z2, vz, 3*natmtot, vz, 3*natmtot, zone, dynq, 3*natmtot )
              end do
            end do
          end do
        end subroutine reciprocal_space_sum

        subroutine real_space_sum( vql, nr, dynq )
          real(dp), intent(in) :: vql(3)
          integer, intent(in) :: nr(3)
          complex(dp), intent(inout) :: dynq(3*natmtot, 3*natmtot)

          integer :: r1, r2, r3, is, ia, ias, js, ja, jas
          real(dp) :: vrc(3), vd(3), vdelta(3), d, H(3, 3), dotp
          complex(dp) :: z1, z2

          z1 = - alpha * lambda**3 * sqrt( det_dielten_inv )
          do r3 = -nr(3), nr(3)
            do r2 = -nr(2), nr(2)
              do r1 = -nr(1), nr(1)
                vrc = matmul( avec, dble( [r1, r2, r3] ) )
                dotp = twopi * dot_product( dble( [r1, r2, r3] ), vql )
                z2 = z1 * cmplx( cos( dotp ), sin( dotp ), dp )

                do is = 1, nspecies
                  do ia = 1, natoms(is)
                    ias = (idxas(ia, is) - 1) * 3 + 1
                    do js = 1, nspecies
                      do ja = 1, natoms(js)
                        jas = (idxas(ja, js) - 1) * 3 + 1
                        vd = vrc + atposc(:, ja, js) - atposc(:, ia, is)
                        vdelta = matmul( dielten_inv, vd )
                        d = sqrt( dot_product( vdelta, vd ) )
                        if( (d < tol) .or. (lambda * d > rmax) ) cycle
                        H = H_fun( lambda*vdelta, lambda*d )
                        H = matmul( transpose( borncharge(:, :, idxas(ia, is)) ), matmul( H, borncharge(:, :, idxas(ja, js)) ) )
                        dynq(ias:ias+2, jas:jas+2) = dynq(ias:ias+2, jas:jas+2) + z2 * H
                      end do
                    end do
                  end do
                end do

              end do
            end do
          end do
        end subroutine real_space_sum

        subroutine limiting_contribution( dynq )
          use constants, only: pi
          complex(dp), intent(inout) :: dynq(3*natmtot, 3*natmtot)

          integer :: ia, ias
          real(dp) :: m(3, 3)
          complex(dp) :: z1

          z1 = - alpha * lambda**3 * 4.0_dp / (3.0_dp * sqrt( pi )) * sqrt( det_dielten_inv )

          do ia = 1, natmtot
            ias = (ia - 1) * 3 + 1
            m = matmul( transpose( borncharge(:, :, ia) ), matmul( dielten_inv, borncharge(:, :, ia) ) )
            dynq(ias:ias+2, ias:ias+2) = dynq(ias:ias+2, ias:ias+2) + z1 * m
          end do
        end subroutine limiting_contribution

        function H_fun( x, y ) result(H)
          use constants, only: pi
          real(dp), intent(in) :: x(3), y
          real(dp) :: H(3, 3)

          real(dp) :: t1, t2, f1, f2

          t1 = erfc( y ) / y**3
          t2 = 2.0_dp / sqrt( pi ) * exp( -y**2 ) / y**2
          f1 = (3.0_dp * t1 + t2 * (3.0_dp + 2.0_dp * y**2)) / y**2
          f2 = - (t1 + t2)
          H = f2 * dielten_inv
          call dgemm( 'n', 't', 3, 3, 1, f1, x, 3, x, 3, 1.0_dp, H, 3 )
        end function H_fun
    end subroutine ph_util_dynmat_lr

end module phonons_util
