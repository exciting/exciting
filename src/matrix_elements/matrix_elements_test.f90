!> title:   matrix elements unit tests
!> summary: Module that implements testing routines
!>          for the [[matrix_elements(module)]] module
!>          and demonstrates its usage
!> author:  Sebastian Tillack
!> date:    August 2022
!> licence: GPL
!>
!> This module provides unit tests for the general matrix elements module [[matrix_elements(module)]].
!>
!> @note The numerical parameters of this tests are set to the minimum to achieve the correct
!> results with the module tolerance `tol`. For higher accuracy, the numerical parameters have
!> to be changed. Especially, `lmaxb` and `lmaxo` (and to some extend also `num_rad`) govern both
!> the accuracy and the cost of the calculations while `rgkmax` has very little impact. @endnote
!> 
!> The idea of the tests is the following
!> 
!> 1. Set up a test system for unittestability. The test system consists of a single atom
!>    at position \({\bf \tau}_\alpha\) with muffin-tin radius \(R_\alpha\) in a unit cell
!>    with lattive vectors \({\bf A} = ({\bf a}_1, {\bf a}_2, {\bf a}_3)\). This is done by the
!>    routine [[init_system(subroutine)]].
!> 2. Compute some auxiliary wavefunctions \(\psi_{n{\bf p}}({\bf r})\) by
!>    expanding plane waves \({\rm e}^{{\rm i}({\bf p+G})\cdot{\bf r}}\) into 
!>    (L)APW and LO basis functions \(\phi_{\nu{\bf p}}({\bf r})\)
!>    \[ \psi_{n{\bf p}}({\bf r}) = \sum\limits_{\nu} C_{\nu n}({\bf p})\, \phi_{\nu{\bf p}}({\bf r})
!>       \approx \frac{1}{\sqrt{\Omega}}\, {\rm e}^{{\rm i}({\bf p+G}_n)\cdot{\bf r}} \]
!>    for a set of test reciprocal lattice vectors \({\bf G}_n\). This is done by the routines
!>    [[gen_obb(subroutine)]], [[gen_obp(subroutine)]], and [[gen_evec(subroutine)]].
!> 3. Compute various matrix elements using the auxiliary wavefunctions and compare them
!>    against their analytic results:
!>    1. Check the orthogonality of the wavefunctions
!>       using [[me_test_orthogonality(subroutine)]].
!>    2. Check the muffin-tin volume integral of the overlap of two wavefunctions
!>       using [[me_test_mt(subroutine)]].
!>    3. Check the muffin-tin surface integral of the overlap of two wavefunctions 
!>       using [[me_test_sf(subroutine)]].
!>    4. Check the muffin-tin volume integral of the overlap of two wavefunction gradients
!>       using [[me_test_grad(subroutine)]].
!>    5. Check the muffin-tin surface integral of the overlap of two wavefunction gradients
!>       using [[me_test_sfgrad(subroutine)]].
!>    6. Check the orthogonality of the matrix elements of different plane-wave operators
!>       and two wavefunctions
!>       using [[me_test_operator(subroutine)]].
!>    7. Check Gauss' theorem by comparing the muffin-tin volume integral of the gradient
!>       of different plane-wave operators and two wavefunctions against the respective
!>       muffin-tin surface integral
!>       using [[me_test_gauss(subroutine)]].
module matrix_elements_test
  use matrix_elements
  use unit_test_framework, only : unit_test_type
  use precision, only: dp
  use constants
  use mod_kpointset

  use muffin_tin_basis, only: mt_basis_type
  use mod_APW_LO, only: apwordmax, apword, nlomax, nlorb, nlotot, lorbl, lolmmax
  use mod_atoms, only: nspecies, natoms, natmtot, idxas, spr, atposc
  use mod_eigensystem, only: idxlo
  use mod_muffin_tin, only: nrmt, rmt, nrmtmax, lmmaxapw, idxlm
  use mod_lattice, only: omega, avec, ainv, bvec
  use modinput

  implicit none
  private

  !> tolerance
  real(dp), parameter :: tol = 1.e-6

  ! test system
  !> lattice vectors
  real(dp), parameter :: latvec(3,3) = reshape( [1.23_dp, 0.04_dp, 0.15_dp, &
                                                 0.09_dp, 1.18_dp, 0.11_dp, &
                                                 0.17_dp, 0.08_dp, 1.09_dp], [3,3])
  !> atom position in lattice coordinates
  real(dp), parameter :: atposl(3) = [0.18_dp, 0.26_dp, 0.37_dp]
  !> muffin-tin radius
  real(dp), parameter :: mt_rad = 0.362_dp
  !> number of radial grid points (inside muffin-tin)
  integer, parameter :: num_rad = 150
  !> number of extra radial grid points (outside muffin-tin)
  integer, parameter :: num_rad_add = 20
  !> maximum angular momentum for basis functions
  integer, parameter :: lmaxb = 10
  !> maximum angular momentum for operator
  integer, parameter :: lmaxo = 8
  !> `rgkmax` determines the cut-off of the (L)APWs
  real(dp), parameter :: rgkmax = 6._dp
  !> wavevector in lattice coordinates
  real(dp), parameter :: vpl(3) = [0.24_dp, 0.18_dp, 0.37_dp]
  !> number of test G-vectors
  integer, parameter :: ng = 3
  !> set of test G-vectors in lattice coordinates
  integer, parameter :: test_gvl(3,ng) = reshape( [0,0,0, &
                                                   1,0,0, &
                                                   1,1,1], [3,ng])
  
  !> fake groundstate variables
  type(groundstate_type), target :: groundstate
  !> fake structure variables
  type(structure_type), target :: structure
  !> wavevector in Cartesian coordinates
  real(dp) :: vpc(3)
  !> set of test G-vectors in Cartesian coordinates
  real(dp) :: test_gvc(3,ng)
  !> basis radial functions
  type(mt_basis_type) :: basis
  !> set of k-vectors for wavefunctions
  type(k_set) :: kset
  !> set of G-vectors for Fourier series
  type(G_set) :: Gset
  !> set of G+k vectors for LAPW expansion
  type(Gk_set) :: Gkset
  !> overlap between basis functions
  complex(dp), allocatable :: obb(:,:)
  !> overlap between basis functions and plane wave
  complex(dp), allocatable :: obp(:,:)
  !> radial integrals times Gaunt coefficients for overlap volume integrals
  complex(dp), allocatable :: rigntov(:,:,:)
  !> radial integrals times Gaunt coefficients for overlap surface integrals
  complex(dp), allocatable :: rigntos(:,:,:,:)
  !> radial integrals times Gaunt coefficients for overlap gradient integrals
  complex(dp), allocatable :: rigntog(:,:,:,:)
  complex(dp), allocatable :: rigntolg(:,:,:,:)
  complex(dp), allocatable :: rigntorg(:,:,:,:)
  complex(dp), allocatable :: rigntomg(:,:,:,:,:)
  !> radial integrals times Gaunt coefficients for mix of overlap surface and gradient integrals
  complex(dp), allocatable :: rigntolgs(:,:,:,:,:)
  complex(dp), allocatable :: rigntorgs(:,:,:,:,:)
  !> characteristic function in reciprocal space
  complex(dp), allocatable :: cfunig(:)
  !> radial integrals times Gaunt coefficients for plane-wave operator volume integrals
  complex(dp), allocatable :: rigntpv(:,:,:,:)
  !> radial integrals times Gaunt coefficients for plane-wave operator surface integrals
  complex(dp), allocatable :: rigntps(:,:,:,:)
  !> radial integrals times Gaunt coefficients for plane-wave operator gradient integrals
  complex(dp), allocatable :: rigntpg(:,:,:,:)
  !> plane-wave operator times characteristic function in reciprocal space
  complex(dp), allocatable :: opig(:,:)
  !> expansion coefficients for plane waves in basis functions
  complex(dp), allocatable :: evec(:,:)

  public :: ng
  public :: me_test_init, me_test_orthogonality, me_test_mt, me_test_sf, me_test_grad, me_test_sfgrad, &
            me_test_operator, me_test_gauss

  contains

    !> Initialize the test system and precompute different quantities for the actual tests.
    subroutine me_test_init
      ! initialize test system
      call init_system
      ! initialize matrix elements module
      call me_init( basis, lmaxo, Gset )
      ! prepare overlap integrals
      call prepare_overlap
      ! prepare plane-wave operator integrals
      call prepare_operator
      ! generate overlap between basis functions
      call gen_obb
      ! generate overlap between basis functions and plane waves
      call gen_obp
      ! find coefficients for expansion of plane waves in basis functions
      call gen_evec
    end subroutine me_test_init

    !> Initialize test module for unit-testability
    !>
    !> This includes the initialization of the test system (single atom
    !> in the unit cell), the setting of necessary input parameters and 
    !> globals that other routines rely on and the creation of suitable
    !> (L)APW and LO radial basis functions.
    subroutine init_system
      use grid_utils, only: linspace
      use mod_symmetry, only: nsymcrys, nsymlat, lsplsymc
      use mod_spin, only: nspnfv
      
      integer :: l, m, lm, ni
      real(dp) :: ri

      real(dp), allocatable :: rgrid(:)

      real(dp), external :: r3mdet

      ! necessary input parameters
      groundstate%lmaxapw = lmaxb
      groundstate%stypenumber = 1
      input%groundstate => groundstate
      structure%epslat = 1e-6_dp
      input%structure => structure
      ! missing constants
      if( allocated( zil ) ) deallocate( zil )
      allocate( zil(0:lmaxb), source=zone )
      do l = 1, lmaxb
        zil(l) = cmplx( -aimag( zil(l-1) ), dble( zil(l-1) ), dp )
      end do
      ! lattice
      avec = latvec
      omega = r3mdet( avec )
      call r3minv( avec, ainv )
      bvec = twopi * transpose( ainv )
      ! single atom
      nspecies = 1
      natmtot = 1
      if( allocated( natoms ) ) deallocate( natoms )
      allocate( natoms(1), source=1 )
      idxas = 0
      idxas(1, 1) = 1
      atposc = 0.0_dp
      atposc(:, 1, 1) = matmul( avec, atposl )
      rmt = 0.0_dp
      rmt(1) = mt_rad
      ! other globals
      nspnfv = 1
      if( allocated( idxlm ) ) deallocate( idxlm )
      allocate( idxlm(0:lmaxb, -lmaxb:lmaxb), source=0 )
      lm = 0
      do l = 0, lmaxb
        do m = -l, l
          lm = lm + 1
          idxlm(l, m) = lm
        end do
      end do
      ! symmetries
      nsymcrys = 1
      nsymlat = 1
      lsplsymc = 0
      lsplsymc(1) = 1
      ! initialize radial grid
      ri = 0.1_dp ! inner part of muffin-tin used for cubic sampling
      if( allocated( spr ) ) deallocate( spr )
      allocate( spr(num_rad+num_rad_add, 1) )
      nrmt = 0
      nrmt(1) = num_rad
      nrmtmax = nrmt(1)
      rgrid = linspace( 1.e-6_dp, 1.0_dp, num_rad_add)
      spr(1:num_rad_add, 1) = ri * rgrid**3
      rgrid = linspace( ri, (1.0_dp-ri*num_rad_add/num_rad)*num_rad/(num_rad-num_rad_add), num_rad+1 )
      spr(num_rad_add+1:, 1) = rgrid(2:)
      spr = spr * mt_rad
      ! set up missing module variables
      vpc = matmul( bvec, vpl )
      test_gvc = matmul( bvec, dble( test_gvl ) )
      call generate_k_vectors( kset, bvec, [1,1,1], vpl, .false., .false. )
      call generate_G_vectors( Gset, bvec, [[0,0,0],[0,0,0]], rgkmax/mt_rad*2, auto_intgv=.true. )
      call generate_Gk_vectors( Gkset, kset, Gset, rgkmax/mt_rad )
      ! set up basis
      call init_basis
    end subroutine init_system

    !> Initialize the (L)APW and LO radial basis functions.
    !> 
    !> Our goal is to expand plane waves in (L)APW and LO basis functions.
    !> This is very difficult using ordinary, atomic like radial functions
    !> in the muffin tins. Instead, inspired by the spherical harmonics 
    !> expansion of a plane, we use spherical Bessel functions \(j_l(x)\) as
    !> radial functions.
    !>
    !> For the (L)APW radial functions, we use 
    !> \[ u^{\alpha}_{l,\xi}(r) \propto \begin{cases}
    !>  j_l(|{\bf p}|\, r) & \xi = 0 \\
    !>  j_l(|{\bf p + G}_1|\, r) & \xi = 1 \end{cases} \;, \]
    !> where \({\b p}\) is the wavevector of the test wavefunction and \(\xi\) is 
    !> the (L)APW order (energy derivative). We use \(\xi = 0,1\).
    !>
    !> For the LO radial functions, we use 
    !> \[f^{\alpha}_l(r) \propto j_l(|{\bf p + G}_n|r) \;. \] 
    !> We add LOs for all vectors \({\bf G}_n\) for which we want to 
    !> expand plane waves into basis functions. Note, that these LOs do not 
    !> go to zero on the muffin-tin boundary.
    subroutine init_basis
      use mod_APW_LO, only: apwfr, lofr, MTBasisInit

      integer :: nrtot, l, ir, ilo, ig
      real(dp) :: vgpc(3)

      real(dp), allocatable :: f(:,:,:), g(:), cf(:,:), jbessel(:)

      nrtot = size( spr, dim=1 )
      allocate( f(nrtot, 0:lmaxb, 3), g(nrtot), cf(3, nrtot), jbessel(0:lmaxb) )

      ! ** set up LAPWs
      apwordmax = 2 ! use LAPWs
      apword = 0
      apword(0:lmaxb, 1) = apwordmax
      lmmaxapw = (lmaxb + 1)**2
      if( allocated( apwfr ) ) deallocate( apwfr )
      allocate( apwfr(nrmtmax, 2, apwordmax, 0:lmaxb, 1), source=0.0_dp )
      do ir = 1, nrtot
        call sbessel( lmaxb, norm2( vpc ) * spr(ir, 1), jbessel )
        f(ir, :, 1) = jbessel
        call sbessel( lmaxb, norm2( vpc + test_gvc(:, 2) ) * spr(ir, 1), jbessel )
        f(ir, :, 2) = jbessel
      end do
      ! normalize and compute radial derivatives
      do l = 0, lmaxb
        f(:, l, 3) = f(:, l, 1)**2 * spr(:, 1)**2
        call fderiv( -1, nrtot, spr(:, 1), f(:, l, 3), g, cf )
        f(:, l, 1) = f(:, l, 1) / sqrt( g(nrmt(1)) )
        call fderiv( 1, nrtot, spr(:, 1), f(:, l, 1), g, cf )
        apwfr(:, 1, 1, l, 1) = f(1:nrmt(1), l, 1)
        apwfr(:, 2, 1, l, 1) = g(1:nrmt(1))

        f(:, l, 3) = f(:, l, 2)**2 * spr(:, 1)**2
        call fderiv( -1, nrtot, spr(:, 1), f(:, l, 3), g, cf )
        f(:, l, 2) = f(:, l, 2) / sqrt( g(nrmt(1)) )
        call fderiv( 1, nrtot, spr(:, 1), f(:, l, 2), g, cf )
        apwfr(:, 1, 2, l, 1) = f(1:nrmt(1), l, 2)
        apwfr(:, 2, 2, l, 1) = g(1:nrmt(1))
      end do

      ! ** set up LOs
      nlotot = 0
      nlomax = (lmaxb + 1) * ng
      lolmmax = (lmaxb + 1)**2
      if( allocated( lofr ) ) deallocate( lofr )
      allocate( lofr(nrmtmax, 2, nlomax, 1), source=0.0_dp )
      nlorb = 0
      nlorb(1) = nlomax
      lorbl = 0
      ilo = 0
      do ig = 1, ng
        vgpc = vpc + test_gvc(:, ig)
        do ir = 1, nrtot
          call sbessel( lmaxb, norm2( vgpc ) * spr(ir, 1), jbessel )
          f(ir, :, 1) = jbessel
        end do
        do l = 0, lmaxb
          ilo = ilo + 1
          lorbl(ilo, 1) = l
          nlotot = nlotot + 2 * l + 1
          ! normalize and compute radial derivatives
          f(:, l, 3) = f(:, l, 1)**2 * spr(:, 1)**2
          call fderiv( -1, nrtot, spr(:, 1), f(:, l, 3), g, cf )
          f(:, l, 1) = f(:, l, 1) / sqrt( g(nrmt(1)) )
          call fderiv( 1, nrtot, spr(:, 1), f(:, l, 1), g, cf )
          lofr(:, 1, ilo, 1) = f(1:nrmt(1), l, 1)
          lofr(:, 2, ilo, 1) = g(1:nrmt(1))
        end do
      end do
      call genidxlo

      deallocate( f, g, cf, jbessel )

      basis = mt_basis_type( spr(:, 1:1), nrmt(1:1), apwfr, lofr, lmaxb, &
        apword(:, 1:1), nlorb(1:1), lorbl(:, 1:1) )
    end subroutine init_basis

    !> Prepare the calculation of overlap matrix elements.
    !>
    !> Use the routine [[me_mt_prepare(subroutine)]] to precompute the radial integrals
    !> times Gaunt coefficients for the overlap of basis functions / wavefuncions, i.e.,
    !> for the identity operator. We prepare the integrals for muffin-tin volume integrals
    !> of both basis functions / wavefunctions 
    !> \[ \langle \phi_{\mu{\bf p}} | \phi_{\nu{\bf p}} \rangle_{\alpha} \]
    !> and their gradients 
    !> \[ \langle \nabla \phi_{\mu{\bf p}} | \phi_{\nu{\bf p}} \rangle_{\alpha}
    !>       + \langle \phi_{\mu{\bf p}} | \nabla \phi_{\nu{\bf p}} \rangle_{\alpha} \]
    !> as well as for muffin-tin surface integrals of basis functions / wavefunctions
    !> \[ \langle \phi_{\mu{\bf p}} | \phi_{\nu{\bf p}} \rangle_{\partial\alpha} \;. \]
    !>
    !> Use the routine [[me_ir_prepare(subroutine)]] to precompute the reciprocal space
    !> representation of the product of the identity operator and the characteristic
    !> function. 
    subroutine prepare_overlap
      integer :: ip, ipp

      real(dp), allocatable :: opir(:), opmt(:,:,:)
      complex(dp), allocatable :: tmp(:,:,:,:)

      allocate( opmt(1, nrmtmax, 1), source=1.0_dp/y00 )
      allocate( opir(Gset%ngrtot), source=1.0_dp )
      call me_mt_alloc( rigntov )
      call me_mt_alloc( rigntos, 3 )
      call me_mt_alloc( rigntog, 3 )
      call me_mt_alloc( rigntolg, 3 )
      call me_mt_alloc( rigntorg, 3 )
      call me_mt_alloc( tmp, 9 )
      rigntomg = reshape(tmp, [size(tmp, dim=1), size(tmp, dim=2), size(tmp, dim=3), 3, 3])
      rigntolgs = reshape(tmp, [size(tmp, dim=1), size(tmp, dim=2), size(tmp, dim=3), 3, 3])
      rigntorgs = reshape(tmp, [size(tmp, dim=1), size(tmp, dim=2), size(tmp, dim=3), 3, 3])
      call me_ir_alloc( cfunig )

      call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntov(:, :, 1) )
      do ip = 1, 3
        call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntos(:, :, 1, ip), surface_integral=ip )
        call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntog(:, :, 1, ip), left_gradient=ip )
        call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zone, rigntog(:, :, 1, ip), right_gradient=ip )
        call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntolg(:, :, 1, ip), left_gradient=ip )
        call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntorg(:, :, 1, ip), right_gradient=ip )
        do ipp = 1, 3
          call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntomg(:, :, 1, ip, ipp), left_gradient=ip, right_gradient=ipp )
          call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntolgs(:, :, 1, ip, ipp), left_gradient=ip, surface_integral=ipp )
          call me_mt_prepare( 1, 1, 0, zone, opmt(:, :, 1), zzero, rigntorgs(:, :, 1, ip, ipp), right_gradient=ip, surface_integral=ipp )
        end do
      end do
      call me_ir_prepare( zone, opir, zzero, cfunig )

      deallocate( opir, opmt, tmp )
    end subroutine prepare_overlap

    !> Prepare the calculation of plane-wave operator matrix elements.
    !>
    !> Use the routine [[me_mt_prepare(subroutine)]] to precompute the radial integrals
    !> times Gaunt coefficients for the matrix elements of different plane-wave operators 
    !> \({\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}}\). We prepare the integrals for muffin-tin 
    !> volume integrals of both basis functions / wavefunctions and operators
    !> \[ \langle \phi_{\mu{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \phi_{\nu{\bf p}} \rangle_{\alpha} \]
    !> and their gradients
    !> \[ \langle \nabla \phi_{\mu{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \phi_{\nu{\bf p}} \rangle_{\alpha}
    !>    + \langle \phi_{\mu{\bf p}} | \nabla {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \phi_{\nu{\bf p}} \rangle_{\alpha}
    !>    + \langle \phi_{\mu{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \nabla \phi_{\nu{\bf p}} \rangle_{\alpha} \]
    !> as well as for
    !> muffin-tin surface integrals of basis functions / wavefunctions and operators
    !> \[ \langle \phi_{\mu{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \phi_{\nu{\bf p}} \rangle_{\partial\alpha} \;. \]
    !>
    !> Use the routine [[me_ir_prepare(subroutine)]] to precompute the reciprocal space
    !> representation of the product of the plane-wave operators and the characteristic
    !> function. 
    subroutine prepare_operator
      use math_utils, only: plane_wave_in_spherical_harmonics
      integer :: ip, ig, g(3)
      real(dp) :: dotp

      complex(dp), allocatable :: opir(:), opmt(:,:,:,:)
      complex(dp), allocatable :: zfmt(:,:)

      allocate( opmt((lmaxo+1)**2, nrmtmax, 1, 0:3) )
      allocate( opir(Gset%ngrtot) )
      call me_mt_alloc( rigntpv, ng )
      call me_mt_alloc( rigntps, 3*ng )
      call me_mt_alloc( rigntpg, 3*ng )
      call me_ir_alloc( opig, ng )

      do ig = 1, ng
        g = test_gvl(:, ig)
        dotp = twopi * dot_product( dble(g), atposl )
        ! generate operator in muffin-tin sphere
        call plane_wave_in_spherical_harmonics( test_gvc(:, ig), spr(1:nrmtmax, 1), lmaxo, zfmt )
        opmt(:, :, 1, 0) = zfmt * cmplx( cos(dotp), sin(dotp), dp )
        call gradzfmt( lmaxo, nrmt(1), spr(:, 1), (lmaxo+1)**2, nrmtmax, opmt(:, :, 1, 0), opmt(:, :, 1, 1:3) )
        call me_mt_prepare( 1, 1, lmaxo, zone, opmt(:, :, 1, 0), zzero, rigntpv(:, :, 1, ig) )
        do ip = 1, 3
          call me_mt_prepare( 1, 1, lmaxo, zone, opmt(:, :, 1, 0), zzero, rigntps(:, :, 1, (ig-1)*3+ip), surface_integral=ip )
          call me_mt_prepare( 1, 1, lmaxo, zone, opmt(:, :, 1, 0), zzero, rigntpg(:, :, 1, (ig-1)*3+ip), left_gradient=ip )
          call me_mt_prepare( 1, 1, lmaxo, zone, opmt(:, :, 1, 0), zone, rigntpg(:, :, 1, (ig-1)*3+ip), right_gradient=ip )
          call me_mt_prepare( 1, 1, lmaxo, zone, opmt(:, :, 1, ip), zone, rigntpg(:, :, 1, (ig-1)*3+ip) )
        end do          
        ! generator operator in interstitial region
        opir = zzero
        opir( Gset%igfft( Gset%ivgig(g(1), g(2), g(3)) ) ) = zone
        call zfftifc( 3, Gset%ngrid, 1, opir )
        call me_ir_prepare( zone, opir, zzero, opig(:, ig) )
      end do

      deallocate( opir, opmt, zfmt )
    end subroutine prepare_operator

    !> Generate the overlap between (L)APW+LO basis functions.
    !>
    !> Compute the overlaps
    !> \[ S_{\mu \nu}({\bf q}) = \int\limits_{\Omega} \phi_{\mu{\bf p}}^\ast ({\bf r}) \,
    !>       \phi_{\nu{\bf p}}({\bf r})\, {\rm d}r^3 \;, \]
    !> where \(\phi_{\nu{\bf p}}({\bf r})\) are (L)APW or LO basis function.
    subroutine gen_obb
      use mod_Gkvector, only: ngkmax_ptr

      integer :: nmat
      integer, target :: ngkmax

      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax
      nmat = Gkset%ngk(1, 1) + nlotot
      
      allocate( obb(nmat, nmat) )
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) ) 
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), zone, rigntov(:, :, 1), zzero, obb )
      call me_ir_mat( Gkset, 1, zone, cfunig, zone, obb )

      deallocate( apwalm)
    end subroutine gen_obb

    !> Generate the overlap between (L)APW+LO basis functions and different plane waves.
    !>
    !> Compute the overlaps
    !> \[ T_{\nu n}({\bf p}) = \frac{1}{\sqrt{\Omega}} \int\limits_{\Omega} \phi_{\nu{\bf p}}^\ast ({\bf r}) \, 
    !>       {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}}\, {\rm d}r^3 \;, \]
    !> where \(\phi_{\nu{\bf p}}({\bf r})\) are (L)APW or LO basis function.
    subroutine gen_obp
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: plane_wave_in_spherical_harmonics

      integer :: nmat, ig, igk
      integer :: l, m, lm, lam
      integer, target :: ngkmax
      real(dp) :: vgpc(3), t1

      real(dp), allocatable :: gc(:), vgc(:,:), f1(:), f2(:), g(:,:), cf(:,:)
      complex(dp), allocatable :: pwrfun(:,:), rint(:), apwalm(:,:,:,:), cmatch(:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax
      nmat = Gkset%ngk(1, 1) + nlotot

      allocate( obp(nmat, ng), source=zzero )
      allocate( gc(ngkmax), vgc(3, ngkmax) )
      allocate( g(nrmtmax, 2), cf(3, nrmtmax) )
      allocate( rint(basis%n_basis_fun_max) )
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) ) 
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      call basis%get_basis_transform( 1, Gkset%ngk(1, 1), nlotot, idxlo(:, :, 1), apwalm(:, :, :, 1), cmatch )
      do ig = 1, ng
        ! interstitial region
        do igk = 1, Gkset%ngk(1, 1)
          vgc(:, igk) = Gset%vgc(:, Gkset%igkig(igk, 1, 1)) - test_gvc(:, ig)
          gc(igk) = norm2( vgc(:, igk) )
        end do
        call gencfunig( Gkset%ngk(1, 1), gc, vgc, obp(:, ig) )
        ! muffin-tin
        vgpc = vpc + test_gvc(:, ig)
        t1 = dot_product( vgpc, atposc(:, 1, 1) )
        call plane_wave_in_spherical_harmonics( vgpc, basis%rad_grid(1:basis%n_rad_grid(1), 1), lmaxb, pwrfun )
        pwrfun = pwrfun * cmplx( cos(t1), sin(t1), dp ) / sqrt( omega )
        lm = 0
        do l = 0, lmaxb
          do m = -l, l
            lm = lm + 1
            do lam = 1, basis%n_rad_fun(l, 1)
              f1 = basis%get_rad_fun( l, 1, 1, lam ) *  dble( pwrfun(lm, :) ) * basis%rad_grid(:, 1)**2
              f2 = basis%get_rad_fun( l, 1, 1, lam ) * aimag( pwrfun(lm, :) ) * basis%rad_grid(:, 1)**2
              call fderiv( -1, nrmt(1), basis%rad_grid(:, 1), f1, g(:, 1), cf )
              call fderiv( -1, nrmt(1), basis%rad_grid(:, 1), f2, g(:, 2), cf )
              rint(basis%idx_basis_fun(lm, lam, 1)) = cmplx( g(basis%n_rad_grid(1), 1), g(basis%n_rad_grid(1), 2), dp )
            end do
          end do
        end do
        call zgemv( 'c', basis%n_basis_fun(1), nmat, zone, &
               cmatch, basis%n_basis_fun_max, &
               rint, 1, zone, &
               obp(:, ig), 1 )
      end do

      deallocate( gc, vgc, f1, f2, g, cf, rint, apwalm, cmatch )
    end subroutine gen_obp

    !> Generate the wavefunction expansion coefficients.
    !>
    !> Compute the expansion coefficients \(C_{\nu n}({\bf p})\) in the expansion
    !> \[ \psi_{n{\bf p}}({\bf r}) = \sum\limits_{\nu} C_{\nu n}({\bf p})\, \phi_{\nu{\bf p}}({\bf r})
    !>    \approx \frac{1}{\sqrt{\Omega}}\, {\rm e}^{{\rm i}({\bf p+G}_n)\cdot{\bf r}} \]
    !> by solving the linear system
    !> \[ {\bf S}({\bf p}) \cdot {\bf C}({\bf p}) = {\bf T}({\bf p}) \;. \]
    !> See also [[gen_obb(subroutine)]] and [[gen_obp(subroutine)]].
    subroutine gen_evec
      integer :: nmat, info

      integer, allocatable :: ipiv(:)
      complex(dp), allocatable :: tmp(:,:)

      nmat = Gkset%ngk(1, 1) + nlotot

      allocate( evec, source=obp )
      allocate( tmp, source=obb )
      allocate( ipiv(nmat) )
      call zgesv( nmat, ng, tmp, nmat, ipiv, evec, nmat, info )
      deallocate( tmp, ipiv )
    end subroutine gen_evec

    !> Check the orthogonality of the wavefunctions.
    !>
    !> Compute the integrals
    !> \[ \begin{align*}
    !>       \langle \psi_{m{\bf p}} | \psi_{n{\bf p}} \rangle
    !>       &\overset{?}{=} \frac{1}{\Omega} \int\limits_{\Omega}
    !>          {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf r}}\, {\rm d}r^3 \\
    !>       &= \delta_{mn} \;.
    !>    \end{align*} \]
    subroutine me_test_orthogonality( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ig
      integer, target :: ngkmax

      complex(dp), allocatable :: mat_ref(:,:), mat_res(:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_ref(ng, ng), source=zzero )
      allocate( mat_res(ng, ng) )

      ! set up reference
      do ig = 1, ng
        mat_ref(ig, ig) = zone
      end do

      ! compute matrix elements using module
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntov(:, :, 1), zzero, mat_res )
      call me_ir_mat( Gkset, 1, evec, zone, cfunig, zone, mat_res )

      ! assertion
      call test_report%assert( all_close( mat_ref, mat_res, tol ), &
             'Expected: Wave functions are orthonormal.')

      deallocate( apwalm, mat_ref, mat_res )
    end subroutine me_test_orthogonality

    !> Check the muffin-tin volume integral of the overlap of two wavefunctions.
    !> 
    !> Compute the integrals
    !> \[ \begin{align*}
    !>       \langle \psi_{m{\bf p}} | \psi_{n{\bf p}} \rangle_{\alpha}
    !>       &\overset{?}{=} \frac{1}{\Omega} \int\limits_{B_{{\bf \tau}_\alpha}(R_\alpha)}
    !>          {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf r}}\, {\rm d}r^3 \\
    !>       &= \frac{4\pi}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) \;.
    !>    \end{align*} \]
    subroutine me_test_mt( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ig1, ig2, gl(3)
      real(dp) :: gc(3), g, dotp, besselj(0:1)
      complex(dp) :: z1
      integer, target :: ngkmax

      complex(dp), allocatable :: mat_ref(:,:), mat_res(:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_ref(ng, ng) )
      allocate( mat_res(ng, ng) )

      ! set up reference
      do ig2 = 1, ng
        do ig1 = 1, ng
          gl = test_gvl(:, ig2) - test_gvl(:, ig1)
          gc = test_gvc(:, ig2) - test_gvc(:, ig1)
          g = norm2( gc )
          call sbessel( 1, g*mt_rad, besselj )
          dotp = twopi * dot_product( dble(gl), atposl )
          z1 = fourpi * cmplx( cos(dotp), sin(dotp), dp ) / omega
          if( g < 1.d-12 ) then
            mat_ref(ig1, ig2) = z1 * mt_rad**3 / 3._dp
          else
            mat_ref(ig1, ig2) = z1 * mt_rad**2 / g * besselj(1)
          end if
        end do
      end do

      ! compute matrix elements using module
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) ) 
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntov(:, :, 1), zzero, mat_res )

      ! assertion
      call test_report%assert( all_close( mat_ref, mat_res, tol ), &
             'Expected: Muffin-tin volume integrals agree with analytic result.')

      deallocate( apwalm, mat_ref, mat_res )
    end subroutine me_test_mt

    !> Check the muffin-tin surface integral of the overlap of two wavefunctions.
    !> 
    !> Compute the integrals
    !> \[ \begin{align*}
    !>       \langle \psi_{m{\bf p}} | \psi_{n{\bf p}} \rangle_{\partial\alpha}
    !>       & \overset{?}{=} \frac{1}{\Omega} \oint\limits_{S_{{\bf \tau}_\alpha}(R_\alpha)}
    !>          {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf r}}\, \hat{\bf e}\, {\rm d}S \\
    !>       &= \frac{4\pi{\rm i}}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf G}_n - {\bf G}_m) \;. 
    !>    \end{align*} \]
    subroutine me_test_sf( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ip, ig1, ig2, gl(3)
      real(dp) :: gc(3), g, dotp, besselj(0:1)
      complex(dp) :: z1
      integer, target :: ngkmax
      character(256) :: errmsg

      complex(dp), allocatable :: mat_ref(:,:,:), mat_res(:,:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_ref(ng, ng, 3) )
      allocate( mat_res(ng, ng, 3) )

      ! set up reference
      do ig2 = 1, ng
        do ig1 = 1, ng
          gl = test_gvl(:, ig2) - test_gvl(:, ig1)
          gc = test_gvc(:, ig2) - test_gvc(:, ig1)
          g = norm2( gc )
          call sbessel( 1, g*mt_rad, besselj )
          dotp = twopi * dot_product( dble(gl), atposl )
          z1 = fourpi * cmplx( -sin(dotp), cos(dotp), dp ) / omega
          if( g < 1.d-12 ) then
            mat_ref(ig1, ig2, :) = zzero
          else
            mat_ref(ig1, ig2, :) = z1 * mt_rad**2 / g * besselj(1) * gc(:)
          end if
        end do
      end do

      ! compute matrix elements using module
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) ) 
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      do ip = 1, 3
        call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntos(:, :, 1, ip), zzero, mat_res(:, :, ip) )

        ! assertion
        write( errmsg, '("Expected: Muffin-tin surface integrals agree with analytic result &
                          for component ",i1," of surface normal vector.")' ) ip
        call test_report%assert( all_close( mat_ref(:, :, ip), mat_res(:, :, ip), tol ), trim( errmsg ) )
      end do

      deallocate( apwalm, mat_ref, mat_res )
    end subroutine me_test_sf

    !> Check the muffin-tin volume integral of the overlap of two wavefunction gradients.
    !>
    !> Compute the integrals
    !> \[ \begin{align*}
    !>       \langle \nabla \psi_{m{\bf p}} | \psi_{n{\bf p}} \rangle_{\alpha}
    !>          + \langle \psi_{m{\bf p}} | \nabla \psi_{n{\bf p}} \rangle_{\alpha}
    !>       &\overset{?}{=} \frac{1}{\Omega}\, {\rm i}({\bf G}_n - {\bf G}_m) 
    !>          \int\limits_{B_{{\bf \tau}_\alpha}(R_\alpha)}
    !>          {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf r}}\, {\rm d}r^3 \\
    !>       &= \frac{4\pi{\rm i}}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf G}_n - {\bf G}_m) \; ,
    !>    \end{align*} \]
    !> as well as
    !> \[ \langle \nabla \psi_{m{\bf p}} | \psi_{n{\bf p}} \rangle_{\alpha}
    !>        = \frac{4\pi{\rm i}}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf p} + {\bf G}_n) \; , \]
    !> \[ \langle \psi_{m{\bf p}} | \nabla \psi_{n{\bf p}} \rangle_{\alpha}
    !>        = -\frac{4\pi{\rm i}}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf p} + {\bf G}_m) \; , \]
    !> and
    !> \[ \langle \nabla \psi_{m{\bf p}} | \nabla^\top \psi_{n{\bf p}} \rangle_{\alpha}
    !>        = \frac{4\pi}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf p} + {\bf G}_n) \cdot ({\bf p} + {\bf G}_m)^\top \; . \]
    subroutine me_test_grad( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ip, ipp, ig1, ig2, gl(3)
      real(dp) :: gc(3), g, dotp, besselj(0:1)
      complex(dp) :: z1
      integer, target :: ngkmax
      character(256) :: errmsg

      complex(dp), allocatable :: mat_ref(:,:,:), mat_refl(:,:,:), mat_refr(:,:,:), mat_refm(:,:,:,:)   
      complex(dp), allocatable :: mat_res(:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_ref(ng, ng, 3), mat_refl(ng, ng, 3), mat_refr(ng, ng, 3), mat_refm(ng, ng, 3, 3) )
      allocate( mat_res(ng, ng) )

      ! set up reference
      do ig2 = 1, ng
        do ig1 = 1, ng
          gl = test_gvl(:, ig2) - test_gvl(:, ig1)
          gc = test_gvc(:, ig2) - test_gvc(:, ig1)
          g = norm2( gc )
          call sbessel( 1, g*mt_rad, besselj )
          dotp = twopi * dot_product( dble(gl), atposl )
          z1 = fourpi * cmplx( -sin(dotp), cos(dotp), dp ) / omega
          if( g < 1.d-12 ) then
            mat_ref(ig1, ig2, :) = zzero
            mat_refl(ig1, ig2, :) = - z1 * mt_rad**3 / 3.0_dp * (vpc + test_gvc(:, ig1))
            mat_refr(ig1, ig2, :) =   z1 * mt_rad**3 / 3.0_dp * (vpc + test_gvc(:, ig2))
            do ipp = 1, 3
              mat_refm(ig1, ig2, :, ipp) = z1 * mt_rad**3 / 3.0_dp * (vpc + test_gvc(:, ig1)) * (vpc(ipp) + test_gvc(ipp, ig2))
            end do
          else
            mat_ref(ig1, ig2, :) = z1 * mt_rad**2 / g * besselj(1) * gc(:)
            mat_refl(ig1, ig2, :) = - z1 * mt_rad**2 / g * besselj(1) * (vpc + test_gvc(:, ig1))
            mat_refr(ig1, ig2, :) =   z1 * mt_rad**2 / g * besselj(1) * (vpc + test_gvc(:, ig2))
            do ipp = 1, 3
              mat_refm(ig1, ig2, :, ipp) = z1 * mt_rad**2 / g * besselj(1) * (vpc + test_gvc(:, ig1)) * (vpc(ipp) + test_gvc(ipp, ig2))
            end do
          end if
        end do
      end do
      mat_refm = mat_refm * cmplx( 0, -1, dp )

      ! compute matrix elements using module
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      do ip = 1, 3
        ! left and right gradient
        call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntog(:, :, 1, ip), zzero, mat_res )
        write( errmsg, '("Expected: Muffin-tin volume integrals agree with analytic result &
                          for component ",i1," of left and right wavefunction gradient.")' ) ip
        call test_report%assert( all_close( mat_ref(:, :, ip), mat_res, tol ), trim( errmsg ) )
        ! left gradient
        call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntolg(:, :, 1, ip), zzero, mat_res )
        write( errmsg, '("Expected: Muffin-tin volume integrals agree with analytic result &
                          for component ",i1," of left wavefunction gradient.")' ) ip
        call test_report%assert( all_close( mat_refl(:, :, ip), mat_res, tol ), trim( errmsg ) )
        ! right gradient
        call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntorg(:, :, 1, ip), zzero, mat_res )
        write( errmsg, '("Expected: Muffin-tin volume integrals agree with analytic result &
                          for component ",i1," of right wavefunction gradient.")' ) ip
        call test_report%assert( all_close( mat_refr(:, :, ip), mat_res, tol ), trim( errmsg ) )
        ! both-sided gradient
        do ipp = 1, 3
          call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntomg(:, :, 1, ip, ipp), zzero, mat_res )
          write( errmsg, '("Expected: Muffin-tin volume integrals agree with analytic result &
                            for components ",i1," and ",i1," of both-sided wavefunction gradient.")' ) ip, ipp
          call test_report%assert( all_close( mat_refm(:, :, ip, ipp), mat_res, tol ), trim( errmsg ) )
        end do
      end do

      deallocate( apwalm, mat_ref, mat_refl, mat_refr, mat_refm, mat_res )
    end subroutine me_test_grad

    !> Check the muffin-tin surface integral of the overlap of two wavefunction gradients.
    !>
    !> Compute the integrals
    !> \[ \begin{align*}
    !>       \langle \psi_{m{\bf p}} | \nabla^\top \psi_{n{\bf p}} \rangle_{\partial\alpha}
    !>       & \overset{?}{=} \frac{1}{\Omega} \oint\limits_{S_{{\bf \tau}_\alpha}(R_\alpha)}
    !>          {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf r}}\, \hat{\bf e}\, {\rm d}S \\
    !>       &= -\frac{4\pi}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf G}_n - {\bf G}_m) \cdot ({\bf p} + {\bf G}_n)^\top \;, 
    !>    \end{align*} \]
    !> and
    !> \[ \langle \nabla^\top \psi_{m{\bf p}} | \psi_{n{\bf p}} \rangle_{\partial\alpha}
    !>        = \frac{4\pi}{\Omega}\, {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf \tau}_\alpha}\,
    !>          \frac{R_\alpha^2}{|{\bf G}_n - {\bf G}_m|}\, j_1( R_\alpha |{\bf G}_n - {\bf G}_m|) 
    !>          ({\bf G}_n - {\bf G}_m) \cdot ({\bf p} + {\bf G}_m)^\top \;. \]
    subroutine me_test_sfgrad( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ip, ipp, ig1, ig2, gl(3)
      real(dp) :: gc(3), g, dotp, besselj(0:1)
      complex(dp) :: z1
      integer, target :: ngkmax
      character(256) :: errmsg

      complex(dp), allocatable :: mat_refl(:,:,:,:), mat_refr(:,:,:,:)
      complex(dp), allocatable :: mat_res(:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_refl(ng, ng, 3, 3), mat_refr(ng, ng, 3, 3) )
      allocate( mat_res(ng, ng) )

      ! set up reference
      do ig2 = 1, ng
        do ig1 = 1, ng
          gl = test_gvl(:, ig2) - test_gvl(:, ig1)
          gc = test_gvc(:, ig2) - test_gvc(:, ig1)
          g = norm2( gc )
          call sbessel( 1, g*mt_rad, besselj )
          dotp = twopi * dot_product( dble(gl), atposl )
          z1 = - fourpi * cmplx( cos(dotp), sin(dotp), dp ) / omega
          if( g < 1.d-12 ) then
            mat_refl(ig1, ig2, :, :) = zzero
            mat_refr(ig1, ig2, :, :) = zzero
          else
            do ipp = 1, 3
              mat_refl(ig1, ig2, :, ipp) = - z1 * mt_rad**2 / g * besselj(1) * (vpc + test_gvc(:, ig1)) * gc(ipp)
              mat_refr(ig1, ig2, :, ipp) =   z1 * mt_rad**2 / g * besselj(1) * (vpc + test_gvc(:, ig2)) * gc(ipp)
            end do
          end if
        end do
      end do

      ! compute matrix elements using module
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      do ip = 1, 3
        do ipp = 1, 3
          ! left gradient
          call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntolgs(:, :, 1, ip, ipp), zzero, mat_res )
          write( errmsg, '("Expected: Muffin-tin surface integrals agree with analytic result &
                            for component ",i1," of left wavefunction gradient &
                            and component ",i1," of surface normal vector.")' ) ip, ipp
          call test_report%assert( all_close( mat_refl(:, :, ip, ipp), mat_res, tol ), trim( errmsg ) )
          ! right gradient
          call me_mt_mat( 1, 1, Gkset%ngk(1, 1), apwalm(:, :, :, 1), evec, zone, rigntorgs(:, :, 1, ip, ipp), zzero, mat_res )
          write( errmsg, '("Expected: Muffin-tin surface integrals agree with analytic result &
                            for component ",i1," of right wavefunction gradient &
                            and component ",i1," of surface normal vector.")' ) ip, ipp
          call test_report%assert( all_close( mat_refr(:, :, ip, ipp), mat_res, tol ), trim( errmsg ) )
        end do
      end do

      deallocate( apwalm, mat_refl, mat_refr, mat_res )
    end subroutine me_test_sfgrad

    !> Check the orthogonality of the matrix elements of different plane-wave operators
    !> and two wavefunctions.
    !>
    !> Compute the integrals
    !> \[ \begin{align*}
    !>       \langle \psi_{m{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \psi_{0{\bf p}} \rangle
    !>       &\overset{?}{=} \frac{1}{\Omega} \int\limits_{\Omega}
    !>          {\rm e}^{{\rm i}({\bf G}_n - {\bf G}_m)\cdot{\bf r}}\, {\rm d}r^3 \\
    !>       &= \delta_{mn} \;.
    !>    \end{align*} \]
    subroutine me_test_operator( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ig
      integer, target :: ngkmax
      character(256) :: errmsg

      complex(dp), allocatable :: mat_ref(:,:), mat_res(:,:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_ref(ng, ng), source=zzero )
      allocate( mat_res(ng, 1, ng) )

      ! set up reference
      do ig = 1, ng
        mat_ref(ig, ig) = zone
      end do

      ! compute matrix elements using module
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) ) 
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      do ig = 1, ng
        call me_mt_mat( 1, 1, Gkset%ngk(1, 1), Gkset%ngk(1, 1), apwalm(:, :, :, 1), apwalm(:, :, :, 1), &
               evec, evec(:, 1:1), zone, rigntpv(:, :, 1, ig), zzero, mat_res(:, :, ig) )
        call me_ir_mat( Gkset, 1, Gkset, 1, &
               evec, evec(:, 1:1), zone, opig(:, ig), zone, mat_res(:, :, ig) )

        ! assertion
        write( errmsg, '("Expected: Matrix elements of plane-wave operator for &
                          test G-vector ",i1," are orthogonal.")' ) ig
        call test_report%assert( all_close( mat_ref(:, ig), mat_res(:, 1, ig), tol ), trim( errmsg ) )
      end do

      deallocate( apwalm, mat_ref, mat_res )
    end subroutine me_test_operator
    
    !> Check Gauss' theorem by comparing the muffin-tin volume integral of the gradient
    !> of different plane-wave operators and two wavefunctions against the respective
    !> muffin-tin surface integral.
    !>
    !> Compute the integrals
    !> \[ \langle \nabla \psi_{m{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \psi_{0{\bf p}} \rangle_{\alpha}
    !>    + \langle \psi_{m{\bf p}} | \nabla {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \psi_{0{\bf p}} \rangle_{\alpha}
    !>    + \langle \psi_{m{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \nabla \psi_{0{\bf p}} \rangle_{\alpha}
    !>    \overset{?}{=} \langle \psi_{m{\bf p}} | {\rm e}^{{\rm i}{\bf G}_n \cdot {\bf r}} | \psi_{0{\bf p}} \rangle_{\partial\alpha} \;. \]
    subroutine me_test_gauss( test_report )
      use mod_Gkvector, only: ngkmax_ptr
      use math_utils, only: all_close
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: ig, ip
      integer, target :: ngkmax
      character(256) :: errmsg

      complex(dp), allocatable :: mat_grad(:,:,:,:), mat_surf(:,:,:,:)
      complex(dp), allocatable :: apwalm(:,:,:,:)

      ngkmax = Gkset%ngkmax
      ngkmax_ptr => ngkmax

      allocate( mat_grad(ng, 1, 3, ng) )
      allocate( mat_surf(ng, 1, 3, ng) )

      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) ) 
      call match( Gkset%ngk(1, 1), Gkset%gkc(:, 1, 1), Gkset%tpgkc(:, :, 1, 1), Gkset%sfacgk(:, :, 1, 1), apwalm )
      do ig = 1, ng
        do ip = 1, 3
          ! compute muffin-tin volume integrals of the gradient
          call me_mt_mat( 1, 1, Gkset%ngk(1, 1), Gkset%ngk(1, 1), apwalm(:, :, :, 1), apwalm(:, :, :, 1), &
                 evec, evec(:, 1:1), zone, rigntpg(:, :, 1, (ig-1)*3+ip), zzero, mat_grad(:, :, ip, ig) )
          ! compute muffin-tin surface integrals
          call me_mt_mat( 1, 1, Gkset%ngk(1, 1), Gkset%ngk(1, 1), apwalm(:, :, :, 1), apwalm(:, :, :, 1), &
                 evec, evec(:, 1:1), zone, rigntps(:, :, 1, (ig-1)*3+ip), zzero, mat_surf(:, :, ip, ig) )

          ! assertion
          write( errmsg, '("Expected: Gauss theorem holds for test G-vector ",i1," &
                            and Cartesian direction ",i1,".")' ) ig, ip
          call test_report%assert( all_close( mat_grad(:, 1, ip, ig), mat_surf(:, 1, ip, ig), tol ), trim( errmsg ) )
        end do
      end do

      deallocate( apwalm, mat_grad, mat_surf )
    end subroutine me_test_gauss
end module
