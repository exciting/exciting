!> This module provides functionalities to compute the electrostatic potential
!> from a given charge density distribution, i.e., solving Poisson's equation,
!> using Weinert's method [1].
!>
!> [1]: [M. Weinert. Solution of Poisson's equation: Beyond Ewald-type methods. 
!> *J. Math. Phys.* **22**, 2433(1981)](https://doi.org/10.1063/1.524800)
module weinert
  use precision, only: dp

  implicit none
  private

  public :: surface_ir, multipoles_ir, poisson_ir
  public :: poisson_and_multipoles_mt, match_bound_mt

  contains

    !> This subroutine evaluates an interstitial function given by the Fourier series
    !> \[ f({\bf r}) = \sum_{\bf G} \hat{f}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
    !> on the surface of the muffin-tin spheres.
    !>
    !> The function evaluated on the surface of the muffin-tin sphere \(\alpha\) is given by
    !> \[ f^{{\rm SF},\alpha}(\theta,\phi) = f({\bf \tau}_\alpha + R_\alpha \hat{\bf r}(\theta,\phi)) = 
    !> \sum_{l,m} f^{{\rm SF},\alpha}_{lm} \, Y_{lm}(\theta,\phi) \;,\]
    !> with
    !> \[ f^{{\rm SF},\alpha}_{lm} = 4\pi \, {\rm i}^l \sum_{\bf G} {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}
    !> \hat{f}({\bf G+p}) \, j_l(|{\bf G+p}| R_\alpha) \, Y_{lm}^\ast({\bf \widehat{G+p}}) \;,\]
    !> where \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) which is centered at \({\bf \tau}_\alpha\)
    !> and \(j_l(x)\) are the spherical Bessel functions.
    subroutine surface_ir( lmax, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, fsf)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: rmt
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function
      complex(dp), intent(in) :: fig(:)
      !> function values on the muffin-tin sphere surfaces given by \(f^{{\rm SF},\alpha}_{lm}\)
      complex(dp), intent(out) :: fsf(:,:)

      integer :: is, ia, ias
    
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call surface_ir_single_mt( lmax, rmt(is), ngp, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
                                     fig, fsf(:,ias))
        end do
      end do
    end subroutine
    !> Same as [[surface_ir(subroutine)]] but for a single muffin-tin sphere.
    subroutine surface_ir_single_mt( lmax, rmt, ngp, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, fsf)
      use constants, only: zzero, fourpi
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function
      complex(dp), intent(in) :: fig(:)
      !> function values on the muffin-tin sphere surfaces given by \(f^{{\rm SF},\alpha}_{lm}\)
      complex(dp), intent(out) :: fsf(:)

      integer :: l, m, lm, igp, ifg, ig(3)
      complex(dp) :: z1, z2, zil

      real(dp), allocatable :: rl3(:)
    
      fsf = zzero
      
      allocate( rl3(0:lmax))
    
      rl3(0) = rmt**3
      do l = 1, lmax
        rl3(l) = rl3(l-1)*rmt
      end do
!$omp parallel default(shared) private(igp,ig,ifg,z1,z2,zil,l,m,lm) reduction(+:fsf)
!$omp do
      do igp = 1, ngp
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        z1 = fourpi*sfacgp(igp)*fig(ifg)
        zil = cmplx( 1._dp, 0._dp, dp)
        do l = 0, lmax
          z2 = z1*zil*jlgpr( l, igp)
          do m = -l, l
            lm = l*(l+1) + m + 1
            fsf(lm) = fsf(lm) + z2*conjg( ylmgp( lm, igp))
          end do
          zil = cmplx( -aimag(zil), dble(zil), dp)
        end do
      end do
!$omp end do
!$omp end parallel

      deallocate( rl3)
    end subroutine

    !> This subroutine calculates the multipole moments corresponding to the extension of
    !> an interstitial function given by the Fourier series
    !> \[ f({\bf r}) = \sum_{\bf G} \hat{f}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
    !> into the interior of the muffin-tin spheres.
    !>
    !> The multipole moments corresponding to the extension of the function inside the the muffin-tin sphere \(\alpha\)
    !> are given by
    !> \[ q^{{\rm IR},\alpha}_{lm} = 4\pi \, {\rm i}^l \, R_\alpha^{l+3} \sum_{\bf G} {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}
    !> \hat{f}({\bf G+p}) \, \frac{j_{l+1}(|{\bf G+p}| R_\alpha)}{|{\bf G+p}| R_\alpha} \, Y_{lm}^\ast({\bf \widehat{G+p}}) \;,\]
    !> where \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) which is centered at \({\bf \tau}_\alpha\)
    !> and \(j_l(x)\) are the spherical Bessel functions.
    subroutine multipoles_ir( lmax, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, qlm)
      use modinput
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: rmt
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
      complex(dp), intent(in) :: fig(:)
      !> multipole moments of the function's extension inside the muffin-tin spheres
      complex(dp), intent(out) :: qlm(:,:)

      integer :: is, ia, ias
    
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call multipoles_ir_single_mt( lmax, rmt(is), ngp, gpc, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
                                        fig, qlm(:,ias), epslat=input%structure%epslat)
        end do
      end do
    end subroutine
    !> Same as [[multipoles_ir(subroutine)]] but for a single muffin-tin sphere.
    subroutine multipoles_ir_single_mt( lmax, rmt, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, qlm, epslat)
      use constants, only: zzero, fourpi, y00
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
      complex(dp), intent(in) :: fig(:)
      !> multipole moments of the function's extension inside the muffin-tin spheres
      complex(dp), intent(out) :: qlm(:)
      !> threshold below which \(|{\bf G+p}|\) is considered zero
      real(dp), optional, intent(in) :: epslat

      integer :: i, l, m, lm, igp, ngpf, ngpz, ifg, ig(3)
      real(dp) :: eps
      complex(dp) :: z1, z2, zil, figzero

      integer, allocatable :: igp_finite(:), igp_zero(:)
      real(dp), allocatable :: rl3(:)
    
      eps = 1.d-12
      if( present( epslat)) eps = epslat

      qlm = zzero
      figzero = zzero
      
      allocate( rl3(0:lmax))
    
      rl3(0) = rmt**3
      do l = 1, lmax
        rl3(l) = rl3( l-1)*rmt
      end do

      igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > eps, i=1, ngp)])
      ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg,z1,z2,zil,l,m,lm) reduction(+:qlm)
!$omp do
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        z1 = fourpi*sfacgp(igp)*fig(ifg)/(gpc(igp)*rmt)
        zil = cmplx( 1._dp, 0._dp, dp)
        do l = 0, lmax
          z2 = z1*zil*rl3(l)*jlgpr( l+1, igp)
          do m = -l, l
            lm = l*(l+1) + m + 1
            qlm(lm) = qlm(lm) + z2*conjg( ylmgp( lm, igp))
          end do
          zil = cmplx( -aimag(zil), dble(zil), dp)
        end do
      end do
!$omp end do
!$omp end parallel
      igp_zero = pack( [(i, i=1, ngp)], [(gpc(i) <= eps, i=1, ngp)])
      ngpz = size( igp_zero)
      do i = 1, ngpz
        igp = igp_zero(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        qlm(1) = qlm(1) + fourpi/3._dp*rl3(0)*y00*fig(ifg)
      end do

      deallocate( rl3)
    end subroutine

    !> This subroutine solves Poisson's equation for a given complex charge density
    !> contained in an isolated muffin-tin sphere using the Green's function approach
    !> and calculates the multipole moments of the charge distribution.
    !>
    !> The general electrostatic potential arising from an arbitrary charge distribution
    !> \[ n({\bf r}) = n({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} n^\alpha_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> inside the muffin-tin sphere \(\alpha\) is given by
    !> \[ v_{\rm sph}({\bf r}) = v_{\rm sph}({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} v_{\rm sph}[n^\alpha_{lm}](r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> with
    !> \[ v_{\rm sph}[n^\alpha_{lm}](r) = \frac{4\pi}{2l+1} \left[ 
    !> \frac{1}{r^{l+1}} \int_0^r s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s + 
    !> r^l \int_r^{R_\alpha} \frac{n^\alpha_{lm}(s)}{s^{l-1}} \, {\rm d}s - 
    !> \frac{r^l}{R_\alpha^{2l+1}} \int_0^{R_\alpha} s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s \right] \;.\]
    !>
    !> In addition, the multipole moments of the charge distribution in the sphere are calculated as
    !> \[ q^{{\rm MT},\alpha}_{lm} = \int_0^{R_\alpha} s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s \,. \]
    !>
    !> @note \(v_{\rm sph}({\rm r})\) as defined above goes to zero on the muffin-tin surface.@endnote
    subroutine poisson_and_multipoles_mt( lmax, nr, r, zrhomt, zvclmt, qlm)
      use constants, only: fourpi
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> number of radial grid points
      integer, intent(in) :: nr
      !> radial grid
      real(dp), intent(in) :: r(:)
      !> complex charge distribution \(n^\alpha_{lm}(r)\)
      complex(dp), intent(in) :: zrhomt(:,:)
      !> complex electrostatic potential \(v_{\rm sph}[n^\alpha_{lm}](r)\)
      complex(dp), intent(out) :: zvclmt(:,:)
      !> multipole moments of the charge distribution \(q^{{\rm MT},\alpha}_{lm}\)
      complex(dp), intent(out) :: qlm(:)

      integer :: l, m, lm, ir
      real(dp) :: t1, t2, t3, t4, t5

      real(dp), allocatable :: ri(:), rl(:), ril1(:), cf(:,:)
      real(dp), allocatable :: fr(:,:), gr(:,:)

      allocate( ri(nr), rl(nr), ril1(nr), cf(3,nr))
      allocate( fr(nr,4), gr(nr,4))
    
      ! initialise r^l and r^(-l-1)
      do ir = 1, nr
        rl(ir) = 1._dp
        ri(ir) = 1._dp/r(ir)
        ril1(ir) = ri(ir)
      end do
      lm = 0
      do l = 0, lmax
        t1 = fourpi/dble( 2*l+1)
        do m = -l, l
          lm = lm + 1
          do ir = 1, nr
            t2 = rl(ir)*r(ir)*r(ir)    ! r^(2+l)
            t3 = ril1(ir)*r(ir)*r(ir)  ! r^(1-l)
            t4 = dble( zrhomt(lm, ir))
            t5 = aimag( zrhomt(lm, ir))
            fr(ir,1) = t2*t4
            fr(ir,2) = t2*t5
            fr(ir,3) = t3*t4
            fr(ir,4) = t3*t5
          end do
          call fderiv( -1, nr, r, fr(:,1), gr(:,1), cf)
          call fderiv( -1, nr, r, fr(:,2), gr(:,2), cf)
          call fderiv( -1, nr, r, fr(:,3), gr(:,3), cf)
          call fderiv( -1, nr, r, fr(:,4), gr(:,4), cf)
          qlm(lm) = cmplx( gr(nr,1), gr(nr,2), 8)
          do ir = 1, nr
            t2 = ril1(ir)*gr(ir,1) + rl(ir)*(gr(nr,3) - gr(ir,3) - gr(nr,1)*ril1(nr)/rl(nr))
            t3 = ril1(ir)*gr(ir,2) + rl(ir)*(gr(nr,4) - gr(ir,4) - gr(nr,2)*ril1(nr)/rl(nr))
            zvclmt( lm, ir) = t1*cmplx( t2, t3, 8)
          end do
        end do
        ! update r^l and r^(-l-1)
        if( l < lmax) then
          do ir = 1, nr
            rl(ir) = rl(ir)*r(ir)
            ril1(ir) = ril1(ir)*ri(ir)
          end do
        end if
      end do
    
      deallocate( ri, rl, ril1, cf, fr, gr)
      return
    end subroutine

    !> This subroutine solves Poisson's equation for a complex charge density given in the
    !> interstitial region by 
    !> \[ n({\bf r}) = \sum_{\bf G} \hat{n}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
    !> and having multipole moments \(q^\alpha_{lm}\) inside the muffin-tin spheres.
    !>
    !> This is done by constructing a pseudodensity with the given multipole moments and having a 
    !> rapidly converging Fourier series. The Fourier components of the pseudodensity are given by
    !> \[ \hat{n}^{\rm ps}({\bf G+p}) = \hat{n}({\bf G+p}) + \sum_\alpha \hat{n}^{{\rm ps},\alpha}({\bf G+p}) \;.\]
    !> See [[pseudodensity_ir_single_mt(subroutine)]] for further details.
    !>
    !> From the pseudodensity the electrostatic potential in the interstitial region is obtained
    !> by solving Poisson's equation in reciprocal space, i.e., 
    !> \[ V({\bf r}) = \sum_{\bf G} \hat{V}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \;, \]
    !> with
    !> \[ \hat{V}({\bf G+p}) = 4\pi \frac{\hat{n}^{\rm ps}({\bf G+p})}{|{\bf G+p}|^2} \;.\]
    subroutine poisson_ir( lmax, npsden, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, zrhoig, qlm, zvclig)
      use modinput
      use mod_lattice, only: omega
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: rmt
      use constants, only: zzero, fourpi
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> pseudodensity expansion order
      integer, intent(in) :: npsden
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> on entry: Fourier components \(\hat{n}({\bf G+p})\) of the interstitial charge density on the FFT grid
      !> on exit: Fourier components \(\hat{n}^{\rm ps}({\bf G+p})\) of the pseudodensity on the FFT grid
      complex(dp), intent(inout) :: zrhoig(:)
      !> muffin-tin multipole moments \(q^\alpha_{lm}\) of the charge density
      complex(dp), intent(in) :: qlm(:,:)
      !> Fourier components \(\hat{V}({\bf G+p})\) of the interstitial electrostatic potential on the FFT grid
      complex(dp), intent(out) :: zvclig(:)

      integer :: i, is, ia, ias, igp, ngpf, ifg, ig(3)

      integer, allocatable :: igp_finite(:)
    
      ! add Fourier components of pseudodensity from multipole moments
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call pseudodensity_ir_single_mt( lmax, rmt(is), omega, npsden, ngp, gpc, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
                                           zrhoig, qlm(:,ias), epslat=input%structure%epslat)
        end do
      end do
  
      ! solve Poisson's equation in reciprocal space
      zvclig = zzero
      igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > input%structure%epslat, i=1, ngp)])
      ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg)
!$omp do
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        zvclig(ifg) = fourpi*zrhoig(ifg)/(gpc(igp)**2)
      end do
!$omp end do
!$omp end parallel
    end subroutine

    !> This subroutine computes the Fourier components of a quickly converging pseudodensity with multipole moments \(q_{lm}\)
    !> in muffin-tin \(\alpha\) and adds them to the input argument `zrhoig`.
    !>
    !> The Fourier components to be added are given by
    !> \[ \hat{n}^{{\rm ps},\alpha}({\bf G+p}) = \frac{4\pi}{\Omega} {\rm e}^{-{\rm i}({\bf G+p}) \cdot {\bf \tau}_\alpha}
    !> \sum_{l,m} \left( -\frac{\rm i}{R_\alpha} \right)^l \frac{(2l + 2N + 3)!!}{(2l + 1)!!} 
    !> \frac{j_{l+N+1}(|{\bf G+p}|R_\alpha)}{(|{\bf G+p}|R_\alpha)^{N+1}} \, q^\alpha_{lm} \, Y_{lm}(\widehat{\bf G+p}) \;,\]
    !> where \(\Omega\) is the unit cell volume, \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) 
    !> which is centered at \({\bf \tau}_\alpha\) and \(j_l(x)\) are the spherical Bessel functions.
    !> \(N\) is the expansion order of the pseudodensity.
    subroutine pseudodensity_ir_single_mt( lmax, rmt, omega, npsden, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, zrhoig, qlm, epslat)
      use constants, only: zzero, fourpi, zi, y00
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> unit cell volume \(\Omega\)
      real(dp), intent(in) :: omega
      !> pseudodensity expansion order
      integer, intent(in) :: npsden
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{n}^{\rm ps}({\bf G+p})\) of the pseudodensity on the FFT grid
      complex(dp), intent(inout) :: zrhoig(:)
      !> muffin-tin multipole moments \(q^\alpha_{lm}\) of the charge density
      complex(dp), intent(in) :: qlm(:)
      !> threshold below which \(|{\bf G+p}|\) is considered zero
      real(dp), optional, intent(in) :: epslat

      integer :: i, l, m, lm, igp, ngpf, ngpz, ifg, ig(3)
      real(dp) :: eps, fpo, t1, t2
      complex(dp) :: z1, z2
    
      integer, allocatable :: igp_finite(:), igp_zero(:)
      complex(dp), allocatable :: zrp(:)

      real(dp), external :: factnm
      
      eps = 1.d-12
      if( present( epslat)) eps = epslat

      fpo = fourpi/omega
      allocate( zrp( (lmax+1)**2))
    
      ! add Fourier components of pseudodensity from multipole moments
      t1 = 1.d0
      do l = 0, lmax
        t2 = t1*factnm( 2*l + 2*npsden + 3, 2)/factnm( 2*l+1, 2)
        t1 = t1/rmt
        do m = -l, l
          lm = l*(l+1) + m + 1
          zrp(lm) = conjg( zi**l)*t2*qlm(lm)
        end do
      end do
    
      igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > eps, i=1, ngp)])
      ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg,l,m,lm,z1,z2,t1,t2)
!$omp do
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        z1 = fpo*conjg( sfacgp(igp))
        t1 = (gpc(igp)*rmt)**(npsden + 1)
        do l = 0, lmax
          t2 = jlgpr( l+npsden+1, igp)/t1
          z2 = zzero
          do m = -l, l
            lm = l*(l+1) + m + 1
            z2 = z2 + zrp(lm)*ylmgp( lm, igp)
          end do
          zrhoig(ifg) = zrhoig(ifg) + z1*t2*z2
        end do
      end do
!$omp end do
!$omp end parallel
      igp_zero = pack( [(i, i=1, ngp)], [(gpc(i) <= eps, i=1, ngp)])
      ngpz = size( igp_zero)
      do i = 1, ngpz
        igp = igp_zero(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        zrhoig(ifg) = zrhoig(ifg) + fpo/factnm( 2*npsden + 3, 2)*zrp(1)*y00
      end do
      deallocate( zrp)
    end subroutine

    !> Given a complex muffin-tin function
    !> \[ f({\bf r}) = f({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} f^\alpha_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> and an interstitial function given on the muffin-tin sphere surface by
    !> \(f^{{\rm SF},\alpha}_{lm}\) according to subroutine [[surface_ir(subroutine)]],
    !> this subroutine matches the muffin-tin function to the interstial function by
    !> adding the homogeneous function
    !> \[ f^{\rm hom}({\bf r}) = f^{\rm hom}({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} \left[ f^{{\rm SF},\alpha}_{lm} - f^\alpha_{lm}(R_\alpha) \right] \,
    !> \left( \frac{r_\alpha}{R_\alpha} \right)^l \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> to the muffin-tin function.
    subroutine match_bound_mt( lmax, nr, r, rmt, firsf, fmt)
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> number of radial grid points
      integer, intent(in) :: nr
      !> radial grid
      real(dp), intent(in) :: r(:)
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> interstitial function on muffin-tin sphere surface \(f^{{\rm SF},\alpha}_{lm}\)
      complex(dp), intent(in) :: firsf(:)
      !> muffin-tin function \(f^\alpha_{lm}(r)\)
      complex(dp), intent(inout) :: fmt(:,:)

      integer :: l, m, lm, ir
      complex(dp) :: df

      real(dp), allocatable :: rr(:,:)

      allocate( rr(nr,2), source=1._dp)

      do ir = 1, nr
        rr(ir,2) = r(ir)/rmt
      end do

      lm = 0
      do l = 0, lmax
        do m = -l, l
          lm = lm + 1
          df = firsf(lm) - fmt(lm,nr)
          do ir = 1, nr
            fmt(lm,ir) = fmt(lm,ir) + df*rr(ir,1)
          end do
        end do
        do ir = 1, nr
          rr(ir,1) = rr(ir,1)*rr(ir,2)
        end do
      end do

      deallocate( rr)
    end subroutine

end module weinert
