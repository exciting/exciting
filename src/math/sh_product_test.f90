module sh_product_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  
  implicit none
  private

  !> effective zero for comparisons
  real(dp), parameter :: TOL = 1.e-16_dp
  !> sphere radius
  real(dp), parameter :: RAD = 2._dp
  !> number of radial points
  integer, parameter :: NUMRAD = 500
  !> maximum \(l\) for product of test functions
  integer, parameter :: LMAX = 9
  !> effective infinite \(l\)
  integer, parameter :: LMAXINF = 20
  !> maximum number of \((l,m)\) pairs
  integer, parameter :: LMMAX = (LMAX + 1)**2
  integer, parameter :: LMMAXINF = (LMAXINF + 1)**2
  !> wavevector of test plane wave
  real(dp), parameter :: VKC(3) = [0.391_dp, -0.482_dp, 0.168_dp]
  !> maximum \(l\) for test functions
  integer :: lmaxfun
  !> radial grid
  real(dp), allocatable :: rgrid(:)
  !> real spherical harmonics test functions
  real(dp) :: dfun(LMMAX,NUMRAD,2)
  !> complex spherical harmonics test function
  complex(dp) :: zfun(LMMAX,NUMRAD)
  !> real spherical harmonics reference functions
  real(dp) :: dfun_ref(LMMAX,NUMRAD,2)
  !> complex spherical harmonics reference function
  complex(dp) :: zfun_ref(LMMAX,NUMRAD)
  !> real forward spherical harmonics transform matrix
  real(dp), allocatable :: dfsht(:,:)
  !> real backward spherical harmonics transform matrix
  real(dp), allocatable :: dbsht(:,:)
  !> complex forward spherical harmonics transform matrix
  complex(dp), allocatable :: zfsht(:,:)
  !> complex backward spherical harmonics transform matrix
  complex(dp), allocatable :: zbsht(:,:)

  public :: sh_product_test_driver

  contains

    !> Run tests for spherical harmonics products
    subroutine sh_product_test_driver( mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure

      !> test object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 14

      ! Initialize test object
      call test_report%init( n_assertions, mpiglobal)

      ! Run the test with expansion cutoff for test function
      ! of LMAX
      lmaxfun = LMAX
      ! initialize tests
      call init_test
      ! test subroutine dshmul
      call product_test_real( test_report)
      ! test subroutine zshmul and zshmulc
      call product_test_complex( test_report)
      ! test real-space product
      call real_space_test( test_report)

      ! Run the test with expansion cutoff for test function
      ! of LMAX/2
      lmaxfun = LMAX/2
      ! initialize tests
      call init_test
      ! test subroutine dshmul
      call product_test_real( test_report)
      ! test subroutine zshmul and zshmulc
      call product_test_complex( test_report)
      ! test real-space product
      call real_space_test( test_report)
  
      ! report results
      if (present(kill_on_failure)) then
        call test_report%report('sh_product', kill_on_failure)
      else
        call test_report%report('sh_product')
      end if

      ! Finalise test object
      call test_report%finalise()
    end subroutine sh_product_test_driver

    !> generate radial grid, test and reference functions, 
    !> and spherical harmonics transform matrices
    !>
    !> The complex test function is the spherical harmonics expansion of a plane wave
    !> with wavevector \({\bf k}\)
    !> \[ f_{\rm test}({\bf r}) 
    !>    = 4\pi \sum_{l=0}^{l_{\rm max}^{\rm fun}} \sum_{m=-l}^l
    !>      {\rm i}^l\, j_{l}(kr)\, Y_{lm}^\ast(\hat{\bf k})\, Y_{lm}(\hat{\bf r}) 
    !>    \approx {\rm e}^{{\rm i}{\bf k}\cdot{\bf r}} \;. \]
    !> The real test functions are its real and imaginary part.
    !> The expansion is truncated at `lmaxfun`.
    !>
    !> The complex reference function is the product of the test function with itself,
    !> i.e.,
    !> \[ f_{\rm ref}({\bf r}) = f_{\rm test}^2({\bf r}) \;. \]
    !> This product is again a spherical harmonics expansion. 
    !> This is truncated at `LMAX`.
    !> To find the exact result for \(f_{\rm ref}\), we write
    !> \(f_{\rm test}({\bf r}) = {\rm e}^{{\rm i}{\bf k}\cdot{\bf r}}
    !> - f_{\rm res}({\bf r})\),
    !> where \(f_{\rm res}\) is the residue from the plane wave expansion
    !> with spherical harmonics of order \(l > l_{\rm max}^{\rm fun}\).
    !> With that, we find
    !> \[ f_{\rm ref}({\bf r}) = {\rm e}^{2{\rm i}{\bf k}\cdot{\bf r}}
    !>    - 2\,{\rm e}^{{\rm i}{\bf k}\cdot{\bf r}}\, f_{\rm res}({\bf r})
    !>    + f_{\rm res}^2({\bf r}) \;. \]
    !> Again, the real reference functions are the real and imaginary part of the complex one.
    subroutine init_test
      use constants, only: fourpi, zi
      use grid_utils, only: linspace
      use mod_SHT, only: gen_rshtmat, gen_zshtmat

      integer :: l, m, lm, ir
      real(dp) :: kc, tp(2)
      complex(dp) :: z1
      real(dp), allocatable :: jl(:,:,:)
      complex(dp), allocatable :: ylm(:), exact(:,:), residue(:,:)

      allocate( jl(0:LMAXINF,NUMRAD,2))
      allocate( ylm(LMMAXINF), exact(LMMAXINF,NUMRAD), residue(LMMAXINF,NUMRAD))

      ! generate radial grid
      rgrid = linspace( 0._dp, RAD, NUMRAD)

      ! generate spherical Bessel functions and spherical harmonics
      call sphcrd( VKC, kc, tp)
      do ir = 1, NUMRAD
        call sbessel( LMAXINF, kc*rgrid(ir), jl(:,ir,1))
        call sbessel( LMAX, 2*kc*rgrid(ir), jl(:,ir,2))
      end do
      call genylm( LMAXINF, tp, ylm)

      ! generate test functions
      dfun = 0._dp; dfun_ref = 0._dp
      zfun = 0._dp; zfun_ref = 0._dp
      residue = 0._dp
      z1 = cmplx(fourpi, 0._dp, dp)
      lm = 0
      do l = 0, lmaxfun
        do m = -l, l
          lm = lm + 1
          do ir = 1, NUMRAD
            exact(lm,ir) = z1 * jl(l,ir,1) * conjg( ylm(lm))
            zfun(lm,ir) = exact(lm,ir)
            if( l <= LMAX) zfun_ref(lm,ir) = z1 * jl(l,ir,2) * conjg( ylm(lm))
          end do
        end do
        z1 = z1 * zi
      end do
      do l = lmaxfun + 1, LMAXINF
        do m = -l, l
          lm = lm + 1
          do ir = 1, NUMRAD
            exact(lm,ir) = z1 * jl(l,ir,1) * conjg( ylm(lm))
            residue(lm,ir) = exact(lm,ir)
            if( l <= LMAX) zfun_ref(lm,ir) = z1 * jl(l,ir,2) * conjg( ylm(lm))
          end do
        end do
        z1 = z1 * zi
      end do
      call zshmul( LMAXINF, LMAXINF, LMAX, NUMRAD, cmplx( -2._dp, 0._dp, dp), &
             exact, LMMAXINF, &
             residue, LMMAXINF, cmplx( 1._dp, 0._dp, dp), &
             zfun_ref, LMMAX)
      call zshmul( LMAXINF, LMAXINF, LMAX, NUMRAD, cmplx( 1._dp, 0._dp, dp), &
             residue, LMMAXINF, &
             residue, LMMAXINF, cmplx( 1._dp, 0._dp, dp), &
             zfun_ref, LMMAX)

      do ir = 1, NUMRAD
        call decompose_zflm( lmaxfun, zfun(:,ir), dfun(:,ir,1), dfun(:,ir,2))
        call decompose_zflm( LMAX, zfun_ref(:,ir), dfun_ref(:,ir,1), dfun_ref(:,ir,2))
      end do

      ! generate spherical harmonics transform matrices
      call gen_rshtmat( LMAX, dfsht, dbsht)
      call gen_zshtmat( LMAX, zfsht, zbsht)
    end subroutine init_test

    !> test subroutine [[dshmul(subroutine)]]
    !>
    !> Here, we compute the differences
    !> \[ D_r({\bf r}) = \Re f_{\rm ref}({\bf r}) - 
    !>      \left( \Re f_{\rm test}^2({\bf r}) - \Im f_{\rm test}^2({\bf r}) \right) \]
    !> and 
    !> \[ D_i({\bf r}) = \Im f_{\rm ref}({\bf r}) - 
    !>      2\, \Re f_{\rm test}({\bf r})\, \Im f_{\rm test}({\bf r}) \;, \]
    !> which should be zero.
    !> The error is evaluated as the integrals over the expansion sphere
    !> \[ \Delta_r = \int_{\rm sphere} D_r^2({\bf r})\, {\rm d}{\bf r}\; \text{ and }\;
    !>    \Delta_i = \int_{\rm sphere} D_i^2({\bf r})\, {\rm d}{\bf r} \; . \]
    !> Note: Both errors should be zero for all values of `lmaxfun` and `LMAX`
    !> given a sufficiently large `LMAXINF`.
    subroutine product_test_real( test_report)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: lm
      real(dp) :: error(2), f(NUMRAD), g(NUMRAD), cf(3,NUMRAD)
      real(dp) :: dprod_gnt(LMMAX,NUMRAD,2)

      ! real part of product
      call dshmul( lmaxfun, lmaxfun, LMAX, NUMRAD, 1._dp, &
             dfun(:,:,1), LMMAX, &
             dfun(:,:,1), LMMAX, 0._dp, &
             dprod_gnt(:,:,1), LMMAX)
      call dshmul( lmaxfun, lmaxfun, LMAX, NUMRAD, -1._dp, &
             dfun(:,:,2), LMMAX, &
             dfun(:,:,2), LMMAX, 1._dp, &
             dprod_gnt(:,:,1), LMMAX)
      ! imaginary part of product
      call dshmul( lmaxfun, lmaxfun, LMAX, NUMRAD, 2._dp, &
             dfun(:,:,1), LMMAX, &
             dfun(:,:,2), LMMAX, 0._dp, &
             dprod_gnt(:,:,2), LMMAX)

      error = 0._dp
      do lm = 1, LMMAX
        ! real part
        f = dprod_gnt(lm,:,1) - dfun_ref(lm,:,1)
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(1) = error(1) + g(NUMRAD)
        ! imaginary part
        f = dprod_gnt(lm,:,2) - dfun_ref(lm,:,2)
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(2) = error(2) + g(NUMRAD)
      end do

      ! this error should always be zero
      call test_report%assert( error(1) < TOL, &
             'Test product of real spherical harmonics expansions (real part) &
             using Gaunt coefficients. &
             Expected: 0.0')
      ! this error should always be zero
      call test_report%assert( error(2) < TOL, &
             'Test product of real spherical harmonics expansions (imaginary part) &
             using Gaunt coefficients. &
             Expected: 0.0')
    end subroutine product_test_real

    !> test subroutines [[zshmul(subroutine)]] and [[zshmulc(subroutine)]]
    !>
    !> Here, we compute the differences
    !> \[ D_c({\bf r}) = f_{\rm ref}({\bf r}) - f_{\rm test}^2({\bf r}) \]
    !> and
    !> \[ D_{cc}({\bf r}) = 1 - f_{\rm test}^\ast({\bf r})\, f_{\rm test}({\bf r}) \;, \]
    !> which should be zero.
    !> The error is evaluated as the integrals over the expansion sphere
    !> \[ \Delta_c = \int_{\rm sphere} |D_c({\bf r})|^2\, {\rm d}{\bf r}\; \text{ and }\;
    !>    \Delta_{cc} = \int_{\rm sphere} |D_{cc}({\bf r})|^2\, {\rm d}{\bf r} \; . \]
    !> Note: The error \(\Delta_c\) should be zero for all values of `lmaxfun` and `LMAX`
    !> given a sufficiently large `LMAXINF`. The error \(\Delta_{cc}\) converges to 
    !> zero for large `lmaxfun`.
    subroutine product_test_complex( test_report)
      use constants, only: y00
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: lm
      real(dp) :: error(2), f(NUMRAD), g(NUMRAD), cf(3,NUMRAD)
      complex(dp) :: zprod_gnt(LMMAX,NUMRAD)

      ! complex product
      call zshmul( lmaxfun, lmaxfun, LMAX, NUMRAD, cmplx( 1._dp, 0._dp, dp), &
             zfun, LMMAX, &
             zfun, LMMAX, cmplx( 0._dp, 0._dp, dp), &
             zprod_gnt, LMMAX)

      error = 0._dp
      do lm = 1, LMMAX
        ! complex
        f = abs( zprod_gnt(lm,:) - zfun_ref(lm,:))
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(1) = error(1) + g(NUMRAD)
      end do

      ! complex conjugate product
      call zshmulc( lmaxfun, lmaxfun, LMAX, NUMRAD, cmplx( 1._dp, 0._dp, dp), &
             zfun, LMMAX, &
             zfun, LMMAX, cmplx( 0._dp, 0._dp, dp), &
             zprod_gnt, LMMAX)
      zprod_gnt(1,:) = zprod_gnt(1,:) - 1._dp/y00
      do lm = 1, LMMAX
        ! complex
        f = abs( zprod_gnt(lm,:))
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(2) = error(2) + g(NUMRAD)
      end do

      ! this error should always be zero
      call test_report%assert( error(1) < TOL, &
             'Test product of complex spherical harmonics expansions &
             using Gaunt coefficients. &
             Expected: 0.0')
      ! this error should be zero for large enough lmaxfun
      if( lmaxfun >= 9) then
        call test_report%assert( error(2) < TOL, &
               'Test product of complex spherical harmonics expansions (conjugated) &
               using Gaunt coefficients. &
               Expected: 0.0')
      else
        call test_report%assert( error(2) >= TOL, &
               'Test product of complex spherical harmonics expansions (conjugated) &
               using Gaunt coefficients. &
               Expected: >0.0')
      end if
    end subroutine product_test_complex

    !> demonstrate necessity of subroutine [[dshmul(subroutine)]], 
    !> [[zshmul(subroutine)]] and [[zshmulc(subroutine)]]
    !>
    !> Here, we show that evaluating products of spherical harmonic
    !> expansions in real space is inaccurate when done as usual.
    !>
    !> The spherical harmonics transform (SHT) transforms a function
    !> given by its expansion coefficients to a function given by
    !> its values on a real-space grid determined by a fixed set of 
    !> pairs \((\theta_i,\phi_i)\)
    !> \[ F({\bf r}_i) = \sum_{l=0}^{l_{\rm max}} \sum_{m=-l}^l f_{lm}(r)\, Y_{lm}(\theta_i,\phi_i) \;, \]
    !> which in can be expressed as a matrix-vector product
    !> \[ {\bf F} = {\bf Y} \cdot {\bf f} \;, \]
    !> where \({\bf Y}\) is a square \((l_{\rm max} + 1)^2\) matrix
    !> with the values of \(Y_{lm}(\theta_i,\phi_i)\).
    !> 
    !> The standard approach to compute the expansion coefficients
    !> \(h_{lm}(r) \equiv {\bf h}\) of the product of two spherical
    !> harmonic expansions \(h({\bf r}) = f({\bf r})\, g({\bf r})\)
    !> is to first transform \(f\) and \(g\) to real space
    !> \[ {\bf F} = {\bf Y} \cdot {\bf f}\; \text{ and }\; 
    !>    {\bf G} = {\bf Y} \cdot {\bf g} \;, \]
    !> taking the pointwise product in real space
    !> \[ {\bf H} = {\bf F} \circ {\bf G} \;, \]
    !> and obtaining the expansion coefficients from the inverse transform
    !> \[ {\bf h} = {\bf Y}^{-1} \cdot {\bf H} \;. \]
    !> This approach, however, only gives the correct expansion coefficients
    !> if the sum of the order of the largest non-zero coefficients
    !> in \(f\) and \(g\) is not larger than \(l_{\rm max}\), i.e., for
    !> \[ l_{\rm max}^f + l_{\rm max}^g \leq l_{\rm max} \;. \]
    subroutine real_space_test( test_report)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      integer :: lm
      real(dp) :: error(3), f(NUMRAD), g(NUMRAD), cf(3,NUMRAD)
      real(dp) :: dprod_rsp(LMMAX,NUMRAD,2), dtmp(LMMAX,NUMRAD,2)
      complex(dp) :: zprod_rsp(LMMAX,NUMRAD), ztmp(LMMAX,NUMRAD)

      ! real part in real space
      call dgemm( 'n', 'n', LMMAX, NUMRAD, LMMAX, 1._dp, &
             dbsht, LMMAX, &
             dfun(:,:,1), LMMAX, 0._dp, &
             dprod_rsp(:,:,1), LMMAX)
      ! imaginary part in real space
      call dgemm( 'n', 'n', LMMAX, NUMRAD, LMMAX, 1._dp, &
             dbsht, LMMAX, &
             dfun(:,:,2), LMMAX, 0._dp, &
             dprod_rsp(:,:,2), LMMAX)
      dtmp(:,:,1) = dprod_rsp(:,:,1)**2 - dprod_rsp(:,:,2)**2
      dtmp(:,:,2) = 2._dp * dprod_rsp(:,:,1) * dprod_rsp(:,:,2)
      ! real part of product
      call dgemm( 'n', 'n', LMMAX, NUMRAD, LMMAX, 1._dp, &
             dfsht, LMMAX, &
             dtmp(:,:,1), LMMAX, 0._dp, &
             dprod_rsp(:,:,1), LMMAX)
      ! imaginary part of product
      call dgemm( 'n', 'n', LMMAX, NUMRAD, LMMAX, 1._dp, &
             dfsht, LMMAX, &
             dtmp(:,:,2), LMMAX, 0._dp, &
             dprod_rsp(:,:,2), LMMAX)
      ! complex function in real space
      call zgemm( 'n', 'n', LMMAX, NUMRAD, LMMAX, cmplx( 1._dp, 0._dp, dp), &
             zbsht, LMMAX, &
             zfun, LMMAX, cmplx( 0._dp, 0._dp, dp), &
             zprod_rsp, LMMAX)
      ztmp = zprod_rsp * zprod_rsp
      ! complex product
      call zgemm( 'n', 'n', LMMAX, NUMRAD, LMMAX, cmplx( 1._dp, 0._dp, dp), &
             zfsht, LMMAX, &
             ztmp, LMMAX, cmplx( 0._dp, 0._dp, dp), &
             zprod_rsp, LMMAX)

      error = 0._dp
      do lm = 1, LMMAX
        ! real part
        f = dprod_rsp(lm,:,1) - dfun_ref(lm,:,1)
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(1) = error(1) + g(NUMRAD)
        ! imaginary part
        f = dprod_rsp(lm,:,2) - dfun_ref(lm,:,2)
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(2) = error(2) + g(NUMRAD)
        ! complex
        f = abs( zprod_rsp(lm,:) - zfun_ref(lm,:))
        f = f**2 * rgrid**2
        call fderiv( -1, NUMRAD, rgrid, f, g, cf)
        error(3) = error(3) + g(NUMRAD)
      end do

      ! this errors should be zero only for lmaxfun <= LMAX/2
      if( lmaxfun <= LMAX/2) then
        call test_report%assert( error(1) < TOL, &
               'Test product of real spherical harmonics expansions (real part) &
               using real-space representation. &
               Expected: 0.0')
        call test_report%assert( error(2) < TOL, &
               'Test product of real spherical harmonics expansions (imaginary part) &
               using real-space representation. &
               Expected: 0.0')
        call test_report%assert( error(3) < TOL, &
               'Test product of complex spherical harmonics expansions &
               using real-space representation. &
               Expected: 0.0')
      else
        call test_report%assert( error(1) >= TOL, &
               'Test product of real spherical harmonics expansions (real part) &
               using real-space representation. &
               Expected: >0.0')
        call test_report%assert( error(2) >= TOL, &
               'Test product of real spherical harmonics expansions (imaginary part) &
               using real-space representation. &
               Expected: >0.0')
        call test_report%assert( error(3) >= TOL, &
               'Test product of complex spherical harmonics expansions &
               using real-space representation. &
               Expected: >0.0')
      end if
    end subroutine real_space_test
end module
