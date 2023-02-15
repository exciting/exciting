!> title:   (L)APW+LO matrix elements
!> summary: Module that implements procedures for the 
!>          calculation of matrix elements of local operators
!>          and (L)APW+LO basis functions / wavefunctions.
!> author:  Sebastian Tillack
!> date:    August 2022
!> licence: GPL
!>
!> This module contains the core procedures underlying the calculation
!> of matrix elements within the (L)APW+LO basis.
!> See [[matrix_elements(module)]] for usage and further documentation.
module matrix_elements_lapw_lo
  use precision, only: dp
  use asserts, only: assert
  use constants, only: zzero, zone
  use modmpi
  use muffin_tin_basis, only: mt_basis_type, generate_non_zero_clebsch_gordan
  use mod_kpointset, only: G_set

  implicit none
  private

  ! local helper arrays
  
  ! MUFFIN TIN BASIS
  !> maximum \(l\) for basis functions
  integer :: lmaxb
  !> maximum \(l\) for operator expansion
  integer :: lmaxo
  !> muffin-tin basis object
  type(mt_basis_type), public :: basis
  !> map from \(l\) and \(m\) to combined \((l,m)\) index
  integer, allocatable :: lm_join(:,:)
  !> map from combined \((l,m)\) index to \(l\) and \(m\)
  integer, allocatable :: lm_split(:,:)
  !> number of non-zero Clebsch-Gordan coefficients for a given \((l,m)\) pair
  integer, allocatable :: cg_num(:,:,:)
  !> \((l',m')\) pair of non-zero Clebsch-Gordan coefficients for a given \((l,m)\) pair
  integer, allocatable :: cg_lm(:,:,:,:)
  !> value non-zero Clebsch-Gordan coefficients for a given \((l,m)\) pair
  complex(dp), allocatable :: cg_val(:,:,:,:)

  ! INTERSTITIAL BASIS

  !> set of \({\bf G}\)-vectors for interstitial basis functions
  !> and real-space representation of operators
  type(G_set), pointer, public :: Gset
  !> set of \({\bf G}\)-vectors with twice the cut-off of `Gset`
  !> corresponding to real-space grid with twice the point density
  type(G_set) :: Gset2
  !> characteristic function on real space grid corresponding to `Gset2`
  real(dp), allocatable, public:: cfunr(:)

  public :: me_lapwlo_gen_helper_arrays, me_lapwlo_del_helper_arrays, &
            me_lapwlo_mt_rignt, me_lapwlo_ir_opig, &
            me_lapwlo_mt_mat, me_lapwlo_ir_mat

  contains

    ! ******************************************
    ! INITIALIZATION AND FINALIZATION

    !> Generate local helper arrays.
    subroutine me_lapwlo_gen_helper_arrays( mt_basis, lmax_operator, Gset0 )
      use constants, only: zzero
      use gaunt
      use mod_kpointset, only: generate_G_vectors
      !> (L)APW and LO radial functions
      type(mt_basis_type), intent(in) :: mt_basis
      !> maximum \(l\) used in operator expansion
      integer, intent(in) :: lmax_operator
      !> set of \({\bf G}\)-vectors for interstitial basis functions
      !> real-space representation of operators
      type(G_set), target, intent(in) :: Gset0

      integer :: lmax, l, m, lm, ig

      complex(dp), allocatable :: cfung(:), zfft(:)

      ! clean up
      call me_lapwlo_del_helper_arrays

      ! asign module variables
      lmaxb = mt_basis%lmax_basis
      lmaxo = lmax_operator
      basis = mt_basis

      lmax = max( lmaxb, lmaxo ) + 1 ! + 1 for gradient and surface terms
      
      ! build l, m <--> (l,m) maps
      allocate( lm_join(0:lmax, -lmax:lmax), source=0 )
      allocate( lm_split(2, (lmax+1)**2) )
      lm = 0
      do l = 0, lmax
        do m = -l, l
          lm = lm + 1
          lm_join(l, m) = lm
          lm_split(:, lm) = [l, m]
        end do
      end do

      ! check if Gaunt coefficients are available and create them if not
      if( .not. gaunt_coeff_yyy%check_bounds( lmax, lmaxo+1, lmaxb+1 ) ) &
        gaunt_coeff_yyy = non_zero_gaunt_yyy( lmax, lmaxo+1, lmaxb+1 ) 
      if( .not. gaunt_coeff_yry%check_bounds( lmax, lmaxo+1, lmaxb+1 ) ) &
        gaunt_coeff_yry = non_zero_gaunt_yry( lmax, lmaxo+1, lmaxb+1 ) 
      if( .not. gaunt_coeff_rrr%check_bounds( lmax, lmaxo+1, lmaxb+1 ) ) &
        gaunt_coeff_rrr = non_zero_gaunt_rrr( lmax, lmaxo+1, lmaxb+1 )

      ! generate Clebsch-Gordan coefficients
      call generate_non_zero_clebsch_gordan( lmaxb, 1e-64_dp, cg_num, cg_lm, cg_val )

      ! generate G-vectors sets
      Gset => Gset0
      call generate_G_vectors( Gset2, Gset%bvec, [[0,0,0],[0,0,0]], 2*Gset%gmaxvr, auto_intgv=.true. )
      
      ! generate characteristic function
      allocate( cfung(Gset2%ngrtot), cfunr(Gset2%ngrtot) )
      allocate( zfft(Gset2%ngrtot), source=zzero )
      call gencfunig( Gset2%ngvec, Gset2%gc, Gset2%vgc, cfung )
      do ig = 1, Gset2%ngvec
        zfft(Gset2%igfft(ig)) = cfung(ig)
      end do
      call zfftifc( 3, Gset2%ngrid, 1, zfft )
      cfunr = dble( zfft )
      deallocate( cfung, zfft )
    end subroutine me_lapwlo_gen_helper_arrays

    !> Delete local helper arrays.
    subroutine me_lapwlo_del_helper_arrays
      use mod_kpointset, only: delete_G_vectors
      if( allocated( lm_join ) ) deallocate( lm_join )
      if( allocated( lm_split ) ) deallocate( lm_split )
      if( allocated( cg_num ) ) deallocate( cg_num )
      if( allocated( cg_val ) ) deallocate( cg_val )
      if( allocated( cfunr ) ) deallocate( cfunr )
      if( associated( Gset ) ) nullify( Gset )
      call delete_G_vectors( Gset2 )
    end subroutine me_lapwlo_del_helper_arrays

    ! INITIALIZATION AND FINALIZATION
    ! ******************************************

    ! ******************************************
    ! MUFFIN-TIN INTEGRALS
    
    !> This subroutine computes the sum of the radial MT integrals between basis functions and an operator times Gaunt coefficients:
    !> \[
    !> O^{\alpha}_{\lambda_1 \lambda_2} = \int\limits_\alpha \phi^{\alpha\ast}_{\lambda_1}({\bf r})\,
    !> O({\bf r})\, \phi^\alpha_{\lambda_2}({\bf r})\, {\rm d}^3 r =
    !> \langle Y_{l_{\lambda_1} m_{\lambda_1}} | X_{lm} | Y_{l_{\lambda_2} m_{\lambda_2}} \rangle
    !> R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_{\lambda_1}, l_{\lambda_2}) \;,
    !> \]
    !> where \(R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_{\lambda_1}, l_{\lambda_2})\) are the 
    !> radial integrals of the opeartor radial functions \(o^\alpha_{lm}(r)\) as obtained from the subroutine
    !> [[me_lapwlo_mt_ri(subroutine)]].
    !> The operator can be given by an expansion in real (\(X_{lm} = R_{lm}\)) or complex (\(X_{lm} = Y_{lm}\))
    !> spherical harmonics with real or complex valued radial functions \(o^{\alpha}_{lm}(r)\).
    !> 
    !> Optionally, instead of \(g_{\tilde{\lambda}}(r)\), also its radial derivatives
    !> \(\frac{{\rm d}^k}{{\rm d}r^k}g_{\tilde{\lambda}}(r)\) can be used on both sides
    !> of the integral individually by the use of the arguments `left_radial_derivative`
    !> and `right_radial_derivative`.
    !>
    !> The integration domain \(\mathcal D\) is determined by the argument `surface_integral`
    !> and is either the intertior of the muffin-tin sphere or one of the Cartesian components
    !> of its surface normal element.
    !>
    !> Whether to use the basis functions \(\phi^\alpha_\lambda({\bf r})\) or one of the Cartesian components 
    !> of their gradients \({\bf \nabla}_i \phi^\alpha_\lambda({\bf r})\)
    !> is determined by the arguments `left_gradient` and `right_gradient`.
    !> For more details on gradients see [[muffin_tin_basis(module):get_gradient_rad_fun(function)]].
    !>
    !> This routine computes
    !> \[ O^{\alpha,{\rm out}}_{\lambda_1 \lambda_2} = 
    !>    a\, O^{\alpha}_{\lambda_1 \lambda_2} + b\, O^{\alpha,{\rm in}}_{\lambda_1 \lambda_2} \; . \]
    !>
    !> See also [[me_lapwlo_mt_gaunt_sum(subroutine)]].
    !
    !> @note
    !>
    !> *   More convenient interfaces for standard usecases are provided by the polymorphic
    !>     subroutine [[me_mt_prepare(subroutine)]].
    !>
    !> @endnote
    subroutine me_lapwlo_mt_rignt( is, ias, lmax_op, alpha, rfun, beta, rignt, &
        left_radial_derivative, right_radial_derivative, &
        left_gradient, right_gradient, surface_integral, real_expansion )
      !> index of the species of the MT \(\alpha\)
      integer, intent(in) :: is
      !> index of the atom of the MT \(\alpha\)
      integer, intent(in) :: ias
      !> maximum l for operator expansion
      integer, intent(in) :: lmax_op
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> radial functions of operator \(o^\alpha_{lm}(r)\)
      class(*), intent(in) :: rfun(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> On input: \(O^{\alpha,{\rm in}}_{\lambda_1 \lambda_2}\); 
      !> On output: \(O^{\alpha,{\rm out}}_{\lambda_1 \lambda_2}\) 
      !> (See also [[me_mt_alloc(subroutine)]].)
      complex(dp), intent(inout) :: rignt(:,:)
      !> order of radial derivative (0 or 1) for left basis functions
      integer, intent(in) :: left_radial_derivative
      !> order of radial derivative (0 or 1) for right basis functions
      integer, intent(in) :: right_radial_derivative
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions. 0, if no gradient should be used.
      integer, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions. 0, if no gradient should be used.
      integer, intent(in) :: right_gradient
      !> Cartesian component (1, 2, or 3) of surface normal for surface integrals. 0, if volume integral should be computed.
      integer, intent(in) :: surface_integral
      !> if `.true.`, operator is given as an real spherical harmonics expansion
      logical, intent(in) :: real_expansion

      integer :: lmax, ng1, ng2, lg1(2), lg2(2)
      integer :: llstart1, llstart2, llstop1, llstop2, ll1, ll2, ig1, ig2
      integer :: l1, l2, m1, m2, lm1, lm2, lam1, lam2, idx1, idx2

      complex(dp), allocatable :: zri(:,:,:)

      ! check prefactors
      if( beta == zzero ) then
        rignt = zzero
      else if( beta /= zone ) then
        rignt = beta * rignt
      end if
      if( alpha == zzero ) return

      ! sanity checks
      call assert( surface_integral >= 0 .and. surface_integral <= 3, &
        'Invalid value for `surface_integral` argument. Allowed values are 0, 1, 2, or 3.' )
      call assert( left_gradient >= 0 .and. left_gradient <= 3, &
        'Invalid value for `left_gradient`. Allowed values are 0, 1, 2, or 3.' )
      call assert( right_gradient >= 0 .and. right_gradient <= 3, &
        'Invalid value for `left_gradient`. Allowed values are 0, 1, 2, or 3.' )
      select type( rfun )
        type is( real(dp) )
        type is( complex(dp) )
        class default
          call assert( .false., &
            'Argument `rfun` must be of type `double real` or `double complex`.' )
      end select

      ! set maximum l
      lmax = lmax_op
      if( surface_integral > 0 ) lmax = lmax_op + 1

      ! set number of gradient terms
      ng1 = 1; lg1 = 0
      if( left_gradient > 0 ) then
        ng1 = 2
        lg1(1) = -1; lg1(2) = 1
      end if
      ng2 = 1; lg2 = 0
      if( right_gradient > 0 ) then
        ng2 = 2
        lg2(1) = -1; lg2(2) = 1
      end if

      allocate( zri(basis%n_rad_fun_max, basis%n_rad_fun_max, (lmax+1)**2) )

      ! loops over gradient contributions lg (lg = +1/-1 term for gradients, lg = 0 if no gradient)
      do ig1 = 1, ng1
        llstart1 = max( 0, lg1(ig1) )         ! limits for ll = l + lg
        llstop1 = lmaxb + lg1(ig1)
        do ig2 = 1, ng2
          llstart2 = max( 0, lg2(ig2) )
          llstop2 = lmaxb + lg2(ig2)

!$omp parallel default( shared ) private( ll1, ll2, l1, l2, m1, m2, lm1, lm2, lam1, lam2, idx1, idx2, zri )
!$omp do collapse(2)
          ! loops over ll = l+lg
          do ll1 = llstart1, llstop1
            do ll2 = llstart2, llstop2
              l1 = ll1 - lg1(ig1)
              l2 = ll2 - lg2(ig2)

              ! compute radial integrals
              call me_lapwlo_mt_ri( is, ias, l1, l2, lmax_op, rfun, zri, &
                     left_radial_derivative, right_radial_derivative, lg1(ig1), lg2(ig2), surface_integral, real_expansion )
              if( alpha /= zone ) zri = alpha * zri

              ! loop over m of basis functions
              do m1 = -l1, l1
                lm1 = lm_join(l1, m1)
                do m2 = -l2, l2
                  lm2 = lm_join(l2, m2)

                  ! loop over all basis functions with given lm
                  do lam1 = 1, basis%n_rad_fun(l1, is)
                    idx1 = basis%idx_basis_fun(lm1, lam1, is)
                    do lam2 = 1, basis%n_rad_fun(l2, is)
                      idx2 = basis%idx_basis_fun(lm2, lam2, is)

                      ! sum over Gaunt coefficients from basis functions and operator
                      call me_lapwlo_mt_gaunt_sum( lmax, lm1, lm2, lam1, lam2, lg1(ig1), lg2(ig2), left_gradient, right_gradient, real_expansion, zri, rignt(idx1, idx2) )

                    end do
                  end do

                end do
              end do

            end do
          end do
!$omp end do
!$omp end parallel

        end do
      end do
    end subroutine me_lapwlo_mt_rignt

    !> Computes radial MT integrals of basis functions with angular momentum
    !> \(l_1\) and \(l_2\), respectively, and a set or radial functions of an operator
    !> \[ R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2,lm}(l_1,l_2) = 
    !>    \int_0^{R_\alpha} g_{\tilde{\lambda}_1}(r) \, o^{\alpha}_{lm}(r) \, g_{\tilde{\lambda}_2}(r) \, r^2 \, {\rm d}r \;. \]
    !> The radial funtions of the operator can be either real or complex valued.
    !>
    !> When `left_gradient=-1` or `right_gradient=-1` instead of the bare basis functions
    !> \[ g^{\alpha,-}_{\tilde{\lambda}}(r) = \left[\frac{l}{2l+1}\right]^\frac{1}{2} \left[ \frac{l+1}{r} + \frac{{\rm d}}{{\rm d}r} \right] 
    !>    g^\alpha_{\tilde{\lambda}}(r) \]
    !> corresponding to the \({\bf Y}_{l\,l-1\,m}(\hat{\bf r})\) term in the equations in 
    !> [[muffin_tin_basis(module):get_gradient_rad_fun(function)]] is used.
    !>
    !> When `left_gradient=1` or `right_gradient=1` instead of the bare basis functions
    !> \[ g^{\alpha,+}_{\tilde{\lambda}}(r) = \left[\frac{l+1}{2l+1}\right]^\frac{1}{2} \left[ \frac{l}{r} - \frac{{\rm d}}{{\rm d}r} \right] 
    !>    g^\alpha_{\tilde{\lambda}}(r) \]
    !> corresponding to the \({\bf Y}_{l\,l+1\,m}(\hat{\bf r})\) term in the equations in
    !> [[muffin_tin_basis(module):get_gradient_rad_fun(function)]] is used.
    !>
    !> When `left_gradient=0` or `right_gradient=0` the bare basis functions \(g^\alpha_\tilde{\lambda}(r)\) are used.
    !> 
    !> For surface integrals, i.e., `surface_integral` = 1,2,3, no actual integral is sloved
    !> but the integrand is evaluated at the muffin-tin radius, i.e.,
    !> \[ R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2,lm}(l_1,l_2) =
    !>    g_{\tilde{\lambda}_1}(R_{\alpha}) \, o^{\alpha}_{lm}(R_{\alpha}) \, g_{\tilde{\lambda}_2}(R_{\alpha}) \, R_{\alpha}^2 \]
    !> will be computed. Further, additional Gaunt coefficients coming from the surface element
    !> will appear. See [[me_lapwlo_mt_ri_surf(subroutine)]] for further explanation.
    subroutine me_lapwlo_mt_ri( is, ias, l1, l2, lmax_op, rfun, ri, &
        left_radial_derivative, right_radial_derivative, &
        left_gradient, right_gradient, surface_integral, real_expansion )
      use constants, only: zzero
      !> index of the species of the MT \(\alpha\)
      integer, intent(in) :: is
      !> index of the atom of the MT \(\alpha\)
      integer, intent(in) :: ias
      !> angular momentum \(l_1\) of left basis functions
      integer, intent(in) :: l1
      !> angular momentum \(l_2\) of right basis functions
      integer, intent(in) :: l2
      !> maximum l for operator expansion
      integer, intent(in) :: lmax_op
      !> radial functions of operator \(o^\alpha_{lm}(r)\)
      class(*), intent(in) :: rfun(:,:)
      !> radial integrals \(R^\alpha_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_1, l_2)\)
      complex(dp), intent(out) :: ri(:,:,:)
      !> order of radial derivative (0 or 1) for left basis functions
      integer, intent(in) :: left_radial_derivative
      !> order of radial derivative (0 or 1) for right basis functions
      integer, intent(in) :: right_radial_derivative
      !> gradient component (-1, 0, or 1) for left basis functions
      integer, intent(in) :: left_gradient
      !> gradient component (-1, 0, or 1) for right basis functions
      integer, intent(in) :: right_gradient
      !> Cartesian component (1, 2, or 3) of surface normal for surface integrals
      integer, intent(in) :: surface_integral
      !> operator is given as an real spherical harmonics expansion
      logical, intent(in) :: real_expansion

      integer :: lmmax, ir, nr, l3, m3, lm3, lam1, lam2

      real(dp), allocatable :: rri_tmp(:,:,:)
      complex(dp), allocatable :: zri_tmp(:,:,:)
      real(dp), allocatable :: f1(:), f2(:), fr(:), gf(:), cf(:,:), r2(:)

      select type( rfun )
        type is( real(dp) )
        type is( complex(dp) )
        class default
          call assert( .false., &
            'Argument `rfun` must be of type `double real` or `double complex`.' )
      end select

      call assert( left_radial_derivative >= 0, &
        'Argument `left_radial_derivative` must not be negative.' )
      call assert( right_radial_derivative >= 0, &
        'Argument `right_radial_derivative` must not be negative.' )
      call assert( abs( left_gradient ) <= 1, &
        'Argument `left_gradient` must be either -1, 0, or 1.' )
      call assert( abs( right_gradient ) <= 1, &
        'Argument `right_gradient` must be either -1, 0, or 1.' )
      call assert( surface_integral >= 0 .and. surface_integral <= 3, &
        'Argument `surface_integral` must be either 0, 1, 2, or 3.' )

      lmmax = (lmax_op + 1)**2

      ri = zzero
      select type( rfun )
        type is( real(dp) )
          allocate( rri_tmp(size(ri,dim=1), size(ri,dim=2), size(ri,dim=3)), source=0.0_dp )
        type is( complex(dp) )
          allocate( zri_tmp(size(ri,dim=1), size(ri,dim=2), size(ri,dim=3)), source=zzero )
      end select

      nr = basis%n_rad_grid(is)
      allocate( r2(nr) )
      do ir = 1, nr
        r2(ir) = basis%rad_grid(ir, is)**2
      end do

      allocate( fr(nr), gf(nr), cf(3,nr) )

      do lam1 = 1, basis%n_rad_fun(l1, is)
        do lam2 = 1, basis%n_rad_fun(l2, is)

          f1 = basis%get_gradient_rad_fun( l1, is, ias, lam1, left_gradient, radial_derivative=left_radial_derivative )
          f2 = basis%get_gradient_rad_fun( l2, is, ias, lam2, right_gradient, radial_derivative=right_radial_derivative )

          do l3 = 0, lmax_op
            do m3 = -l3, l3
              lm3 = lm_join(l3, m3)

              if( allocated( rri_tmp ) ) then
                rri_tmp(lam1, lam2, lm3) = integrate_real( f1, f2 )
              else
                zri_tmp(lam1, lam2, lm3) = integrate_complex( f1, f2 )
              end if

            end do
          end do

        end do
      end do

      if( allocated( f1 ) ) deallocate( f1 )
      if( allocated( f2 ) ) deallocate( f2 )
      deallocate( r2, fr, gf, cf)

      if( allocated( rri_tmp ) ) then
        call me_lapwlo_mt_ri_surf( lmax_op, surface_integral, real_expansion, rri_tmp, ri )
      else
        call me_lapwlo_mt_ri_surf( lmax_op, surface_integral, real_expansion, zri_tmp, ri )
      end if

      contains
        function integrate_real( f1, f2 ) result( res )
          real(dp), intent(in) :: f1(:), f2(:)
          real(dp) :: res

          select type( rfun )
            type is( real(dp) )
              if( surface_integral == 0 ) then
                do ir = 1, nr
                  fr(ir) = f1(ir) * rfun(lm3, ir) * f2(ir) * r2(ir)
                end do
                call fderiv( -1, nr, basis%rad_grid(:, is), fr, gf, cf )
                res = gf(nr)
              else
                res = f1(nr) * rfun(lm3, nr) * f2(nr) * r2(nr)
              end if
            class default
              res = 0.0_dp
          end select
        end function integrate_real

        function integrate_complex( f1, f2 ) result( res )
          real(dp), intent(in) :: f1(:), f2(:)
          complex(dp) :: res

          real(dp) :: int1, int2

          select type( rfun )
            type is( complex(dp) )
              if( surface_integral == 0 ) then
                do ir = 1, nr
                  fr(ir) = f1(ir) * dble( rfun(lm3, ir) ) * f2(ir) * r2(ir)
                end do
                call fderiv( -1, nr, basis%rad_grid(:, is), fr, gf, cf )
                int1 = gf(nr)
                do ir = 1, nr
                  fr(ir) = f1(ir) * aimag( rfun(lm3, ir) ) * f2(ir) * r2(ir)
                end do
                call fderiv( -1, nr, basis%rad_grid(:, is), fr, gf, cf )
                int2 = gf(nr)
                res = cmplx( int1, int2, dp )
              else
                res = f1(nr) * rfun(lm3, nr) * f2(nr) * r2(nr)
              end if
            class default
              res = zzero
          end select
        end function integrate_complex
    end subroutine me_lapwlo_mt_ri

    !> Account for additional Gaunt coefficients in the case of integrals over the muffin-tin surface.
    !>
    !> When calculating integrals over the muffin-tin surface, \(\partial \alpha\), the sum of the radial
    !> integrals times the Gaunt coefficients (see [[me_lapwlo_mt_rignt(subroutine)]]) takes the form
    !> \[ \begin{align*} O^{\alpha}_{\lambda_1 \lambda_2} 
    !>    &= \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^l \sum_{m'=-1}^1 C_{im'} 
    !>    \langle Y_{l_{\lambda_1}m_{\lambda_2}} | X_{lm} | X_{1m'} | Y_{l_{\lambda_2}m_{\lambda_2}} \rangle 
    !>    R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2,lm}(l_{\lambda_1},l_{\lambda_2}) \\
    !>    &= \sum_{l''=0}^{l_{\rm max,o}+1} \sum_{m''=-l''}^{l''}
    !>    \langle Y_{l_{\lambda_1}m_{\lambda_2}} | X_{l''m''} | Y_{l_{\lambda_2}m_{\lambda_2}} \rangle 
    !>    R^{\partial_i \alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2,l''m''}(l_{\lambda_1},l_{\lambda_2}) \;, 
    !>    \end{align*} \]
    !> where the coefficients \(C_{im}\) describe the surface normal in terms of spherical harmonics, 
    !> \[ \hat{\bf e}(\theta,\phi) 
    !>     = \begin{pmatrix} \sin\theta \cos\phi \\ \sin\theta \sin\phi \\ \cos\theta \end{pmatrix}
    !>     = \sum_{m=-1}^1 C_{:m} X_{1m}(\theta,\phi) \;, \]
    !> with 
    !> \[ {\bf C} = \sqrt{\frac{4\pi}{3}} \begin{pmatrix} 1/\sqrt{2} & 0 & -1/\sqrt{2} \\
    !>    {\rm i}/\sqrt{2} & 0 & {\rm i}/\sqrt{2} \\ 0 & 1 & 0 \end{pmatrix} 
    !>    \text{ for } X_{lm} = Y_{lm}\; \text{ and} \; 
    !>    {\bf C} = \sqrt{\frac{4\pi}{3}} \begin{pmatrix} 0 & 0 & -1 \\ -1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} 
    !>    \text{ for } X_{lm} = R_{lm} \;. \]
    !> The product of spherical harmonics \(X_{lm} \, X_{1m'}\) in the above expression is expandend
    !> into a series of single spherical harmonics \(X_{l''m''}\) using additional Gaunt coefficients
    !> in order to recast the summation structure as defined in [[me_lapwlo_mt_rignt(subroutine)]]. 
    !> Hence, this subroutine computes
    !> \[ R^{\partial_i \alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2,l''m''}(l_{\lambda_1},l_{\lambda_2}) 
    !>    = \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^{l} 
    !>    R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2,lm}(l_{\lambda_1},l_{\lambda_2}) 
    !>    \sum_{m'=-1}^1 C_{im'} \langle X_{l''m''} | X_{lm} | X_{1m'} \rangle \;. \]
    subroutine me_lapwlo_mt_ri_surf( lmax_op, surface_integral, real_expansion, ri, ri_surf )
      use constants, only: fourpi, zzero, zone, zi, sqrt_two
      use gaunt
      !> maximum l for operator expansion
      integer, intent(in) :: lmax_op
      !> Cartesian component (1, 2, or 3) of surface normal for surface integrals
      integer, intent(in) :: surface_integral
      !> operator is given as an real spherical harmonics expansion
      logical, intent(in) :: real_expansion
      !> radial integrals \(R^\alpha_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}\)
      class(*), intent(in) :: ri(:,:,:)
      !> surface radial integrals \(R^{\partial_i \alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2, l''m''}\)
      complex(dp), intent(out) :: ri_surf(:,:,:)

      integer :: lmmax, nm, mm(2), i, j, l, m, lm, lm0, lm2
      real(dp) :: cr(2), r1
      complex(dp) :: cz(2), z1
      type(non_zero_gaunt_real), pointer :: gnt

      select type( ri )
        type is( real(dp) )
        type is( complex(dp) )
        class default
          call assert( .false., &
            'Argument `ri` must be of type `double real` or `double complex`.' )
      end select

      call assert( surface_integral >= 0 .and. surface_integral <= 3, &
        'Argument `surface_integral` must be either 0, 1, 2, or 3.' )

      lmmax = (lmax_op + 1)**2

      if( surface_integral == 0 ) then
        select type( ri )
          type is( real(dp) )
            ri_surf = cmplx( ri, 0.0_dp, dp )
          type is( complex(dp) )
            ri_surf = ri
        end select
        return
      else
        ri_surf = zzero
      end if
      
      if( real_expansion ) then
        cr = 0.0_dp
        select case( surface_integral )
          case(1)
            nm = 1
            mm(1) = 1
            cr(1) = -1.0_dp
          case(2)
            nm = 1
            mm(1) = -1
            cr(1) = -1.0_dp
          case(3)
            nm = 1
            mm(1) = 0
            cr(1) = 1.0_dp
        end select
        cr = cr * sqrt( fourpi / 3.0_dp )
      else
        cz = zzero
        select case( surface_integral )
          case(1)
            nm = 2
            mm(1) = -1; mm(2) = 1
            cz(1) = zone / sqrt_two; cz(2) = -zone / sqrt_two
          case(2)
            nm = 2
            mm(1) = -1; mm(2) = 1
            cz(1) = zi / sqrt_two; cz(2) = zi / sqrt_two
          case(3)
            nm = 1
            mm(1) = 0
            cz(1) = zone
        end select
        cz = cz * sqrt( fourpi / 3._dp )
      end if

      do i = 1, nm
        lm0 = lm_join(1, mm(i))
        
        do l = 0, lmax_op+1
          do m = -l, l
            lm = lm_join(l ,m)
            if( real_expansion ) then
              gnt => gaunt_coeff_rrr
            else
              gnt => gaunt_coeff_yyy
            end if

            do j = 1, gnt%num(lm, lm0)
              lm2 = gnt%lm2(j, lm, lm0)
              if( lm2 > lmmax ) exit
              if( real_expansion ) then
                r1 = cr(i) * gnt%val(j, lm, lm0)
                select type( ri )
                  type is( real(dp) )
                    ri_surf(:, :, lm) = ri_surf(:, :, lm) + r1 * ri(:, :, lm2)
                  type is( complex(dp) )
                    ri_surf(:, :, lm) = ri_surf(:, :, lm) + r1 * ri(:, :, lm2)
                end select
              else
                z1 = cz(i) * gnt%val(j, lm, lm0)
                select type( ri )
                  type is( real(dp) )
                    ri_surf(:, :, lm) = ri_surf(:, :, lm) + z1 * ri(:, :, lm2)
                  type is( complex(dp) )
                    ri_surf(:, :, lm) = ri_surf(:, :, lm) + z1 * ri(:, :, lm2)
                end select
              end if
            end do

          end do
        end do

      end do
    end subroutine me_lapwlo_mt_ri_surf

    !> This subroutine compute the sum of radial integrals times Gaunt coefficients, i.e., 
    !> \[
    !> O^{\alpha}_{\lambda_1 \lambda_2} = \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^l 
    !> \langle Y_{l_{\lambda_1} m_{\lambda_1}} | X_{lm} | Y_{l_{\lambda_2} m_{\lambda_2}} \rangle
    !> R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_{\lambda_1}, l_{\lambda_2}) \;.
    !> \]
    !> 
    !> When integrals involving the gradient of basis functions have to be computed, additional sums
    !> of the type \(\sum_{m'=-(l\pm 1)}^{l\pm 1}\) and coefficients \(c^{\pm}_{i,mm',l}\) occur,
    !> where \(\pm\) corresponds to the \(l+1\) or \(l-1\) term in the gradient expression 
    !> given in [[muffin_tin_basis(module):get_gradient_rad_fun(function)]] and is 
    !> given by the input parameters `g1` and `g2`, respectively, and \(i\) is the Cartesian direction
    !> of the gradient, given by `left_gradient` and `right_gradient`, respectively.
    subroutine me_lapwlo_mt_gaunt_sum( lmax_op, lm1, lm2, lam1, lam2, g1, g2, left_gradient, right_gradient, real_expansion, ri, rignt )
      use constants, only: sqrt_two,zzero
      use wigner3j_symbol, only: clebsch_gordan
      use gaunt
      !> maximum \(l\) for operator expansion
      integer, intent(in) :: lmax_op
      !> \((l,m)\) on the left and right
      integer, intent(in) :: lm1, lm2
      !> \(\tilde{\lambda}\) on the left and right
      integer, intent(in) :: lam1, lam2
      !> part of left radial functions in case of gradients (-1/+1 for gradient, 0 if no gradient)
      integer, intent(in) :: g1
      !> part of right radial functions in case of gradients (-1/+1 for gradient, 0 if no gradient)
      integer, intent(in) :: g2
      !> gradient direction on the left (1,2,3 for gradient, 0 if no gradient)
      integer, intent(in) :: left_gradient
      !> gradient direction on the right (1,2,3 for gradient, 0 if no gradient)
      integer, intent(in) :: right_gradient
      !> `.true.` if operator is expanded in real spherical harmonics
      logical, intent(in) :: real_expansion
      !> radial integrals \(R^\alpha_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_1, l_2)\)
      complex(dp), intent(in) :: ri(:,:,:)
      !> radial integrals times gaunt coefficients \(O^\alpha_{\lambda_1 \lambda_2}\)
      complex(dp), intent(inout) :: rignt

      integer :: i, j, k, lm, lmmax, lmm1, lmm2
      complex(dp) :: cg1, cg2, z1
      type(non_zero_gaunt_real), pointer :: gntr
      type(non_zero_gaunt_complex), pointer :: gntz

      lmmax = (lmax_op + 1)**2

      do i = 1, cg_num(left_gradient, g1, lm1)
        lmm1 = cg_lm(i, left_gradient, g1, lm1)
        cg1 = cg_val(i, left_gradient, g1, lm1)
        do j = 1, cg_num(right_gradient, g2, lm2)
          lmm2 = cg_lm(j, right_gradient, g2, lm2)
          cg2 = cg_val(j, right_gradient, g2, lm2)

          if( real_expansion ) then
            gntz => gaunt_coeff_yry
            do k = 1, gntz%num(lmm1, lmm2)
              lm = gntz%lm2(k, lmm1, lmm2)
              if( lm > lmmax ) exit
              z1 = conjg( cg1 ) * cg2 * gntz%val(k, lmm1, lmm2)
              rignt = rignt + z1 * ri(lam1, lam2, lm)
            end do
          else
            gntr => gaunt_coeff_yyy
            do k = 1, gntr%num(lmm1, lmm2)
              lm = gntr%lm2(k, lmm1, lmm2)
              if( lm > lmmax ) exit
              z1 = conjg( cg1 ) * cg2 * gntr%val(k, lmm1, lmm2)
              rignt = rignt + z1 * ri(lam1, lam2, lm)
            end do
          end if

        end do
      end do
    end subroutine me_lapwlo_mt_gaunt_sum

    !> Having prepared the \({\bf p}\) independent radial integrals times Gaunt coefficients
    !> \(O^\alpha_{\lambda_1 \lambda_2}\) (`rignt`) using [[me_lapwlo_mt_rignt(subroutine)]], this subroutine 
    !> computes the contribution to the matrix elements coming from the given MT \(\alpha\)
    !> \[ M^\alpha_{\mu \nu} = \langle \chi_\mu | O | \chi_\nu \rangle_\mathcal{D} \; . \]
    !>
    !> If the optional arguments `left_evec` or `right_evec` are given, 
    !> the corresponding \(\chi({\bf r})\) are taken to be wavefunctions expressed as linear combinations
    !> of (L)APW+LO basis functions, where the eigenvectors give the corresponding expansion coefficients.
    !> In this case, the matrix elements are given by
    !> \[ M^\alpha_{mn} = \langle \psi_{m{\bf p}_1} | O | \psi_{n{\bf p}_2} \rangle_\mathcal{D} = 
    !>    \sum\limits_{\lambda_1,\lambda_2} {C^{m{\bf p}_1\alpha}_{\lambda_1}}^\ast\, O^\alpha_{\lambda_1 \lambda_2}\, 
    !>    C^{n{\bf p}_2\alpha}_{\lambda_2} \; , \]
    !> where \(C^{n{\bf p}\alpha}_\lambda\) are the muffin-tin eigenvectors as obtained from 
    !> [[muffin_tin_basis(module):transform_evec(subroutine)]].
    !>
    !> Otherwise, \(\chi({\bf r})\) is taken to be (L)APW+LO basis functions and the corresponding
    !> matrix elements are given by
    !> \[ M^\alpha_{\mu \nu} = \langle \phi_{\mu{\bf p}_1} | O | \phi_{\nu{\bf p}_2} \rangle_\mathcal{D} = 
    !>    \sum\limits_{\lambda_1,\lambda_2} {T^\alpha_{\lambda_1 \mu}}^\ast\, O^\alpha_{\lambda_1 \lambda_2}\, 
    !>    T^\alpha_{\lambda_2 \nu} \; , \]
    !> where \(T^\alpha_{\lambda \mu}\) is the basis transformation matrix as obtained from
    !> [[muffin_tin_basis(module):get_basis_transform(subroutine)]].
    !>
    !> Use the optional arguments `left_local_orbitals` and `right_local_orbitals` to include or exclude LOs from the calculation.
    !>
    !> This routine computes
    !> \[ M^{\alpha,{\rm out}}_{\mu \nu} = a\, M^\alpha_{\mu \nu} + b\, M^{\alpha,{\rm in}}_{\mu \nu} \; . \]
    !>
    !
    !> @note
    !>
    !> *   More convenient interfaces for standard usecases are provided by the polymorphic
    !>     subroutine [[me_mt_mat(subroutine)]].
    !> *   The selection whether \(\chi({\bf r})\) or its spatial derivative \(\partial_i \chi({\bf r})\)
    !>     is used or whether the integration domain \(\mathcal{D}\) is the interior or the surface of the MT
    !>     had already been made in the calculation of `rignt`! See [[me_lapwlo_mt_rignt(subroutine)]].
    !>
    !> @endnote
    subroutine me_lapwlo_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, rignt, beta, mat, &
        left_evec, right_evec, diagonal_only, left_local_orbitals, right_local_orbitals )
      use mod_APW_LO, only: nlotot
      use mod_eigensystem, only: idxlo
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> number of \({\bf G+p}\)-vectors at left \({\bf p}\)-point
      integer, intent(in) :: ngp1
      !> number of \({\bf G+p}\)-vectors at right \({\bf p}\)-point
      integer, intent(in) :: ngp2
      !> (L)APW matching coefficients at left \({\bf p}\)-point
      complex(dp), intent(in) :: apwalm1(:,:,:)
      !> (L)APW matching coefficients at right \({\bf p}\)-point
      complex(dp), intent(in) :: apwalm2(:,:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> radial integrals times Gaunt coefficients \(O^\alpha_{\lambda_1 \lambda_2}\)
      class(*), intent(in) :: rignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements \(M^\alpha_{\mu \nu}\)
      complex(dp), intent(inout) :: mat(:,:)
      !> eigenvectors on the left
      complex(dp), optional, intent(in) :: left_evec(:,:)
      !> eigenvectors on the right
      complex(dp), optional, intent(in) :: right_evec(:,:)
      !> compute only diagonal_onlyonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2
      integer :: dimm(2), dimv1(2), dimv2(2)
      integer :: i, n1, n2
      complex(dp), allocatable :: evecmt(:,:), auxmat(:,:)
      complex(dp), external :: zdotc

      ! check prefactors
      if( beta == zzero ) then
        mat = zzero
      else if( beta /= zone ) then
        mat = beta * mat
      end if
      if( alpha == zzero ) return

      select type( rignt )
        type is( real(dp) )
        type is( complex(dp) )
        class default
          call assert( .false., &
            'Argument `rignt` must be of type `double real` or `double complex`.' )
      end select

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      dimm = shape( mat )
      if( present( left_evec ) ) dimv1 = shape( left_evec )
      if( present( right_evec ) ) dimv2 = shape( right_evec )

      n1 = ngp1 + nlotot
      if( present( left_evec ) ) then
        call assert( dimv1(1) >= n1, &
          'First dimension of left eigenvector too small.' )
        n1 = dimv1(2)
      end if
      n2 = ngp2 + nlotot
      if( present( right_evec ) ) then
        call assert( dimv2(1) >= n2, &
          'First dimension of right eigenvector too small.' )
        n2 = dimv2(2)
      end if
      call assert( dimm(1) >= n1 .and. (diag .or. dimm(2) >= n2), &
        'Dimensions of output matrix too small.' )

      allocate( auxmat(basis%n_basis_fun(is), n2) )

      if( present( right_evec ) ) then
        call basis%transform_evec( is, ngp2, nlotot, idxlo(:, :, ias), apwalm2, right_evec, evecmt, &
          use_local_orbitals=lo2 )
      else
        call basis%get_basis_transform( is, ngp2, nlotot, idxlo(:, :, ias), apwalm2, evecmt, &
          use_local_orbitals=lo2 )
      end if

      select type( rignt )
        type is( real(dp) )
          do i = 1, n2
            call dgemv( 'n', basis%n_basis_fun(is), basis%n_basis_fun(is), 1.0_dp, &
              rignt, size(rignt, dim=1), evecmt(1, i), 2, 0.0_dp, auxmat(1, i), 2 )
            call dgemv( 'n', basis%n_basis_fun(is), basis%n_basis_fun(is), 1.0_dp, &
              rignt, size(rignt, dim=1), evecmt(2, i), 2, 0.0_dp, auxmat(2, i), 2 )
          end do
        type is( complex(dp) )
          call zgemm( 'n', 'n', basis%n_basis_fun(is), n2, basis%n_basis_fun(is), zone, &
                 rignt, size(rignt, dim=1), &
                 evecmt, size(evecmt, dim=1), zzero, &
                 auxmat, basis%n_basis_fun(is) )
      end select
      if( alpha /= zone ) auxmat = alpha * auxmat

      if( present( left_evec ) ) then
        call basis%transform_evec( is, ngp1, nlotot, idxlo(:, :, ias), apwalm1, left_evec, evecmt, &
          use_local_orbitals=lo1 )
      else
        call basis%get_basis_transform( is, ngp1, nlotot, idxlo(:, :, ias), apwalm1, evecmt, &
          use_local_orbitals=lo1 )
      end if

      if( diag ) then
        do i = 1, min( n1, n2 )
          mat(i, 1) = mat(i, 1) + zdotc( basis%n_basis_fun(is), evecmt(1, i), 1, auxmat(1, i), 1 )
        end do
      else
        call zgemm( 'c', 'n', n1, n2, basis%n_basis_fun(is), zone, &
               evecmt, size(evecmt, dim=1), &
               auxmat, basis%n_basis_fun(is), zone, &
               mat, dimm(1) )
      end if

      if( allocated( evecmt ) ) deallocate( evecmt )
      deallocate( auxmat )
    end subroutine me_lapwlo_mt_mat

    ! END MUFFIN-TIN INTEGRALS
    ! ******************************************

    ! ******************************************
    ! INTERSTITIAL INTEGRALS
    
    !> This subroutine computes the reciprocal space representation \(\hat{O}_f({\bf G})\) of the product 
    !> of an operator \(O({\bf r})\) and a second interstitial functions \(f({\bf r})\)
    !>
    !> \[
    !> O({\bf r}) \, f({\bf r}) = 
    !> \left( \sum_{\bf G} \hat{O}({\bf G}) \, {\rm e}^{{\rm i} {\bf G} \cdot {\bf r}} \right)
    !> \left( \sum_{\bf G'} \hat{f}({\bf G'}) \, {\rm e}^{{\rm i} {\bf G'} \cdot {\bf r}} \right)
    !> = \sum_{\bf G} \left( \sum_{\bf G'} \hat{O}({\bf G'}) \hat{f}({\bf G-G'}) \right) {\rm e}^{{\rm i} {\bf G} \cdot {\bf r}} 
    !> = \sum_{\bf G} \hat{O}_f({\bf G})\, {\rm e}^{{\rm i}{\bf G} \cdot {\bf r}} \; , \]
    !> with
    !> \[ \hat{O}_f({\bf G}) := \sum_{\bf G'} \hat{O}({\bf G'}) \hat{f}({\bf G-G'}) \; . \]
    !> 
    !> \(\hat{O}({\bf G})\) will be given on the \({\bf G}\)-vector set given by the module variable `Gset`. 
    !> In practice, this is computed by first interpolating \(O({\bf r})\) on a FFT grid with twice the point density. 
    !> Then the product \(O({\bf r}) \, f({\bf r})\) is evaluated in real space and transformed back to reciprocal space
    !> and truncated to the coarse grid `Gset`.
    !>
    !> This routine computes
    !> \[ \hat{O}_f^{\rm out}({\bf G}) = a\, \hat{O}_f({\bf G}) + b\, \hat{O}_f^{\rm in}({\bf G}) \; . \]
    subroutine me_lapwlo_ir_opig( alpha, opir, Gset_op, fir, beta, opig, Gset_prod )
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> set of \({\bf G}\)-vectors used in the Fourier series of the operator \(O({\bf r}\)
      !> and defining the real-space FFT grid, `opir` is given on
      type(G_set), intent(in) :: Gset_op
      !> operator \(O({\bf r})\) in real space 
      !> (on FFT grid corresponding to input argument `Gset_op`)
      class(*), intent(in) :: opir(:)
      !> second function \(f({\bf r})\) in real space 
      !> (on doubled FFT grid corresponding to module variable `Gset2`)
      class(*), intent(in) :: fir(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> set of \({\bf G}\)-vectors used in the Fourier series of the product \(O({\bf r} \, f({\bf r})\)
      type(G_set), intent(in) :: Gset_prod
      !> reciprocal space representation of the product \(\hat{O}_f({\bf G})\)
      complex(dp), intent(inout) :: opig(:)
      
      integer :: ig, ig2, ivg(3)

      complex(dp), allocatable :: zfft(:), zfft2(:)

      ! check array sizes
      call assert( size( opir ) == Gset_op%ngrtot, '(me_lapwlo_ir_opig) &
        Array `opir` must have size `Gset_op%ngrtot`.' )
      call assert( size( fir ) == Gset2%ngrtot, '(me_lapwlo_ir_opig) &
        Array `fir` must have size `Gset2%ngrtot`.' )
      call assert( size( opig ) == Gset_prod%ngvec, '(me_lapwlo_ir_opig) &
        Array `opig` must have size `Gset_prod%ngvec`.' )

      ! check prefactors
      if( beta == zzero ) then
        opig = zzero
      else if( beta /= zone ) then
        opig = beta * opig
      end if
      if( alpha == zzero ) return

      select type( opir )
        type is( real(dp) )
        type is( complex(dp) )
        class default
          call assert( .false., &
            'Argument `opir` must be of type `double real` or `double complex`.' )
      end select

      select type( fir )
        type is( real(dp) )
        type is( complex(dp) )
        class default
          call assert( .false., &
            'Argument `fir` must be of type `double real` or `double complex`.' )
      end select

      if( 2*Gset_op%gmaxvr > Gset2%gmaxvr + 1e-16_dp ) then
        if( mpiglobal%rank == 0 ) then
          write( *, * )
          write( *, '("Warning (me_lapwlo_ir_opig): Maximum G for operator expansion exceeds maximum G of module G-vector set.")' )
          write( *, '("This might lead to a loss of precision.")' )
        end if
      end if

      ! interpolate operator onto doubled FFT grid
      allocate( zfft(Gset_op%ngrtot) )
      allocate( zfft2(Gset2%ngrtot) )
      select type( opir )
        type is( real(dp) )
          zfft = cmplx( opir, 0.0_dp, dp )
        type is( complex(dp) )
          zfft = opir
      end select
      if( alpha /= zone ) zfft = alpha * zfft
      call me_lapwlo_interpolate_ir( Gset_op, zfft, Gset2, zfft2 )

      ! multiply with second function on doubled FFT grid
      select type( fir )
        type is( real(dp) )
          zfft2 = zfft2 * cmplx( fir, 0.0_dp, dp )
        type is( complex(dp) )
          zfft2 = zfft2 * fir
      end select

      ! transform to reciprocal space
      call zfftifc( 3, Gset2%ngrid, -1, zfft2 )

      ! map FFT points to G-grid `Gset`
!$omp parallel default( shared ) private( ig, ig2, ivg )
!$omp do
      do ig = 1, Gset_prod%ngvec
        ivg = modulo( Gset_prod%ivg(:, ig) - Gset2%intgv(:, 1), Gset2%ngrid ) + Gset2%intgv(:, 1)
        ig2 = Gset2%ivgig(ivg(1), ivg(2), ivg(3))
        opig(ig) = opig(ig) + zfft2(Gset2%igfft(ig2))
      end do
!$omp end do
!$omp end parallel
    end subroutine me_lapwlo_ir_opig

    !> Interpolate a function given on the real space FFT grid corresponding to `Gset1`
    !> to the real space FFT grid corresponding to `Gset2`.
    subroutine me_lapwlo_interpolate_ir( Gset1, zfir1, Gset2, zfir2 )
      !> set of \({\bf G}\)-vectors the input function is defined on
      type(G_set), intent(in) :: Gset1
      !> on input: interstitial function on real-space FFT grid corresponding to `Gset1`; 
      !> on output: its Fourier components on reciprocal-space FFT grid
      complex(dp), intent(inout) :: zfir1(Gset1%ngrtot)
      !> set of \({\bf G}\)-vectors the output function is defined on
      type(G_set), intent(in) :: Gset2
      !> interstitial function on real-space FFT grid corresponding to `Gset2`
      complex(dp), intent(out) :: zfir2(Gset2%ngrtot)

      integer :: ig1, ig2, ivg(3)

      ! transform function to reciprocal space
      call zfftifc( 3, Gset1%ngrid, -1, zfir1 )

      ! copy Fourier components to dense grid
      zfir2 = zzero
!$omp parallel default( shared ) private( ig1, ig2, ivg )
!$omp do
      do ig1 = 1, Gset1%ngvec
        ivg = modulo( Gset1%ivg(:, ig1) - Gset2%intgv(:, 1), Gset2%ngrid ) + Gset2%intgv(:, 1)
        ig2 = Gset2%ivgig(ivg(1), ivg(2), ivg(3))
        zfir2(Gset2%igfft(ig2)) = zfir1(Gset1%igfft(ig1))
      end do
!$omp end do
!$omp end parallel

      ! transform function to real space
      call zfftifc( 3, Gset2%ngrid, 1, zfir2 )
    end subroutine me_lapwlo_interpolate_ir

    !> Having prepared the \({\bf p}\) independent interstitial representation of the operator
    !> \(\hat{O}_\Theta({\bf G})\) (`opig`) using [[me_lapwlo_ir_opig(subroutine)]], this subroutine 
    !> computes the contribution to the matrix elements coming from the interstitial region
    !> \[ M^{\rm IR}_{\mu \nu} = \langle \chi_\mu | O | \chi_\nu \rangle_{\rm IR} \; . \]
    !>
    !> If no eigenvectors are given, \(\chi({\bf r})\) is taken to be plane wave basis functions 
    !> \(\phi_{{\bf G}+{\bf p}}({\bf r})\) and the matrix elements are given by
    !> \[ \begin{align*}
    !> M^{\rm IR}_{{\bf G}_1 {\bf G}_2} 
    !> &= \int_{\rm IR} \phi_{{\bf G}_1 + {\bf p}_1}^\ast({\bf r}) \, O({\bf r}) \, \phi_{{\bf G}_2 + {\bf p}_2}({\bf r}) \, {\rm d}^3 r \\
    !> &= \sum_{\bf G} \hat{O}({\bf G}) \, \frac{1}{\Omega} \int_{\Omega}
    !> {\rm e}^{-{\rm i} ({\bf G}_1 + {\bf p}_1) \cdot {\bf r}} \,
    !> {\rm e}^{{\rm i} ({\bf G} + {\bf p}_1 - {\bf p}_2) \cdot {\bf r}} \,
    !> {\rm e}^{{\rm i} ({\bf G}_2 + {\bf p}_2) \cdot {\bf r}} \, \Theta_{\rm IR}({\bf r}) \, {\rm d}^3 r \\
    !> &= \sum_{\bf G} \hat{O}({\bf G}) \, \hat{\Theta}_{\rm IR}({\bf G}_1 - {\bf G}_2 - {\bf G}) \\
    !> &= \hat{O}_\Theta({\bf G}_1 - {\bf G}_2) \;.
    !> \end{align*}\]
    !>
    !> If the optional arguments `left_evec` or `right_evec` are given, 
    !> the corresponding \(\chi({\bf r})\) are taken to be wavefunctions expressed as linear combinations
    !> of plane wave basis functions, where the eigenvectors give the corresponding expansion coefficients.
    !> In this case, the matrix elements are given by
    !> \[ M^{\rm IR}_{mn} = \sum_{{\bf G}_1, {\bf G}_2} {C^{m{\bf p}_1}_{{\bf G}_1+{\bf p}_1}}^\ast \, 
    !> M^{\rm IR}_{{\bf G}_1 {\bf G}_2} \, C^{n{\bf p}_2}_{{\bf G}_2+{\bf p}_2} \; . \]
    !>
    !> Use the optional arguments `left_gradient` and `right_gradient` to use gradient components
    !> \({\bf \nabla}_i \chi({\bf r})\) instead.
    !>
    !> Use the optional argument `G_shift` to use \(\hat{O}_\Theta({\bf G}_1 - {\bf G}_2 + {\bf G}_0)\) instead,
    !> which is equivalent to a multiplication of the operator with \({\rm e}^{{\rm i}{\bf G}_0\cdot{\bf r}}\).
    !>
    !> Use the optional arguments `non_reduced_p1` and `non_reduced_p2` to indicate that `ip1` and `ip2`
    !> refer to non-reduced points within `Gpset1` and `Gpset2`.
    !>
    !> This routine computes
    !> \[ M^{{\rm IR},{\rm out}}_{\mu \nu} = a\, M^{\rm IR}_{\mu \nu} + b\, M^{{\rm IR},{\rm in}}_{\mu \nu} \; . \]
    !
    !> @note 
    !>
    !> *   More convenient interfaces for standard usecases are provided by the polymorphic
    !>     subroutine [[me_ir_mat(subroutine)]].
    !> *   For physical reasons, the operator \(O({\bf r})\) has to be a Bloch wave with wavevector
    !>     \({p_1 - p_2}\). This has to be ensured by the user.
    !> *   In general, also a different function than the interstitial characteristic function \(\Theta_{\rm IR}({\bf r})\)
    !>     can be used. See therefore [[me_lapwlo_ir_opig(subroutine)]].
    !>
    !> @endnote
    subroutine me_lapwlo_ir_mat( Gpset1, ip1, Gpset2, ip2, Gset_op, alpha, opig, beta, mat, &
        left_evec, right_evec, diagonal_only, G_shift, &
        left_gradient, right_gradient, non_reduced_p1, non_reduced_p2 )
      use mod_kpointset, only: Gk_set
      !> set of \({\bf G+p}\) vectors for basis functions on the left
      type(Gk_set), intent(in) :: Gpset1
      !> index of \({\bf p}\)-point in the left
      integer, intent(in) :: ip1
      !> set of \({\bf G+p}\) vectors for basis functions on the right
      type(Gk_set), intent(in) :: Gpset2
      !> index of \({\bf p}\)-point in the right
      integer, intent(in) :: ip2
      !> set of \({\bf G}\)-vectors input argument `opig` is defined on
      !> (`=Gset` if `opig` was created using [[me_lapwlo_ir_opig(subroutine)]])
      type(G_set), intent(in) :: Gset_op
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> operator (times characteristic function) \(\hat{O}_\Theta({\bf G})\) in reciprocal space 
      !> (see [[me_lapwlo_ir_opig(subroutine)]])
      complex(dp), intent(in) :: opig(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements \(M^\alpha_{\mu \nu}\)
      complex(dp), intent(inout) :: mat(:,:)
      !> eigenvectors on the left
      complex(dp), optional, intent(in) :: left_evec(:,:)
      !> eigenvectors on the right
      complex(dp), optional, intent(in) :: right_evec(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> constant shift lattice vector \({\bf G}_0\) to be applied to the operator or, equivalently,
      !> multiplication of operator with \({\rm e}^{{\rm i}{\bf G}_0\cdot{\bf r}}\) (default: `[0,0,0]`)
      integer, optional, intent(in) :: G_shift(3)
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient
      !> `ip1` is the index of a non-reduced \({\bf p}\)-point (default: `.false.`)
      logical, optional, intent(in) :: non_reduced_p1
      !> `ip2` is the index of a non-reduced \({\bf p}\)-point (default: `.false.`)
      logical, optional, intent(in) :: non_reduced_p2

      integer :: dimm(2), dimv1(2), dimv2(2), lgrad, rgrad
      integer :: igp1, igp2, ivg0(3), ivg(3), ig, i, n1, n2, ngp1, ngp2
      logical :: diag, nr1, nr2

      integer, allocatable :: igpig1(:), igpig2(:)
      real(dp), allocatable :: vgpc1(:,:), vgpc2(:,:)
      complex(dp), allocatable :: auxmat1(:,:), auxmat2(:,:)

      complex(dp), external :: zdotc

      ! check prefactors
      if( beta == zzero ) then
        mat = zzero
      else if( beta /= zone ) then
        mat = beta * mat
      end if
      if( alpha == zzero ) return

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lgrad = 0; if( present( left_gradient ) ) lgrad = left_gradient
      rgrad = 0; if( present( right_gradient ) ) rgrad = right_gradient
      ivg0 = 0; if( present( G_shift ) ) ivg0 = G_shift
      nr1 = .false.; if( present( non_reduced_p1 ) ) nr1 = non_reduced_p1
      nr2 = .false.; if( present( non_reduced_p2 ) ) nr2 = non_reduced_p2

      ngp1 = Gpset1%ngk(1, ip1)
      if( nr1 ) ngp1 = Gpset1%ngknr(1, ip1)
      ngp2 = Gpset2%ngk(1, ip2)
      if( nr2 ) ngp2 = Gpset2%ngknr(1, ip2)
      dimm = shape( mat )
      if( present( left_evec ) ) dimv1 = shape( left_evec )
      if( present( right_evec ) ) dimv2 = shape( right_evec )

      call assert( lgrad >= 0 .and. lgrad <= 3, &
        'Argument `left_gradient` must be either 0, 1, 2, or 3.' )
      call assert( rgrad >= 0 .and. rgrad <= 3, &
        'Argument `right_gradient` must be either 0, 1, 2, or 3.' )

      n1 = ngp1
      if( present( left_evec ) ) then
        call assert( dimv1(1) >= n1, &
          'First dimension of left eigenvector too small.' )
        n1 = dimv1(2)
      end if
      n2 = ngp2
      if( present( right_evec ) ) then
        call assert( dimv2(1) >= n2, &
          'First dimension of right eigenvector too small.' )
        n2 = dimv2(2)
      end if
      call assert( dimm(1) >= n1 .and. (diag .or. dimm(2) >= n2), &
        'Dimension of output matrix too small.' )

      if( nr1 ) then
        allocate( igpig1, source=Gpset1%igknrig(:, 1, ip1) )
        allocate( vgpc1, source=Gpset1%vgknrc(:, :, 1, ip1) )
      else
        allocate( igpig1, source=Gpset1%igkig(:, 1, ip1) )
        allocate( vgpc1, source=Gpset1%vgkc(:, :, 1, ip1) )
      end if
      if( nr2 ) then
        allocate( igpig2, source=Gpset2%igknrig(:, 1, ip2) )
        allocate( vgpc2, source=Gpset2%vgknrc(:, :, 1, ip2) )
      else
        allocate( igpig2, source=Gpset2%igkig(:, 1, ip2) )
        allocate( vgpc2, source=Gpset2%vgkc(:, :, 1, ip2) )
      end if

      ! no eigenvectors, diagonal only
      if( .not. present( left_evec ) .and. .not. present( right_evec ) .and. diag ) then
        allocate( auxmat1(min(ngp1, ngp2), 1) )
!$omp parallel default( shared ) private( igp1, ivg, ig )
!$omp do
        do igp1 = 1, min(ngp1, ngp2)
          ivg = Gset%ivg(:, igpig1(igp1)) - Gset%ivg(:, igpig2(igp1)) + ivg0
          ivg = modulo( ivg - Gset_op%intgv(:, 1), Gset_op%ngrid ) + Gset_op%intgv(:, 1)
          ig = Gset_op%ivgig(ivg(1), ivg(2), ivg(3))
          if( ig > Gset_op%ngvec ) then
            auxmat1(igp1, 1) = zzero
          else
            auxmat1(igp1, 1) = opig(ig)
            if( lgrad > 0 ) &
              auxmat1(igp1, 1) = auxmat1(igp1, 1) * cmplx( 0.0_dp, -vgpc1(lgrad, igp1), dp )
            if( rgrad > 0 ) &
              auxmat1(igp1, 1) = auxmat1(igp1, 1) * cmplx( 0.0_dp,  vgpc2(rgrad, igp1), dp )
          end if
        end do
!$omp end do
!$omp end parallel
        if( alpha /= zone ) auxmat1 = alpha * auxmat1
        mat(1:min(ngp1, ngp2), 1) = mat(1:min(ngp1, ngp2), 1) + auxmat1(:, 1)
        deallocate( auxmat1 )
      else
        allocate( auxmat1(ngp1, ngp2) )
!$omp parallel default( shared ) private( igp2, igp1, ivg, ig )
!$omp do collapse(2)
        do igp2 = 1, ngp2
          do igp1 = 1, ngp1
            ivg = Gset%ivg(:, igpig1(igp1)) - Gset%ivg(:, igpig2(igp2)) + ivg0
            ivg = modulo( ivg - Gset_op%intgv(:, 1), Gset_op%ngrid) + Gset_op%intgv(:, 1)
            ig = Gset_op%ivgig(ivg(1), ivg(2), ivg(3))
            if( ig > Gset_op%ngvec ) then
              auxmat1(igp1, igp2) = zzero
            else
              auxmat1(igp1, igp2) = opig(ig)
              if( lgrad > 0 ) &
                auxmat1(igp1, igp2) = auxmat1(igp1, igp2) * cmplx( 0.0_dp, -vgpc1(lgrad, igp1), dp )
              if( rgrad > 0 ) &
                auxmat1(igp1, igp2) = auxmat1(igp1, igp2) * cmplx( 0.0_dp,  vgpc2(rgrad, igp2), dp )
            end if
          end do
        end do
!$omp end do
!$omp end parallel
        if( alpha /= zone ) auxmat1 = alpha * auxmat1
        ! no eigenvectors, full matrix
        if( .not. present( left_evec ) .and. .not. present( right_evec ) ) then
          mat(1:ngp1, 1:ngp2) = mat(1:ngp1, 1:ngp2) + auxmat1
        else if( present( right_evec ) ) then
          allocate( auxmat2(ngp1, dimv2(2)) )
          call zgemm( 'n', 'n', ngp1, n2, ngp2, zone, &
                 auxmat1, ngp1, &
                 right_evec, dimv2(1), zzero, &
                 auxmat2, ngp1 )
          if( present( left_evec ) ) then
            ! left and right eigenvectors, diagonal only
            if( diag ) then
              do i = 1, min(n1, n2)
                mat(i, 1) = mat(i, 1) + zdotc( ngp1, left_evec(1, i), 1, auxmat2(1, i), 1 )
              end do
            ! left and right eigenvectors, all elements
            else
              call zgemm( 'c', 'n', n1, n2, ngp1, zone, &
                     left_evec, dimv1(1), &
                     auxmat2, ngp1, zone, &
                     mat, dimm(1) )
            end if
          else
            do igp1 = 1, ngp1
              ! right eigenvectors, diagnonal only
              if( diag ) then
                if( igp1 <= n2 ) mat(igp1, 1) = mat(igp1, 1) + auxmat2(igp1, igp1)
              ! right eigenvectors, all elements
              else
                do i = 1, n2
                  mat(igp1, i) = mat(igp1, i) + auxmat2(igp1, i)
                end do
              end if
            end do
          end if
          deallocate( auxmat2 )
        else if( present( left_evec ) ) then
          allocate( auxmat2(n1, ngp2) )
          call zgemm( 'c', 'n', n1, ngp2, ngp1, zone, &
                 left_evec, dimv1(1), &
                 auxmat1, ngp1, zzero, &
                 auxmat2, n1 )
          do igp2 = 1, ngp2
            ! left eigenvectors, diagnonal only
            if( diag ) then
              if( igp2 <= n1 ) mat(igp2, 1) = mat(igp2, 1) + auxmat2(igp2, igp2)
            ! left eigenvectors, all elements
            else
              do i = 1, n1
                mat(i, igp2) = mat(i, igp2) + auxmat2(i, igp2)
              end do
            end if
          end do
          deallocate( auxmat2 )
        end if
        deallocate( auxmat1 )
      end if

      deallocate( igpig1, igpig2, vgpc1, vgpc2 )
    end subroutine me_lapwlo_ir_mat

    ! END INTERSTITIAL INTEGRALS
    ! ******************************************
end module matrix_elements_lapw_lo
