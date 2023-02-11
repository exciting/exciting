!> title:   matrix elements
!> summary: Module that provides procedures to compute
!>          different types of matrix elements / integrals
!>          of local multiplicative operators
!> author:  Sebastian Tillack
!> date:    August 2022
!> licence: GPL
!>
!># Matrix elements in the (L)APW+LO basis
!>
!> This module provides functionalities to calculate matrix elements
!> of generic (local) operators in the (L)APW+LO basis.
!>
!> Given a local multiplicative operator \(O({\bf r})\), matrix elements of the type
!> \[ \langle \chi_{\mu} | O | \chi_{\nu} \rangle_{\mathcal{D}} =
!> \int\limits_\mathcal{D} \chi_{\mu}^\ast({\bf r}) \, O({\bf r}) \, \chi_{\nu}({\bf r}) \, {\rm d}\mathcal{D} \]
!> can be calculated.
!>
!> Here, \(\chi_{\mu}({\bf r})\) can be
!> 
!> *   (L)APW+LO basis functions, i.e., (L)APWs \(\phi_{\bf G+k}({\bf r})\) and local orbitals \(\phi_{L}({\bf r})\)
!> *   wavefunctions \(\psi_{n{\bf k}}({\bf r})\) expressed as linear combinations of (L)APW+LO basis functions
!> *   or any of their spatial derivatives \(\partial_i \chi_{\mu}({\bf r})\)
!>
!> The integration domain \(\mathcal{D}\) can be
!>
!> *   the interstitial region (IR)
!> *   the interior of an individual muffin-tin sphere, \(\alpha\)
!> *   the surface of an individual muffin-tin sphere, \(\partial \alpha\)
!>
!> @note 
!> This module only provides convenient user interfaces.
!> The underlying procedures are contained in [[matrix_elements_lapw_lo(module)]].
!> @endnote
!
! See 'matrix_elements.md' for further documentation and usage.
!>{!../src/matrix_elements/matrix_elements.md!}
!>
module matrix_elements
  use precision, only: dp

  implicit none
  private

  !> Allocate the wavevector independent part of muffin-tin matrix elements,
  !> i.e., the array that carries the radial integrals times Gaunt coefficients
  !>
  !> \[
  !> O^{\alpha}_{\lambda_1 \lambda_2} = \int\limits_\alpha \phi^{\alpha\ast}_{\lambda_1}({\bf r})\,
  !> O({\bf r})\, \phi^\alpha_{\lambda_2}({\bf r})\, {\rm d}^3 r =
  !> \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^l 
  !> \langle Y_{l_{\lambda_1} m_{\lambda_1}} | X_{lm} | Y_{l_{\lambda_2} m_{\lambda_2}} \rangle
  !> R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_{\lambda_1}, l_{\lambda_2}) \; .
  !> \]
  !> The array size will be determined by the number of muffin-tin basis functions \(\phi^\alpha_\lambda({\bf r})\)
  !> and the total number of atoms. The array will be initialized with zero. See also [[muffin_tin_basis(module)]].
  !>
  !> The interface is 
  !>
  !> * `me_mt_alloc( rignt )` for a single operator
  !> * `me_mt_alloc( rignt, n )` for `n` operators
  interface me_mt_alloc
    procedure :: mt_alloc_single, mt_alloc_multi
  end interface

  !> Allocate the wavevector independent part of interstitial matrix elements,
  !> i.e., the reciprocal space representation of the product of the operator
  !> and the characteristic function
  !>
  !> \[
  !> O({\bf r})\, \Theta_{\rm IR}({\bf r}) = \sum\limits_{\bf G} \hat{O}_\Theta({\bf G})\, 
  !> {\rm e}^{{\rm i} {\bf G} \cdot {\bf r}} \; .
  !> \]
  !> The array size will be determined by the number of \({\bf G}\)-vectors with 
  !> \(|{\bf G}|<G_{\rm max}\) corresponding to the input argument `Gset` passed to [[me_init(subroutine)]].
  !> The array will be initialized with zero.
  !>
  !> The interface is 
  !>
  !> * `me_ir_alloc( opig )` for a single operator
  !> * `me_ir_alloc( opig, n )` for `n` operators
  interface me_ir_alloc
    procedure :: ir_alloc_single, ir_alloc_multi
  end interface

  !> Prepare the wavevector independent part of muffin-tin matrix elements,
  !> i.e., compute the radial integrals times Gaunt coefficients for a given
  !> operator and muffin-tin \(\alpha\)
  !>
  !> \[
  !> O^{\alpha}_{\lambda_1 \lambda_2} = \int\limits_\mathcal{D} \phi^{\alpha\ast}_{\lambda_1}({\bf r})\,
  !> O({\bf r})\, \phi^\alpha_{\lambda_2}({\bf r})\, {\rm d}\mathcal{D} =
  !> \sum_{l=0}^{l_{\rm max,o}} \sum_{m=-l}^l 
  !> \langle Y_{l_{\lambda_1} m_{\lambda_1}} | X_{lm} | Y_{l_{\lambda_2} m_{\lambda_2}} \rangle
  !> R^{\alpha}_{\tilde{\lambda}_1 \tilde{\lambda}_2, lm}(l_{\lambda_1}, l_{\lambda_2}) \; .
  !> \]
  !> Use the optional arguments `left_radial_derivative` and `right_radial_derivative` to use
  !> \(\frac{{\rm d}^k}{{\rm d}r^k} g^\alpha_\tilde{\lambda}(r)\) instead.
  !>
  !> Use the optional arguments `left_gradient` and `right_gradient` to use 
  !> \({\bf \nabla}_i \phi^\alpha_\lambda({\bf r})\) instead.
  !>
  !> Use the optional argument `surface_integral` to compute integrals over the muffin-tin surface
  !> element \({\rm d}\mathcal{D}=\hat{\bf e}_i {\rm d}S\) instead.
  !>
  !> This routine computes
  !> \[ O^{\alpha,{\rm out}}_{\lambda_1 \lambda_2} = 
  !>    a\, O^{\alpha}_{\lambda_1 \lambda_2} + b\, O^{\alpha,{\rm in}}_{\lambda_1 \lambda_2} \; . \]
  !>
  !> The interface is
  !> `me_mt_prepare( is, ias, lmax, a, rfun, b, rignt 
  !>   [, left_radial_derivative, right_radial_derivative, left_gradient, right_gradient, 
  !>    surface integral] )`, 
  !> where `rfun` and `rignt` can be real or complex valued.
  !>
  !> @note Fore more flexible use, employ [[me_lapwlo_mt_rignt(subroutine)]] directly. @endnote
  interface me_mt_prepare
    procedure :: mt_rignt_real, mt_rignt_complex
  end interface

  !> Prepare the wavevector independent part of interstitial matrix elements,
  !> i.e., compute the reciprocal space representation of the product of the 
  !> operator and the characteristic function 
  !>
  !> \[
  !> O({\bf r})\, \Theta_{\rm IR}({\bf r}) = \sum\limits_{\bf G} \hat{O}_\Theta({\bf G})\, 
  !> {\rm e}^{{\rm i} {\bf G} \cdot {\bf r}} \; .
  !> \]
  !>
  !> This routine computes
  !> \[ \hat{O}^{\rm out}_\Theta({\bf G}) = a\, \hat{O}_\Theta({\bf G}) + b\, \hat{O}^{\rm in}_\Theta({\bf G}) \;. \]
  !>
  !> The interface is 
  !> `me_ir_prepare( a, opir, b, opig [, Gset_op] )`, 
  !> where `opir` can be real or complex valued.
  !>
  !> @note Fore more flexible use, employ [[me_lapwlo_ir_opig(subroutine)]] directly. @endnote
  interface me_ir_prepare
    procedure :: ir_opig_real, ir_opig_complex
  end interface

  !> Having prepared the \({\bf p}\) independent radial integrals times Gaunt coefficients
  !> (`rignt`) using [[me_mt_prepare(subroutine)]], this subroutine computes the 
  !> contribution to the matrix elements coming from the given MT
  !>
  !> \[ M^\alpha_{\mu \nu} = \langle \chi_\mu | O | \chi_\nu \rangle_\mathcal{D} \; .\]
  !>
  !> This routine computes
  !> \[ M^{\alpha,{\rm out}}_{\mu \nu} = a\, M^\alpha_{\mu \nu} + b\, M^{\alpha,{\rm in}}_{\mu \nu} \; . \]
  !>
  !> The following interfaces are allowed:
  !>
  !> *   On both sides \(\chi({\bf r})\) are (L)APW+LO basis functions.  
  !>     *   Same \({\bf p}\)-point on both sides:  
  !>         `me_mt_mat( is, ias, ngp, apwalm, a, rignt, b, mat [, diagonal_only, left_local_orbitals, right_local_orbitals] )`
  !>     *   Different \({\bf p}\)-points on both sides:  
  !>         `me_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, a, rignt, b, mat [, diagonal_only, left_local_orbitals, right_local_orbitals] )`
  !> *   On both sides \(\chi({\bf r})\) are wavefunctions.  
  !>     *   Same \({\bf p}\)-point on both sides:  
  !>         `me_mt_mat( is, ias, ngp, apwalm, evec, a, rignt, b, mat [, diagonal_only, left_local_orbitals, right_local_orbitals] )`
  !>     *   Different \({\bf p}\)-points on both sides:
  !>         `me_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, evec1, evec2, a, rignt, b, mat [, diagonal_only, left_local_orbitals, right_local_orbitals] )`
  !>
  !> Use the optional argument `diagonal_only` to compute only diagonal matrix elements \(M^\alpha_{\mu \mu}\).
  !>
  !> Use the optional arguments `left_local_orbitals` and `right_local_orbitals` to include or exclude LOs from the calculation.
  !>
  !> @note For more flexible usage, employ [[me_lapwlo_mt_mat(subroutine)]] directly. @endnote
  interface me_mt_mat
    procedure :: mt_mat_apwalm_real, mt_mat_apwalm_complex, &
                 mt_mat_evec_real, mt_mat_evec_complex, &
                 mt_mat_apwalm_real_single, mt_mat_apwalm_complex_single, &
                 mt_mat_evec_real_single, mt_mat_evec_complex_single
  end interface

  !> Having prepared the \({\bf p}\) independent interstitial representation of the operator
  !> (`opig`) using [[me_ir_prepare(subroutine)]], this subroutine computes the 
  !> contribution to the matrix elements coming from the interstitial region
  !>
  !> \[ M^{\rm IR}_{\mu \nu} = \langle \chi_\mu | O | \chi_\nu \rangle_{\rm IR} \; . \]
  !>
  !> This routine computes
  !> \[ M^{{\rm IR},{\rm out}}_{\mu \nu} = a\, M^{\rm IR}_{\mu \nu} + b\, M^{{\rm IR},{\rm in}}_{\mu \nu} \; . \]
  !>
  !> The following interfaces are allowed:
  !>
  !> *   On both sides \(\chi({\bf r})\) are plane wave basis functions.  
  !>     *   Same \({\bf p}\)-point on both sides:  
  !>         `me_ir_mat( Gpset, ip, a, opig, b, mat [, Gset_op, diagonal_only] )`
  !>     *   Different \({\bf p}\)-points on both sides:  
  !>         `me_ir_mat( Gpset1, ip1, Gpset2, ip2, a, opig, b, mat [, Gset_op, diagonal_only] )`
  !> *   On both sides \(\chi({\bf r})\) are wavefunctions.  
  !>     *   Same \({\bf p}\)-point on both sides:  
  !>         `me_ir_mat( Gpset, ip, evec, a, opig, b, mat [, Gset_op, diagonal_only] )`
  !>     *   Different \({\bf p}\)-points on both sides:  
  !>         `me_ir_mat( Gpset1, ip1, Gpset2, ip2, evec1, evec2, a, opig, b, mat [, Gset_op, diagonal_only] )`
  !>
  !> Use the optional argument `diagonal_only` to compute only diagonal matrix elements \(M^\alpha_{\mu \mu}\).
  !>
  !> Use the optional arguments `left_gradient` and `right_gradient` to use 
  !> \({\bf \nabla}_i \chi_\mu({\bf r})\) instead.
  !>
  !> @note For more flexible usage, employ [[me_lapwlo_ir_mat(subroutine)]] directly. @endnote
  interface me_ir_mat
    procedure :: ir_mat_basis, ir_mat_evec, &
                 ir_mat_basis_single, ir_mat_evec_single
  end interface

  public :: me_init, me_finit, &
            me_mt_alloc, me_ir_alloc, &
            me_mt_prepare, me_ir_prepare, &
            me_mt_mat, me_ir_mat

  contains

    !> Initialize the module.
    !> Sets module variables and Gaunt coefficients (if necessary).
    subroutine me_init( basis, lmax_operator, Gset )
      use matrix_elements_lapw_lo, only: me_lapwlo_gen_helper_arrays
      use muffin_tin_basis, only: mt_basis_type
      use mod_kpointset, only: G_set
      !> (L)APW and LO radial functions
      type(mt_basis_type), intent(in) :: basis
      !> maximum \(l\) used in operator expansion
      integer, intent(in) :: lmax_operator
      !> set of \({\bf G}\)-vectors for Fourier series
      type(G_set), target, intent(in) :: Gset

      call me_lapwlo_gen_helper_arrays( basis, lmax_operator, Gset )
    end subroutine me_init
    
    !> Finalize the module and free memory.
    subroutine me_finit
      use matrix_elements_lapw_lo, only: me_lapwlo_del_helper_arrays
      call me_lapwlo_del_helper_arrays
    end subroutine me_finit

    ! ******************************************
    ! INTERFACE ME_MT_ALLOC

    !> See [[me_mt_alloc]].
    subroutine mt_alloc_single( rignt )
      use matrix_elements_lapw_lo, only: basis
      use constants, only: zzero
      use mod_atoms, only: natmtot
      !> array for radial integrals times Gaunt coefficients
      !> for a single operator given in complex spherical harmonics
      complex(dp), allocatable, intent(inout) :: rignt(:,:,:)

      if( allocated( rignt ) ) deallocate( rignt )
      allocate( rignt(basis%n_basis_fun_max, basis%n_basis_fun_max, natmtot), source=zzero )
    end subroutine mt_alloc_single

    !> See [[me_mt_alloc]].
    subroutine mt_alloc_multi( rignt, n )
      use matrix_elements_lapw_lo, only: basis
      use constants, only: zzero
      use mod_atoms, only: natmtot
      !> array for radial integrals times Gaunt coefficients
      !> for multiple operators given in complex spherical harmonics
      complex(dp), allocatable, intent(inout) :: rignt(:,:,:,:)
      !> number of operators
      integer, intent(in) :: n

      if( allocated( rignt ) ) deallocate( rignt )
      allocate( rignt(basis%n_basis_fun_max, basis%n_basis_fun_max, natmtot, n), source=zzero )
    end subroutine mt_alloc_multi

    ! END INTERFACE ME_MT_ALLOC
    ! ******************************************

    ! ******************************************
    ! INTERFACE ME_IR_ALLOC

    !> See [[me_ir_alloc]].
    subroutine ir_alloc_single( opig )
      use matrix_elements_lapw_lo, only: Gset
      use constants, only: zzero
      !> array for reciprocal space representation of operator times
      !> characteristic function for a single operator
      complex(dp), allocatable, intent(inout) :: opig(:)

      if( allocated( opig ) ) deallocate( opig )
      allocate( opig(Gset%ngvec), source=zzero )
    end subroutine ir_alloc_single

    !> See [[me_ir_alloc]].
    subroutine ir_alloc_multi( opig, n )
      use matrix_elements_lapw_lo, only: Gset
      use constants, only: zzero
      !> array for reciprocal space representation of operator times
      !> characteristic function for multiple operators
      complex(dp), allocatable, intent(inout) :: opig(:,:)
      !> number of operators
      integer, intent(in) :: n

      if( allocated( opig ) ) deallocate( opig )
      allocate( opig(Gset%ngvec, n), source=zzero )
    end subroutine ir_alloc_multi
    
    ! END INTERFACE ME_IR_ALLOC
    ! ******************************************

    ! ******************************************
    ! INTERFACE ME_MT_PREPARE

    !> See [[me_mt_prepare]].
    subroutine mt_rignt_real( is, ias, lmax, alpha, rrfun, beta, zrignt, &
        left_radial_derivative, right_radial_derivative, left_gradient, right_gradient, surface_integral )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_rignt
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> maximum l for operator expansion
      integer, intent(in) :: lmax
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> real radial functions of operator
      real(dp), intent(in) :: rrfun(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> complex result array. 
      !> See also [[me_mt_alloc(subroutine)]].
      complex(dp), intent(inout) :: zrignt(:,:)
      !> order of radial derivative (0 or 1) for left basis functions (default: 0)
      integer, optional, intent(in) :: left_radial_derivative
      !> order of radial derivative (0 or 1) for right basis functions (default: 0)
      integer, optional, intent(in) :: right_radial_derivative
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient
      !> Cartesian component (1, 2, or 3) of surface normal for surface integrals (default: 0 - volume integral)
      integer, optional, intent(in) :: surface_integral

      integer :: lr, rr, lg, rg, sf

      lr = 0; if( present( left_radial_derivative ) ) lr = left_radial_derivative
      rr = 0; if( present( right_radial_derivative ) ) rr = right_radial_derivative
      lg = 0; if( present( left_gradient ) ) lg = left_gradient
      rg = 0; if( present( right_gradient ) ) rg = right_gradient
      sf = 0; if( present( surface_integral ) ) sf = surface_integral

      call me_lapwlo_mt_rignt( is, ias, lmax, alpha, rrfun, beta, zrignt, lr, rr, lg, rg, sf, .true. )
    end subroutine mt_rignt_real

    !> See [[me_mt_prepare]].
    subroutine mt_rignt_complex( is, ias, lmax, alpha, zrfun, beta, zrignt, &
        left_radial_derivative, right_radial_derivative, left_gradient, right_gradient, surface_integral, real_expansion )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_rignt
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> maximum l for operator expansion
      integer, intent(in) :: lmax
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> complex radial functions of operator
      complex(dp), intent(in) :: zrfun(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> complex result array. 
      !> See also [[me_mt_alloc(subroutine)]].
      complex(dp), intent(inout) :: zrignt(:,:)
      !> order of radial derivative (0 or 1) for left basis functions (default: 0)
      integer, optional, intent(in) :: left_radial_derivative
      !> order of radial derivative (0 or 1) for right basis functions (default: 0)
      integer, optional, intent(in) :: right_radial_derivative
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient
      !> Cartesian component (1, 2, or 3) of surface normal for surface integrals (default: 0 - volume integral)
      integer, optional, intent(in) :: surface_integral
      !> operator is given as an real spherical harmonics expansion (default: `.false.`)
      logical, optional, intent(in) :: real_expansion

      integer :: lr, rr, lg, rg, sf
      logical :: real

      lr = 0; if( present( left_radial_derivative ) ) lr = left_radial_derivative
      rr = 0; if( present( right_radial_derivative ) ) rr = right_radial_derivative
      lg = 0; if( present( left_gradient ) ) lg = left_gradient
      rg = 0; if( present( right_gradient ) ) rg = right_gradient
      sf = 0; if( present( surface_integral ) ) sf = surface_integral
      real = .false.; if( present( real_expansion ) ) real = real_expansion

      call me_lapwlo_mt_rignt( is, ias, lmax, alpha, zrfun, beta, zrignt, lr, rr, lg, rg, sf, real )
    end subroutine mt_rignt_complex

    ! END INTERFACE ME_MT_PREPARE
    ! ******************************************

    ! ******************************************
    ! INTERFACE ME_IR_PREPARE
      
    !> See [[me_ir_prepare]].
    subroutine ir_opig_real( alpha, ropir, beta, opig, Gset_op )
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_opig, Gset, cfunr
      use mod_kpointset, only: G_set
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> real operator in real space 
      !> (on FFT grid corresponding to input argument `Gset_op`)
      real(dp), intent(in) :: ropir(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> reciprocal space representation of the product
      !> (on \({\bf G}\)-vectors in `Gset_op`)
      complex(dp), intent(inout) :: opig(:)
      !> set of \({\bf G}\)-vectors the operator is defined on 
      !> (default: input argument `Gset` passed to [[me_init(subroutine)]])
      type(G_set), optional, intent(in) :: Gset_op
      
      if( present( Gset_op ) ) then
        call me_lapwlo_ir_opig( alpha, ropir, Gset_op, cfunr, beta, opig, Gset_op )
      else
        call me_lapwlo_ir_opig( alpha, ropir, Gset, cfunr, beta, opig, Gset )
      end if
    end subroutine ir_opig_real

    !> See [[me_ir_prepare]].
    subroutine ir_opig_complex( alpha, zopir, beta, opig, Gset_op )
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_opig, Gset, cfunr
      use mod_kpointset, only: G_set
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> complex operator in real space 
      !> (on FFT grid corresponding to input argument `Gset_op`)
      complex(dp), intent(in) :: zopir(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> reciprocal space representation of the product
      !> (on \({\bf G}\)-vectors in `Gset_op`)
      complex(dp), intent(inout) :: opig(:)
      !> set of \({\bf G}\)-vectors the operator is defined on 
      !> (default: input argument `Gset` passed to [[me_init(subroutine)]])
      type(G_set), optional, intent(in) :: Gset_op
      
      if( present( Gset_op ) ) then
        call me_lapwlo_ir_opig( alpha, zopir, Gset_op, cfunr, beta, opig, Gset_op )
      else
        call me_lapwlo_ir_opig( alpha, zopir, Gset, cfunr, beta, opig, Gset )
      end if
    end subroutine ir_opig_complex

    ! END INTERFACE ME_IR_PREPARE
    ! ******************************************
      
    ! ******************************************
    ! INTERFACE ME_MT_MAT

    !> See [[me_mt_mat]].
    subroutine mt_mat_apwalm_real( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, rrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
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
      !> real radial integrals times Gaunt coefficients
      real(dp), intent(in) :: rrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, rrignt, beta, mat, &
             diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_apwalm_real

    !> See [[me_mt_mat]].
    subroutine mt_mat_apwalm_complex( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, zrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
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
      !> complex radial integrals times Gaunt coefficients
      complex(dp), intent(in) :: zrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, zrignt, beta, mat, &
             diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_apwalm_complex

    !> See [[me_mt_mat]].
    subroutine mt_mat_evec_real( is, ias, ngp1, ngp2, apwalm1, apwalm2, evec1, evec2, alpha, rrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals)
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
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
      !> eigenvectors on the left
      complex(dp), intent(in) :: evec1(:,:)
      !> eigenvectors on the right
      complex(dp), intent(in) :: evec2(:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> real radial integrals times Gaunt coefficients
      real(dp), intent(in) :: rrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, rrignt, beta, mat, &
             left_evec=evec1, right_evec=evec2, diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_evec_real

    !> See [[me_mt_mat]].
    subroutine mt_mat_evec_complex( is, ias, ngp1, ngp2, apwalm1, apwalm2, evec1, evec2, alpha, zrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
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
      !> eigenvectors on the left
      complex(dp), intent(in) :: evec1(:,:)
      !> eigenvectors on the right
      complex(dp), intent(in) :: evec2(:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> complex radial integrals times Gaunt coefficients
      complex(dp), intent(in) :: zrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp1, ngp2, apwalm1, apwalm2, alpha, zrignt, beta, mat, &
             left_evec=evec1, right_evec=evec2, diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_evec_complex

    !> See [[me_mt_mat]].
    subroutine mt_mat_apwalm_real_single( is, ias, ngp, apwalm, alpha, rrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> number of \({\bf G+p}\)-vectors
      integer, intent(in) :: ngp
      !> (L)APW matching coefficients
      complex(dp), intent(in) :: apwalm(:,:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> real radial integrals times Gaunt coefficients
      real(dp), intent(in) :: rrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp, ngp, apwalm, apwalm, alpha, rrignt, beta, mat, &
             diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_apwalm_real_single

    !> See [[me_mt_mat]].
    subroutine mt_mat_apwalm_complex_single( is, ias, ngp, apwalm, alpha, zrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> number of \({\bf G+p}\)-vectors
      integer, intent(in) :: ngp
      !> (L)APW matching coefficients
      complex(dp), intent(in) :: apwalm(:,:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> complex radial integrals times Gaunt coefficients
      complex(dp), intent(in) :: zrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp, ngp, apwalm, apwalm, alpha, zrignt, beta, mat, &
             diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_apwalm_complex_single

    !> See [[me_mt_mat]].
    subroutine mt_mat_evec_real_single( is, ias, ngp, apwalm, evec, alpha, rrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> number of \({\bf G+p}\)-vectors
      integer, intent(in) :: ngp
      !> (L)APW matching coefficients
      complex(dp), intent(in) :: apwalm(:,:,:)
      !> eigenvectors
      complex(dp), intent(in) :: evec(:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> real radial integrals times Gaunt coefficients
      real(dp), intent(in) :: rrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp, ngp, apwalm, apwalm, alpha, rrignt, beta, mat, &
             left_evec=evec, right_evec=evec, diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_evec_real_single

    !> See [[me_mt_mat]].
    subroutine mt_mat_evec_complex_single( is, ias, ngp, apwalm, evec, alpha, zrignt, beta, mat, &
        diagonal_only, left_local_orbitals, right_local_orbitals )
      use matrix_elements_lapw_lo, only: me_lapwlo_mt_mat
      !> index of the species of the MT
      integer, intent(in) :: is
      !> index of the atom of the MT
      integer, intent(in) :: ias
      !> number of \({\bf G+p}\)-vectors
      integer, intent(in) :: ngp
      !> (L)APW matching coefficients
      complex(dp), intent(in) :: apwalm(:,:,:)
      !> eigenvectors
      complex(dp), intent(in) :: evec(:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> complex radial integrals times Gaunt coefficients
      complex(dp), intent(in) :: zrignt(:,:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> include local orbitals on the left (default: `.true.`)
      logical, optional, intent(in) :: left_local_orbitals
      !> include local orbitals on the right (default: `.true.`)
      logical, optional, intent(in) :: right_local_orbitals

      logical :: diag, lo1, lo2

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lo1 = .true.; if( present( left_local_orbitals ) ) lo1 = left_local_orbitals
      lo2 = .true.; if( present( right_local_orbitals ) ) lo2 = right_local_orbitals

      call me_lapwlo_mt_mat( is, ias, ngp, ngp, apwalm, apwalm, alpha, zrignt, beta, mat, &
             left_evec=evec, right_evec=evec, diagonal_only=diag, left_local_orbitals=lo1, right_local_orbitals=lo2 )
    end subroutine mt_mat_evec_complex_single

    ! END INTERFACE ME_MT_MAT
    ! ******************************************

    ! ******************************************
    ! INTERFACE ME_IR_MAT
    
    !> See [[me_ir_mat]].
    subroutine ir_mat_basis( Gpset1, ip1, Gpset2, ip2, alpha, opig, beta, mat, &
        Gset_op, diagonal_only, left_gradient, right_gradient )
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_mat, Gset
      use mod_kpointset, only: G_set, Gk_set
      !> set of \({\bf G+p}\) vectors for basis functions on the left
      type(Gk_set), intent(in) :: Gpset1
      !> index of \({\bf p}\)-point on the left
      integer, intent(in) :: ip1
      !> set of \({\bf G+p}\) vectors for basis functions on the right
      type(Gk_set), intent(in) :: Gpset2
      !> index of \({\bf p}\)-point on the right
      integer, intent(in) :: ip2
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> operator (times characteristic function) in reciprocal space (see [[me_ir_prepare(subroutine)]])
      complex(dp), intent(in) :: opig(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> set of \({\bf G}\)-vectors on which the operator is defined
      type(G_set), optional, intent(in) :: Gset_op
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient

      logical :: diag
      integer :: lg, rg

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lg = 0; if( present( left_gradient ) ) lg = left_gradient
      rg = 0; if( present( right_gradient ) ) rg = right_gradient

      if( present( Gset_op ) ) then
        call me_lapwlo_ir_mat( Gpset1, ip1, Gpset2, ip2, Gset_op, alpha, opig, beta, mat, &
          diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      else
        call me_lapwlo_ir_mat( Gpset1, ip1, Gpset2, ip2, Gset, alpha, opig, beta, mat, &
          diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      end if
    end subroutine ir_mat_basis

    !> See [[me_ir_mat]].
    subroutine ir_mat_evec( Gpset1, ip1, Gpset2, ip2, evec1, evec2, alpha, opig, beta, mat, &
        Gset_op, diagonal_only, left_gradient, right_gradient )
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_mat, Gset
      use mod_kpointset, only: G_set, Gk_set
      !> set of \({\bf G+p}\) vectors for basis functions on the left
      type(Gk_set), intent(in) :: Gpset1
      !> index of \({\bf p}\)-point on the left
      integer, intent(in) :: ip1
      !> set of \({\bf G+p}\) vectors for basis functions on the right
      type(Gk_set), intent(in) :: Gpset2
      !> index of \({\bf p}\)-point on the right
      integer, intent(in) :: ip2
      !> eigenvectors on the left
      complex(dp), intent(in) :: evec1(:,:)
      !> eigenvectors on the right
      complex(dp), intent(in) :: evec2(:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> operator (times characteristic function) in reciprocal space (see [[me_ir_prepare(subroutine)]])
      complex(dp), intent(in) :: opig(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> set of \({\bf G}\)-vectors on which the operator is defined
      type(G_set), optional, intent(in) :: Gset_op
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient

      logical :: diag
      integer :: lg, rg

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lg = 0; if( present( left_gradient ) ) lg = left_gradient
      rg = 0; if( present( right_gradient ) ) rg = right_gradient

      if( present( Gset_op ) ) then
        call me_lapwlo_ir_mat( Gpset1, ip1, Gpset2, ip2, Gset_op, alpha, opig, beta, mat, &
               left_evec=evec1, right_evec=evec2, diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      else
        call me_lapwlo_ir_mat( Gpset1, ip1, Gpset2, ip2, Gset, alpha, opig, beta, mat, &
               left_evec=evec1, right_evec=evec2, diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      end if
    end subroutine ir_mat_evec

    !> See [[me_ir_mat]].
    subroutine ir_mat_basis_single( Gpset, ip, alpha, opig, beta, mat, &
        Gset_op, diagonal_only, left_gradient, right_gradient )
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_mat, Gset
      use mod_kpointset, only: G_set, Gk_set
      !> set of \({\bf G+p}\) vectors for basis functions
      type(Gk_set), intent(in) :: Gpset
      !> index of \({\bf p}\)-point
      integer, intent(in) :: ip
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> operator (times characteristic function) in reciprocal space (see [[me_ir_prepare(subroutine)]])
      complex(dp), intent(in) :: opig(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> set of \({\bf G}\)-vectors on which the operator is defined
      type(G_set), optional, intent(in) :: Gset_op
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient

      logical :: diag
      integer :: lg, rg

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lg = 0; if( present( left_gradient ) ) lg = left_gradient
      rg = 0; if( present( right_gradient ) ) rg = right_gradient

      if( present( Gset_op ) ) then
        call me_lapwlo_ir_mat( Gpset, ip, Gpset, ip, Gset_op, alpha, opig, beta, mat, &
          diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      else
        call me_lapwlo_ir_mat( Gpset, ip, Gpset, ip, Gset, alpha, opig, beta, mat, &
          diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      end if
    end subroutine ir_mat_basis_single

    !> See [[me_ir_mat]].
    subroutine ir_mat_evec_single( Gpset, ip, evec, alpha, opig, beta, mat, &
        Gset_op, diagonal_only, left_gradient, right_gradient )
      use matrix_elements_lapw_lo, only: me_lapwlo_ir_mat, Gset
      use mod_kpointset, only: G_set, Gk_set
      !> set of \({\bf G+p}\) vectors for basis functions
      type(Gk_set), intent(in) :: Gpset
      !> index of \({\bf p}\)-point
      integer, intent(in) :: ip
      !> eigenvectors
      complex(dp), intent(in) :: evec(:,:)
      !> prefactor \(a\)
      complex(dp), intent(in) :: alpha
      !> operator (times characteristic function) in reciprocal space (see [[me_ir_prepare(subroutine)]])
      complex(dp), intent(in) :: opig(:)
      !> prefactor \(b\)
      complex(dp), intent(in) :: beta
      !> matrix elements
      complex(dp), intent(inout) :: mat(:,:)
      !> set of \({\bf G}\)-vectors on which the operator is defined
      type(G_set), optional, intent(in) :: Gset_op
      !> compute only diagonal matrix elements, stored in first index of `mat` (default: `.false.`)
      logical, optional, intent(in) :: diagonal_only
      !> Cartesian component (1, 2, or 3) of gradient for left basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: left_gradient
      !> Cartesian component (1, 2, or 3) of gradient for right basis functions (default: 0 - no gradient)
      integer, optional, intent(in) :: right_gradient

      logical :: diag
      integer :: lg, rg

      diag = .false.; if( present( diagonal_only ) ) diag = diagonal_only
      lg = 0; if( present( left_gradient ) ) lg = left_gradient
      rg = 0; if( present( right_gradient ) ) rg = right_gradient

      if( present( Gset_op ) ) then
        call me_lapwlo_ir_mat( Gpset, ip, Gpset, ip, Gset_op, alpha, opig, beta, mat, &
               left_evec=evec, right_evec=evec, diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      else
        call me_lapwlo_ir_mat( Gpset, ip, Gpset, ip, Gset, alpha, opig, beta, mat, &
               left_evec=evec, right_evec=evec, diagonal_only=diag, left_gradient=lg, right_gradient=rg )
      end if
    end subroutine ir_mat_evec_single

    ! END INTERFACE ME_IR_MAT
    ! ******************************************
end module matrix_elements
