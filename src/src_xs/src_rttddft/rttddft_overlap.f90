!> Module that manages the overlap matrix in RT-TDDFT
module rttddft_Overlap
#include "asserts.fpp"
  use constants, only: y00, zi, zone, zzero
  use math_utils, only: all_zero
  use matrix_elements, only: me_mt_alloc, me_mt_prepare, me_mt_mat, me_ir_mat
  use mod_atoms, only: nspecies, natoms, idxas
  use mod_gvector, only: cfunig
  use mod_kpointset, only: Gk_set
  use mod_muffin_tin, only: nrmtmax
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_arrays, only: positive_matrix_set
  use rttddft_timings, only: timesec_RTTDDFT, Timing_RTTDDFT_overlap
  use rttddft_VectorPotential, only: Vector_Potential_Field

  implicit none

  private
  integer(i32), parameter :: n_cartesian = 3

  !> Type to encapsulate the set of overlap matrices
  type, public, extends(positive_matrix_set) :: overlap_set
    private
    logical :: identity = .false.
  contains
    private
    procedure, public :: allocate => overlap_set_allocate
    procedure, public :: assert_is_identity => overlap_set_assert_is_identity
    procedure, public :: calculate => overlap_set_calculate
    procedure, public :: initialize => overlap_set_initialize
    procedure, public :: is_identity => overlap_set_is_identity
    procedure, public :: is_not_identity => overlap_set_is_not_identity
    procedure :: initialize_as_identity => overlap_set_initialize_as_identity
    final :: destructor
  end type
    
contains
  !> Wrapper for calling `this%allocate_array`
  pure subroutine overlap_set_allocate( this, use_lapwlo_basis, m, ki, kf )
    class(overlap_set), intent(inout) :: this
    !> If `.true.`, LAPW+LO basis is used (allocate matrix). Otherwise, overlap is identity
    logical, intent(in) :: use_lapwlo_basis
    !> Overlap matrix dimension
    integer(i32), intent(in) :: m
    !> First \( \mathbf{k} \)-point
    integer(i32), intent(in) :: ki
    !> Last \( \mathbf{k} \)-point
    integer(i32), intent(in) :: kf

    this%identity = .not. use_lapwlo_basis 
    call this%allocate_array( [1, 1 , ki], [m, m, kf] )
  end subroutine

  subroutine overlap_set_assert_is_identity( this )
    class(overlap_set), intent(in) :: this

    CALL_ASSERT( this%is_identity(), "Overlap matrix is not identity" )
  end subroutine

  pure elemental logical function overlap_set_is_identity( this )
    class(overlap_set), intent(in) :: this

    overlap_set_is_identity = this%identity
  end function

  pure elemental logical function overlap_set_is_not_identity( this )
    class(overlap_set), intent(in) :: this

    overlap_set_is_not_identity = .not. this%identity
  end function

  pure elemental subroutine destructor( this )
    type(overlap_set), intent(inout) :: this

    call this%deallocate_if_allocated()
  end subroutine

  !> Initialize the overlap matrix set
  subroutine overlap_set_initialize( this, apwalm, Gkset, p_MT, a_tot, t_overlap, mathcalH, mathcalB )
    class(overlap_set), intent(inout) :: this
    !> See [[overlap_set_calculate]]
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> See [[overlap_set_calculate]]
    type(Gk_set), intent(in) :: Gkset
    !> See [[overlap_set_calculate]]
    complex(dp), contiguous, optional, intent(in) :: p_MT(:, :, :, :, :)
    !> See [[overlap_set_calculate]]
    type(Timing_RTTDDFT_overlap), optional, intent(inout) :: t_overlap
    !> See [[overlap_set_calculate]]
    type(Vector_Potential_Field), optional, intent(in) :: a_tot
    !> See [[overlap_set_calculate]]
    complex(dp), contiguous, optional, intent(out) :: mathcalH(:, :, :, :, :)
    !> See [[overlap_set_calculate]]
    complex(dp), contiguous, optional, intent(out) :: mathcalB(:, :, :, :, :)

    if( this%is_not_identity() ) then
      call this%calculate( apwalm, Gkset, p_MT, a_tot, t_overlap, mathcalH, mathcalB )
    else
      call this%initialize_as_identity()
    end if
  end subroutine

  subroutine overlap_set_initialize_as_identity( this )
    class(overlap_set), intent(inout) :: this

    integer(i32) :: i, k

    this%array = zzero
    ! TODO: use DO CONCURRENT here after ifort2021 support is dropped
    do k = lbound( this%array, 3 ), ubound( this%array, 3 )
      do i = 1, size( this%array, 1 )
        this%array(i, i, k) = zone
      end do
    end do
  end subroutine

  !> Calculate the overlap matrix in the (L)APW + LO basis
  subroutine overlap_set_calculate( S, apwalm, Gkset, p_MT, a_tot, t_overlap, mathcalH, mathcalB )
    !> Overlap matrix (of basis functions)
    class(overlap_set), intent(inout) :: S
    !> Matching coefficients of the (L)APWs (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> Set of \( \mathbf{G} + \mathbf{k} \) vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Muffin-tin part of the Momentum matrix (nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: p_MT(:, :, :, :, :)
    !> Object that packs information about timings to update the overlap matrix
    type(Timing_RTTDDFT_overlap), optional, intent(inout) :: t_overlap
    !> Total vector potential (needed if `mathcalH` is present)
    type(Vector_Potential_Field), optional, intent(in) :: a_tot
    !> Auxiliary matrix needed to evaluate force corrections in MD calculations, see [[rttddft_GlobalMDVariables]]
    complex(dp), contiguous, optional, intent(out) :: mathcalH(:, :, :, :, :)
    !> Auxiliary matrix needed to evaluate force corrections in MD calculations, see [[rttddft_GlobalMDVariables]]
    complex(dp), contiguous, optional, intent(out) :: mathcalB(:, :, :, :, :)

    integer(i32) :: ik, is, ia, ias, first_kpt, last_kpt, shift
    integer(i32), parameter :: l_max = 0 ! Overlap operator is equal to (1.0/y00)*Y_{00}, it has only l=0 component
    integer(i32), parameter :: lm_max = (l_max+1)**2
    real(dp) :: ti
    real(dp), allocatable :: rfun(:, :)
    complex(dp), allocatable :: mt_part(:, :, :)

    if( S%is_identity() ) return
    first_kpt = lbound( S%array, 3 )
    last_kpt = ubound( S%array, 3 )
    shift = first_kpt - 1

    ! Check optional arguments
    if( present( mathcalH ) ) then
      CALL_ASSERT( present(a_tot), 'a_tot must be present' )
    end if
    if( present( mathcalH ) .or. present( mathcalB ) ) then
      CALL_ASSERT( present(p_MT), 'p_MT must be present' )
    end if
    if( present( t_overlap ) ) call timesec( ti ) 

    call me_mt_alloc( mt_part )
    allocate( rfun(lm_max, nrmtmax), source = 1._dp / y00 ) ! Overlap operator = (1.0/y00)*Y_{00}
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call me_mt_prepare( is, ias, l_max, zone, rfun, zzero, mt_part(:, :, ias) )
      end do
    end do

    ! N.B.: Avoiding if/else inside a loop can improve performance
    if( present( mathcalB ) .and. present( mathcalH ) ) then
      do ik = first_kpt, last_kpt
        call calculate_overlap_ik( ik, Gkset, apwalm(:, :, :, :, ik-shift), mt_part, &
          S%array(:, :, ik), mathcalB(:, :, :, :, ik-shift), &
          mathcalH(:, :, :, :, ik-shift), a_tot%components, p_MT(:, :, :, :, ik-shift) )
      end do
    else if( present( mathcalB ) ) then
      do ik = first_kpt, last_kpt
        call calculate_overlap_ik( ik, Gkset, apwalm(:, :, :, :, ik-shift), mt_part, &
          S%array(:, :, ik), mathcalB(:, :, :, :, ik-shift), &
          p_MT_ik=p_MT(:, :, :, :, ik-shift) )
      end do
    else if( present( mathcalH ) ) then
      do ik = first_kpt, last_kpt
        call calculate_overlap_ik( ik, Gkset, apwalm(:, :, :, :, ik-shift), mt_part, &
          S%array(:, :, ik), mathcalH_ik=mathcalH(:, :, :, :, ik-shift), &
          a_tot=a_tot%components, p_MT_ik=p_MT(:, :, :, :, ik-shift) )
      end do
    else
      do ik = first_kpt, last_kpt
        call calculate_overlap_ik( ik, Gkset, apwalm(:, :, :, :, ik-shift), mt_part, &
          S%array(:, :, ik) )
      end do
    end if

    if( present( t_overlap ) ) call timesec_RTTDDFT( ti, t_overlap%total )
  end subroutine

  !> Calculate the overlap matrix for a given k-point.
  subroutine calculate_overlap_ik( ik, Gkset, apwalm, mt_part, &
      overlap_ik, mathcalB_ik, mathcalH_ik, a_tot, p_MT_ik )
    !> ik: the index of the k-point considered
    integer(i32), intent(in) :: ik
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Matching coefficients of the (L)APWs at the current k-point (ngkmax, apwordmax, lmmaxapw, natmtot)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :)
    !> MT part of the overlap matrix
    complex(dp), contiguous, intent(in) :: mt_part(:, :, :)
    !> Overlap matrix (of basis functions) at the current k-point (nmatmax, nmatmax)
    complex(dp), contiguous, intent(out) :: overlap_ik(:, :)
    !> mathcalB matrix, given a k-point 
    complex(dp), contiguous, optional, intent(out) :: mathcalB_ik(:, :, :, :)
    !> mathcalH matrix, given a k-point 
    complex(dp), contiguous, optional, intent(out) :: mathcalH_ik(:, :, :, :)
    !> `x`, `y`, and `z` components of the (total) vector potential
    real(dp), optional, intent(in) :: a_tot(n_cartesian)
    !> Muffin-tin part of the Momentum matrix (nmatmax, nmatmax, 3, natmtot) for the k-point `ik`
    complex(dp), contiguous, intent(in), optional :: p_MT_ik(:, :, :, :)
    
    integer(i32) :: i, is, ia, ias, n_planewaves
    complex(dp), allocatable :: tmp(:, :)

    CALL_ASSERT( size( overlap_ik, 1 ) == size( overlap_ik, 2), "overlap_ik must be a square matrix" )
    if ( present( mathcalB_ik ) .or. present( mathcalH_ik ) ) then
      CALL_ASSERT( present(p_MT_ik),  'p_MT_ik must be present' )
    end if
    if ( present( mathcalH_ik ) ) then
      CALL_ASSERT( present(a_tot), 'a_tot must be passed as argument when calculate_mathcalH is .True.' )
    end if
    if ( present( mathcalH_ik ) ) mathcalH_ik(:, :, :, :) = zzero
    if ( present( mathcalB_ik ) ) mathcalB_ik(:, :, :, :) = -zi*p_MT_ik

    n_planewaves = Gkset%ngk(1, ik)
    overlap_ik = zzero
    
    ! It is better to split the two cases (even with some code duplication): 
    ! - to avoid an "if" inside the double loop over atoms
    ! - to avoid allocating "tmp" when not needed (this can be a large array)
    if( present( mathcalB_ik ) .or. present( mathcalH_ik ) ) then
      allocate( tmp, source = overlap_ik )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, n_planewaves, apwalm(:, :, :, ias), zone, &
            mt_part(:, :, ias), zzero, tmp )
          overlap_ik = overlap_ik + tmp
          if( present( mathcalB_ik ) ) call update_mathcalB_ik_ias( n_planewaves, &
            tmp, Gkset%vgkc(:, :, 1, ik), mathcalB_ik(:, :, :, ias) )
          if( present( mathcalH_ik ) ) call update_mathcalH_ik_ias( n_planewaves, &
            tmp, Gkset%vgkc(:, :, 1, ik), p_MT_ik(:, :, :, ias), a_tot, mathcalH_ik(:, :, :, ias) )
        end do
      end do
    else 
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, n_planewaves, apwalm(:, :, :, ias), zone, &
            mt_part(:, :, ias), zone, overlap_ik )
        end do
      end do
    end if
    call me_ir_mat( Gkset, ik, Gkset, ik, zone, cfunig, zone, overlap_ik )
    ! Fill other elements, so that the rest of overlap is the identity matrix
    do i = size( overlap_ik, 1 ), 1, -1
      if( all_zero(overlap_ik(i, i)) ) then
        overlap_ik(i, i) = zone
      else
        exit
      end if
    end do
    
  end subroutine

  !> (private) Update the `mathcalB` matrix for a given atom and k-point
  subroutine update_mathcalB_ik_ias( n_pw, overlap_ik_ias, gplusk_cart, mathcalB_ik_ias )
    !> Number of plane waves
    integer(i32), intent(in) :: n_pw
    !> Overlap MT part for a given atom and k-point
    complex(dp), contiguous, intent(in) :: overlap_ik_ias(:, :)
    !> G+k vector in Cartesian coordinates (3, ngp)
    real(dp), contiguous, intent(in) :: gplusk_cart(:, :)
    !> mathcalB matrix, given an atom and k-point 
    complex(dp), contiguous, intent(inout) :: mathcalB_ik_ias(:, :, :)

    integer(i32) :: i, i_cart
    associate( m => size(overlap_ik_ias, 1) )
      CALL_ASSERT( size(overlap_ik_ias, 2) == m, "overlap_ik_ias must be square" )
      CALL_ASSERT( size(gplusk_cart, 1) == n_cartesian, "gplusk_cart must have n_cartesian elements along 1st dim" )
      CALL_ASSERT( size(gplusk_cart, 2) >= n_pw, "gplusk_cart must at least n_pw elements along 2nd dim" )
      CALL_ASSERT( all(shape(mathcalB_ik_ias) == [m, m, n_cartesian]), "mathcalB_ik_ias must have shape [m, m, n_cartesian]" )
      CALL_ASSERT( n_pw <= m, "n_pw must be <= m" )
      do i_cart = 1, n_cartesian
        do i = 1, n_pw
          mathcalB_ik_ias(1:m, i, i_cart) = mathcalB_ik_ias(1:m, i, i_cart) + &
            zi*gplusk_cart(i_cart, i)*overlap_ik_ias(1:m, i)
        end do
      end do
    end associate
  end subroutine

  !> (private) Update the mathcalH matrix for a given atom and k-point
  subroutine update_mathcalH_ik_ias( n_pw, overlap_ik_ias, gplusk_cart, p_MT_ias_ik, a_tot, mathcalH_ik_ias )
    !> Number of plane waves
    integer(i32), intent(in) :: n_pw
    !> Overlap MT part for a given atom and k-point
    complex(dp), contiguous, intent(in) :: overlap_ik_ias(:, :)
    !> G+k vector in Cartesian coordinates (3, ngp)
    real(dp), contiguous, intent(in) :: gplusk_cart(:, :)
    !> MT-part of the momentum matrix, for a given atom and k-point
    complex(dp), contiguous, intent(in) :: p_MT_ias_ik(:, :, :)
    !> `x`, `y`, and `z` components of the (total) vector potential
    real(dp), intent(in) :: a_tot(n_cartesian)
    !> mathcalB matrix, given an atom and k-point 
    complex(dp), intent(inout) :: mathcalH_ik_ias(:, :, :)

    integer(i32) :: i, j, i_cart
    real(dp) :: fact, aux, a_scaled(n_cartesian)

    fact = dot_product( a_tot, a_tot ) / (2._dp * c**2)
    a_scaled = a_tot / c
    associate( m => size(overlap_ik_ias, 1) )
      CALL_ASSERT( size(overlap_ik_ias, 2) == m, "overlap_ik_ias must be square" )
      CALL_ASSERT( size(gplusk_cart, 1) == n_cartesian, "gplusk_cart must have n_cartesian elements along 1st dim" )
      CALL_ASSERT( size(gplusk_cart, 2) >= n_pw, "gplusk_cart must at least n_pw elements along 2nd dim" )
      CALL_ASSERT( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "mathcalH_ik_ias must have shape [m, m, n_cartesian]" )
      CALL_ASSERT( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "p_MT_ias_ik must have shape [m, m, n_cartesian]" )
      CALL_ASSERT( n_pw <= m, "n_pw must be <= m" )
      do i_cart = 1, n_cartesian
        do i = 1, n_pw
          do j = 1, n_pw
            aux = gplusk_cart(i_cart, i)-gplusk_cart(i_cart, j)
            mathcalH_ik_ias(j, i, i_cart) = mathcalH_ik_ias(j, i, i_cart) + &
              zi*aux*( fact*overlap_ik_ias(j, i) + ( &
              a_scaled(1)*p_MT_ias_ik(j, i, 1) + &
              a_scaled(2)*p_MT_ias_ik(j, i, 2) + &
              a_scaled(3)*p_MT_ias_ik(j, i, 3) ) )
          end do
          aux = gplusk_cart(i_cart, i)
          do j = n_pw + 1, m
            mathcalH_ik_ias(j, i, i_cart) = mathcalH_ik_ias(j, i, i_cart) + &
              zi*aux*( fact*overlap_ik_ias(j, i) + ( &
              a_scaled(1)*p_MT_ias_ik(j, i, 1) + &
              a_scaled(2)*p_MT_ias_ik(j, i, 2) + &
              a_scaled(3)*p_MT_ias_ik(j, i, 3) ) )
          end do
        end do
        do i = n_pw + 1, m
          do j = 1, n_pw
            aux = -gplusk_cart(i_cart, j)
            mathcalH_ik_ias(j, i, i_cart) = mathcalH_ik_ias(j, i, i_cart) + &
              zi*aux*( fact*overlap_ik_ias(j, i) + ( &
              a_scaled(1)*p_MT_ias_ik(j, i, 1) + &
              a_scaled(2)*p_MT_ias_ik(j, i, 2) + &
              a_scaled(3)*p_MT_ias_ik(j, i, 3) ) )
          end do
        end do
      end do
    end associate
  end subroutine
end module rttddft_Overlap
