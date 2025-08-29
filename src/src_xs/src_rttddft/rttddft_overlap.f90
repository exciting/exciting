!> Module that manages the overlap matrix in RT-TDDFT
module rttddft_Overlap
  use asserts, only: assert
  use constants, only: y00, zi, zone, zzero
  use matrix_elements, only: me_mt_alloc, me_mt_prepare, me_mt_mat, me_ir_mat
  use mod_atoms, only: nspecies, natoms, idxas
  use mod_eigensystem, only: nmat
  use mod_gvector, only: cfunig
  use mod_kpointset, only: Gk_set
  use mod_muffin_tin, only: nrmtmax
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_GlobalMDVariables, only: mathcalH, mathcalB
  use rttddft_timings, only: timesec_RTTDDFT, Timing_RTTDDFT_overlap
  use rttddft_VectorPotential, only: Vector_Potential_Field

  implicit none

  private
  integer(i32), parameter :: n_cartesian = 3

  public :: update_overlap_lapw

contains
  !> In `update_overlap_lapw`, we obtain the overlap of the (L)APWs at time \( t \).
  subroutine update_overlap_lapw( first_kpt, overlap, apwalm, Gkset, t_overlap, pmatmt, &
    a_tot, update_mathcalH, update_mathcalB )
    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Overlap matrix (of basis functions) (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), contiguous, intent(out) :: overlap(:, :, first_kpt :)
    !> Matching coefficients of the (L)APWs (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, first_kpt :)
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Object that packs information about timings to update the overlap matrix
    type(Timing_RTTDDFT_overlap), optional, intent(inout) :: t_overlap
    !> Muffin-tin part of the Momentum matrix (nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(inout) :: pmatmt(:, :, :, :, first_kpt :)
    !> Total vector potential (needed if `update_mathcalH` is `.True.`)
    type(Vector_Potential_Field), optional, intent(in) :: a_tot
    !> if `.True.`, update `mathcalH`
    logical, intent(in), optional :: update_mathcalH
    !> if `.True.`, update `mathcalB`
    logical, intent(in), optional :: update_mathcalB

    integer(i32) :: ik, last_kpt, is, ia, ias
    integer(i32), parameter :: l_max = 0 ! Overlap operator is equal to (1.0/y00)*Y_{00}, it has only l=0 component
    integer(i32), parameter :: lm_max = (l_max+1)**2
    real(dp) :: ti
    real(dp), allocatable :: rfun(:, :)
    logical :: get_mathcalH, get_mathcalB
    complex(dp), allocatable :: mt_part(:, :, :)

    last_kpt = ubound( overlap, 3 )
    get_mathcalH = .False.
    if( present( update_mathcalH ) ) get_mathcalH = update_mathcalH
    if( get_mathcalH ) call assert( present(a_tot), 'a_tot must be passed as argument when update_mathcalH is .True.' )
    get_mathcalB = .False.
    if( present( update_mathcalB ) ) get_mathcalB = update_mathcalB
    if( present( t_overlap ) ) call timesec( ti ) 

    call me_mt_alloc( mt_part )
    allocate( rfun(lm_max, nrmtmax), source = 1._dp / y00 ) ! Overlap operator = (1.0/y00)*Y_{00}
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call me_mt_prepare( is, ias, l_max, zone, rfun, zzero, mt_part(:, :, ias) )
      end do
    end do

    do ik = first_kpt, last_kpt
      if ( get_mathcalB .or. get_mathcalH ) then
        call calculate_overlap_ik( ik, Gkset, apwalm(:, :, :, :, ik), mt_part, &
          nmat(1, ik), overlap(:, :, ik), get_mathcalB, get_mathcalH, a_tot%components, pmatmt(:, :, :, :, ik) )
      else
        call calculate_overlap_ik( ik, Gkset, apwalm(:, :, :, :, ik), mt_part, &
          nmat(1, ik), overlap(:, :, ik), .false., .false. )
      end if
    end do

    if( present( t_overlap ) ) call timesec_RTTDDFT( ti, t_overlap%total )
  end subroutine update_overlap_lapw

  !> Calculate the overlap matrix for a given k-point.
  subroutine calculate_overlap_ik( ik, Gkset, apwalm, mt_part, dim_overlap_ik, &
      overlap_ik, calculate_mathcalB, calculate_mathcalH, a_tot, p_MT_ik )
    !> ik: the index of the k-point considered
    integer(i32), intent(in) :: ik
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Matching coefficients of the (L)APWs at the current k-point (ngkmax, apwordmax, lmmaxapw, natmtot)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :)
    !> MT part of the overlap matrix
    complex(dp), contiguous, intent(in) :: mt_part(:, :, :)
    !> Number of LAPW+LOs for this `k-point` (equal to the "used" dimension of the overlap matrix)
    integer(i32), intent(in) :: dim_overlap_ik
    !> Overlap matrix (of basis functions) at the current k-point (nmatmax, nmatmax)
    complex(dp), contiguous, intent(out) :: overlap_ik(:, :)
    !> If `.True.`, calculate the MT contributions to the auxiliary matrix mathcalB
    logical, intent(in) :: calculate_mathcalB
    !> If `.True.`, calculate the MT contributions to the auxiliary matrix mathcalH
    logical, intent(in) :: calculate_mathcalH
    !> `x`, `y`, and `z` components of the (total) vector potential
    real(dp), optional, intent(in) :: a_tot(n_cartesian)
    !> Muffin-tin part of the Momentum matrix (nmatmax, nmatmax, 3, natmtot) for the k-point `ik`
    complex(dp), intent(in), optional :: p_MT_ik(:, :, :, :)
    
    integer(i32) :: i, is, ia, ias, n_planewaves
    complex(dp), allocatable :: tmp(:, :)

    if ( calculate_mathcalB .or. calculate_mathcalH ) call assert( present(p_MT_ik), &
      'p_MT_ik must be passed as argument when calculate_mathcalB or calculate_mathcalH is .True.' )
    if ( calculate_mathcalH ) call assert( present(a_tot), 'a_tot must be passed as argument when calculate_mathcalH is .True.' )
    if ( calculate_mathcalH ) mathcalH(:, :, :, :, ik) = zzero
    if ( calculate_mathcalB ) mathcalB(:, :, :, :, ik) = -zi*p_MT_ik

    n_planewaves = Gkset%ngk(1, ik)
    overlap_ik = zzero
    ! It is better to split the two cases (even with some code duplication): 
    ! - to avoid an "if" inside the double loop over atoms
    ! - to avoid allocating "tmp" when not needed (this can be a large array)
    if( calculate_mathcalB .or. calculate_mathcalH ) then
      allocate( tmp, source = overlap_ik )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, n_planewaves, apwalm(:, :, :, ias), zone, &
            mt_part(:, :, ias), zzero, tmp )
          overlap_ik = overlap_ik + tmp
          if( calculate_mathcalB ) call update_mathcalB_ik_ias( n_planewaves, &
            tmp, Gkset%vgkc(:, :, 1, ik), mathcalB(:, :, :, ias, ik) )
          if( calculate_mathcalH ) call update_mathcalH_ik_ias( n_planewaves, &
            tmp, Gkset%vgkc(:, :, 1, ik), p_MT_ik(:, :, :, ias), a_tot, mathcalH(:, :, :, ias, ik) )
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
    do concurrent( i = dim_overlap_ik+1:size( overlap_ik, 1 ) )
      overlap_ik(i, i) = zone
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
      call assert( size(overlap_ik_ias, 2) == m, "overlap_ik_ias must be square" )
      call assert( size(gplusk_cart, 2) == n_cartesian, "gplusk_cart must have n_cartesian elements along 1st dim" )
      call assert( size(gplusk_cart, 2) <= n_pw, "gplusk_cart must at least n_pw elements along 2nd dim" )
      call assert( all(shape(mathcalB_ik_ias) == [m, m, n_cartesian]), "mathcalB_ik_ias must have shape [m, m, n_cartesian]" )
      call assert( n_pw <= m, "n_pw must be <= m" )
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
      call assert( size(overlap_ik_ias, 2) == m, "overlap_ik_ias must be square" )
      call assert( size(gplusk_cart, 2) == n_cartesian, "gplusk_cart must have n_cartesian elements along 1st dim" )
      call assert( size(gplusk_cart, 2) <= n_pw, "gplusk_cart must at least n_pw elements along 2nd dim" )
      call assert( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "mathcalH_ik_ias must have shape [m, m, n_cartesian]" )
      call assert( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "p_MT_ias_ik must have shape [m, m, n_cartesian]" )
      call assert( n_pw <= m, "n_pw must be <= m" )
      do i_cart = 1, 3
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
