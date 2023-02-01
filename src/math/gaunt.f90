!> Objects to calculate and loop over non-zero Gaunt coefficients
module gaunt
  use precision, only: dp
  use asserts, only: assert

  implicit none
  private

  ! ******************************************
  ! ** OBJECTS FOR NON-ZERO GAUNT COEFFICIENTS
  ! ******************************************

  !> Object carrying non-zero gaunt coefficients.
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type, abstract :: non_zero_gaunt
    private
    !> maximum value of \(l_1\), \(l_2\), and \(l_3\) for which coefficients are contained
    integer :: lmax(3) = [0,0,0]
    !> number of non-zero coefficients for given \((l_1,m_1)\) and \((l_3,m_3)\)
    integer, allocatable, public :: num(:,:)
    !> combined \((l_2,m_2)\) indices of non-zero coefficients for given \((l_1,m_1)\) and \((l_3,m_3)\)
    integer, allocatable, public :: lm2(:,:,:)

    contains
      !> check if objects contains coefficients for a given range of \(l_1\), \(l_2\), and \(l_3\)
      procedure, public, non_overridable :: check_bounds
  end type non_zero_gaunt

  !> Object carrying real non-zero gaunt coefficients.
  !> Extends [[non_zero_gaunt(type)]].
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type, extends(non_zero_gaunt), public :: non_zero_gaunt_real
    private
    !> value of non-zero coefficients for given \((l_1,m_1)\) and \((l_3,m_3)\)
    real(dp), allocatable, public :: val(:,:,:)
  end type non_zero_gaunt_real

  !> Object carrying complex non-zero gaunt coefficients.
  !> Extends [[non_zero_gaunt(type)]].
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type, extends(non_zero_gaunt), public :: non_zero_gaunt_complex
    private
    !> value of non-zero coefficients for given \((l_1,m_1)\) and \((l_3,m_3)\)
    complex(dp), allocatable :: val(:,:,:)
  end type non_zero_gaunt_complex

  abstract interface
    function real_coeff( l1, l2, l3, m1, m2, m3 ) result( c )
      integer, intent(in) :: l1, l2, l3
      integer, intent(in) :: m1, m2, m3
      real(8) :: c
    end function
    function complex_coeff( l1, l2, l3, m1, m2, m3 ) result( c )
      integer, intent(in) :: l1, l2, l3
      integer, intent(in) :: m1, m2, m3
      complex(8) :: c
    end function
  end interface

  ! ******************************************
  ! ** MODULE VARIABLES CARRYING NON-ZERO GAUNT COEFFICIENTS
  ! ******************************************

  !> Objects of non-zero gaunt coefficients of the form
  !> \( \langle Y_{l_1 m_1} | Y_{l_2 m_2} | Y_{l_3 m_3} \rangle \).
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type(non_zero_gaunt_real), public :: gaunt_coeff_yyy
  !> Objects of non-zero gaunt coefficients of the form
  !> \( \langle R_{l_1 m_1} | R_{l_2 m_2} | R_{l_3 m_3} \rangle \).
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type(non_zero_gaunt_real), public :: gaunt_coeff_rrr
  !> Objects of non-zero gaunt coefficients of the form
  !> \( \langle Y_{l_1 m_1} | R_{l_2 m_2} | Y_{l_3 m_3} \rangle \).
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type(non_zero_gaunt_complex), public :: gaunt_coeff_yry
  !> Objects of non-zero gaunt coefficients of the form
  !> \( \langle Y_{l_1 m_1} | R_{l_2 m_2} | R_{l_3 m_3} \rangle \).
  !> See [[non_zero_gaunt_yyy(function)]] for construction and usage.
  type(non_zero_gaunt_complex), public :: gaunt_coeff_yrr

  public :: non_zero_gaunt_yyy, non_zero_gaunt_yry, non_zero_gaunt_yrr, non_zero_gaunt_rrr

  contains

    ! ------------------------------------------
    ! -- NON-ZERO GAUNT OBJECT PROCEDURES
    ! ------------------------------------------

    ! constructors

    function non_zero_gaunt_real_gen( lmax1, lmax2, lmax3, coeff, tolerance ) result( gaunt )
      integer, intent(in) :: lmax1
      integer, intent(in) :: lmax2
      integer, intent(in) :: lmax3
      procedure(real_coeff) :: coeff
      real(dp), optional, intent(in) :: tolerance
      type(non_zero_gaunt_real) :: gaunt

      integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      real(dp) :: eps
      real(dp) :: gnt

      eps = 1.d-20
      if( present( tolerance ) ) eps = tolerance

      call assert( lmax1 >= 0,  'non_zero_gaunt_real_gen: lmax1 is negative.' )
      call assert( lmax2 >= 0,  'non_zero_gaunt_real_gen: lmax2 is negative.' )
      call assert( lmax3 >= 0,  'non_zero_gaunt_real_gen: lmax3 is negative.' )
      call assert( lmax1 <= 50, 'non_zero_gaunt_real_gen: lmax1 is too large (>50).' )
      call assert( lmax2 <= 50, 'non_zero_gaunt_real_gen: lmax2 is too large (>50).' )
      call assert( lmax3 <= 50, 'non_zero_gaunt_real_gen: lmax3 is too large (>50).' )

      gaunt%lmax = [lmax1, lmax2, lmax3]
      if( allocated( gaunt%num ) ) deallocate( gaunt%num )
      allocate( gaunt%num((gaunt%lmax(1)+1)**2, (gaunt%lmax(3)+1)**2), source=0 )
      if( allocated( gaunt%lm2 ) ) deallocate( gaunt%lm2 )
      allocate( gaunt%lm2(gaunt%lmax(2)+2, (gaunt%lmax(1)+1)**2, (gaunt%lmax(3)+1)**2), source=0 )
      if( allocated( gaunt%val ) ) deallocate( gaunt%val )
      allocate( gaunt%val( gaunt%lmax(2)+2, (gaunt%lmax(1)+1)**2, (gaunt%lmax(3)+1)**2), source=0._dp )

      lm3 = 0
      do l3 = 0, gaunt%lmax(3)
        do m3 = -l3, l3
          lm3 = lm3 + 1

          lm1 = 0
          do l1 = 0, gaunt%lmax(1)
            do m1 = -l1, l1
              lm1 = lm1 + 1

              lm2 = 0
              do l2 = 0, gaunt%lmax(2)
                do m2 = -l2, l2
                  lm2 = lm2 + 1

                  gnt = coeff( l1, l2, l3, m1, m2, m3 )
                  if( abs( gnt ) >= eps ) then
                    gaunt%num(lm1, lm3) = gaunt%num(lm1, lm3) + 1
                    gaunt%lm2(gaunt%num(lm1, lm3), lm1, lm3) = lm2
                    gaunt%val(gaunt%num(lm1, lm3), lm1, lm3) = gnt
                  end if

                end do
              end do

            end do
          end do

        end do
      end do
    end function non_zero_gaunt_real_gen

    function non_zero_gaunt_complex_gen( lmax1, lmax2, lmax3, coeff, tolerance ) result( gaunt )
      integer, intent(in) :: lmax1
      integer, intent(in) :: lmax2
      integer, intent(in) :: lmax3
      procedure(complex_coeff) :: coeff
      real(dp), optional, intent(in) :: tolerance
      type(non_zero_gaunt_complex) :: gaunt

      integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      real(dp) :: eps
      complex(dp) :: gnt

      eps = 1.d-20
      if( present( tolerance ) ) eps = tolerance

      call assert( lmax1 >= 0,  'non_zero_gaunt_complex_gen: lmax1 is negative.' )
      call assert( lmax2 >= 0,  'non_zero_gaunt_complex_gen: lmax2 is negative.' )
      call assert( lmax3 >= 0,  'non_zero_gaunt_complex_gen: lmax3 is negative.' )
      call assert( lmax1 <= 50, 'non_zero_gaunt_complex_gen: lmax1 is too large (>50).' )
      call assert( lmax2 <= 50, 'non_zero_gaunt_complex_gen: lmax2 is too large (>50).' )
      call assert( lmax3 <= 50, 'non_zero_gaunt_complex_gen: lmax3 is too large (>50).' )

      gaunt%lmax = [lmax1, lmax2, lmax3]
      if( allocated( gaunt%num ) ) deallocate( gaunt%num )
      allocate( gaunt%num((gaunt%lmax(1)+1)**2, (gaunt%lmax(3)+1)**2), source=0 )
      if( allocated( gaunt%lm2 ) ) deallocate( gaunt%lm2 )
      allocate( gaunt%lm2(gaunt%lmax(2)+2, (gaunt%lmax(1)+1)**2, (gaunt%lmax(3)+1)**2), source=0 )
      if( allocated( gaunt%val ) ) deallocate( gaunt%val )
      allocate( gaunt%val( gaunt%lmax(2)+2, (gaunt%lmax(1)+1)**2, (gaunt%lmax(3)+1)**2), source=cmplx( 0.0_dp, 0.0_dp, dp ) )

      lm3 = 0
      do l3 = 0, gaunt%lmax(3)
        do m3 = -l3, l3
          lm3 = lm3 + 1

          lm1 = 0
          do l1 = 0, gaunt%lmax(1)
            do m1 = -l1, l1
              lm1 = lm1 + 1

              lm2 = 0
              do l2 = 0, gaunt%lmax(2)
                do m2 = -l2, l2
                  lm2 = lm2 + 1

                  gnt = coeff( l1, l2, l3, m1, m2, m3 )
                  if( abs( gnt ) >= eps ) then
                    gaunt%num(lm1, lm3) = gaunt%num(lm1, lm3) + 1
                    gaunt%lm2(gaunt%num(lm1, lm3), lm1, lm3) = lm2
                    gaunt%val(gaunt%num(lm1, lm3), lm1, lm3) = gnt
                  end if

                end do
              end do

            end do
          end do

        end do
      end do
    end function non_zero_gaunt_complex_gen

    !> Generates object carrying non-zero gaunt coefficients of the type
    !> \(\langle Y_{l_1 m_1} | Y_{l_2 m_2} | Y_{l_3 m_3} \rangle\)
    !> and the combined \((l_2,m_2)\) index for which the coefficients are not zero.
    !>
    !> Often times sums of the form
    !> \[
    !>   g_{l_1 m_1, l_3 m_3} = \sum_{l_2, m_2} \langle Y_{l_1 m_1} | Y_{l_2 m_2} | Y_{l_3 m_3} \rangle f_{l_2 m_2}
    !> \]
    !> have to be calculated. However, most of the Gaunt coefficients are zero 
    !> and don't contribute to the sum. Therefore, sums like this can be evaluated as follows
    !>
    !>     type(non_zero_gaunt_real) :: gaunt
    !>     ...
    !>     do i = 1, gaunt%num(lm1, lm3)
    !>         g(lm1, lm3) = g(lm1, lm3) + gaunt%val(i, lm1, lm3) * f(gnt%lm2(i, lm1, lm3))
    !>     end do
    !>
    !> Before usage, one can check if the object carries the Gaunt coefficients
    !> in a requested range for \(l_1\), \(l_2\), and \(l_3\) by invoking the
    !> logical function `gaunt%check_bounds`. If the object does not contain
    !> the requested coefficiets, it can be recreated using the constructor:
    !>
    !>     if( .not. gaunt%check_bounds(lmax1,lmax2,lmax3) ) &
    !>         gaunt = non_zero_gaunt_yyy( lmax1, lmax2, lmax3 )
    !>
    function non_zero_gaunt_yyy( lmax1, lmax2, lmax3, tolerance ) result( gaunt )
      use wigner3j_symbol, only: gaunt_yyy
      !> maximum l for 1st spherical harmonic
      integer, intent(in) :: lmax1
      !> maximum l for 2nd spherical harmonic
      integer, intent(in) :: lmax2
      !> maximum l for 3rd spherical harmonic
      integer, intent(in) :: lmax3
      !> coefficients smaller than tolerance are considered zero (default: 1.d-20)
      real(dp), optional, intent(in) :: tolerance
      !> object carrying non-zero Gaunt coefficients
      type(non_zero_gaunt_real) :: gaunt

      if( present( tolerance ) ) then
        gaunt = non_zero_gaunt_real_gen( lmax1, lmax2, lmax3, gaunt_yyy, tolerance=tolerance )
      else
        gaunt = non_zero_gaunt_real_gen( lmax1, lmax2, lmax3, gaunt_yyy )
      end if
    end function non_zero_gaunt_yyy

    !> Generates object carrying non-zero gaunt coefficients of the type
    !> \(\langle Y_{l_1 m_1} | R_{l_2 m_2} | Y_{l_3 m_3} \rangle\)
    !> and the combined \((l_2,m_2)\) index for which the coefficients are not zero.
    !> See [[non_zero_gaunt_yyy(function)]] for details.
    function non_zero_gaunt_yry( lmax1, lmax2, lmax3, tolerance ) result( gaunt )
      use wigner3j_symbol, only: gaunt_yry
      !> maximum l for 1st spherical harmonic
      integer, intent(in) :: lmax1
      !> maximum l for 2nd spherical harmonic
      integer, intent(in) :: lmax2
      !> maximum l for 3rd spherical harmonic
      integer, intent(in) :: lmax3
      !> coefficients smaller than tolerance are considered zero (default: 1.d-20)
      real(dp), optional, intent(in) :: tolerance
      !> object carrying non-zero Gaunt coefficients
      type(non_zero_gaunt_complex) :: gaunt

      if( present( tolerance ) ) then
        gaunt = non_zero_gaunt_complex_gen( lmax1, lmax2, lmax3, gaunt_yry, tolerance=tolerance )
      else
        gaunt = non_zero_gaunt_complex_gen( lmax1, lmax2, lmax3, gaunt_yry )
      end if
    end function non_zero_gaunt_yry

    !> Generates object carrying non-zero gaunt coefficients of the type
    !> \(\langle Y_{l_1 m_1} | R_{l_2 m_2} | R_{l_3 m_3} \rangle\)
    !> and the combined \((l_2,m_2)\) index for which the coefficients are not zero.
    !> See [[non_zero_gaunt_yyy(function)]] for details.
    function non_zero_gaunt_yrr( lmax1, lmax2, lmax3, tolerance ) result( gaunt )
      use wigner3j_symbol, only: gaunt_yrr
      !> maximum l for 1st spherical harmonic
      integer, intent(in) :: lmax1
      !> maximum l for 2nd spherical harmonic
      integer, intent(in) :: lmax2
      !> maximum l for 3rd spherical harmonic
      integer, intent(in) :: lmax3
      !> coefficients smaller than tolerance are considered zero (default: 1.d-20)
      real(dp), optional, intent(in) :: tolerance
      !> object carrying non-zero Gaunt coefficients
      type(non_zero_gaunt_complex) :: gaunt

      if( present( tolerance ) ) then
        gaunt = non_zero_gaunt_complex_gen( lmax1, lmax2, lmax3, gaunt_yrr, tolerance=tolerance )
      else
        gaunt = non_zero_gaunt_complex_gen( lmax1, lmax2, lmax3, gaunt_yrr )
      end if
    end function non_zero_gaunt_yrr

    !> Generates object carrying non-zero gaunt coefficients of the type
    !> \(\langle R_{l_1 m_1} | R_{l_2 m_2} | R_{l_3 m_3} \rangle\)
    !> and the combined \((l_2,m_2)\) index for which the coefficients are not zero.
    !> See [[non_zero_gaunt_yyy(function)]] for details.
    function non_zero_gaunt_rrr( lmax1, lmax2, lmax3, tolerance ) result( gaunt )
      use wigner3j_symbol, only: gaunt_rrr
      !> maximum l for 1st spherical harmonic
      integer, intent(in) :: lmax1
      !> maximum l for 2nd spherical harmonic
      integer, intent(in) :: lmax2
      !> maximum l for 3rd spherical harmonic
      integer, intent(in) :: lmax3
      !> coefficients smaller than tolerance are considered zero (default: 1.d-20)
      real(dp), optional, intent(in) :: tolerance
      !> object carrying non-zero Gaunt coefficients
      type(non_zero_gaunt_real) :: gaunt

      if( present( tolerance ) ) then
        gaunt = non_zero_gaunt_real_gen( lmax1, lmax2, lmax3, gaunt_rrr, tolerance=tolerance )
      else
        gaunt = non_zero_gaunt_real_gen( lmax1, lmax2, lmax3, gaunt_rrr )
      end if
    end function non_zero_gaunt_rrr

    ! type-bound procedures

    logical function check_bounds( this, l1, l2, l3 )
      class(non_zero_gaunt) :: this
      integer, intent(in) :: l1, l2, l3
      check_bounds = (l1 <= this%lmax(1)) .and. (l2 <= this%lmax(2)) .and. (l3 <= this%lmax(3))
    end function check_bounds

end module gaunt
