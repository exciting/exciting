!> Functions and utilities for the Wigner \(3j\)-symbol and related
!> quantities such as Clebsch-Gordan or Gaunt coefficients.
module wigner3j_symbol
  use precision, only: dp
  use asserts, only: assert
  use constants, only: sqrt_two, fourpi

  implicit none
  private

  real(dp), parameter :: sqrthalf = sqrt_two/2
  real(dp), parameter :: isqrtfourpi = 1.d0/sqrt(fourpi)

  public :: wigner3j, clebsch_gordan, &
            gaunt_yyy, gaunt_yry, gaunt_yrr, gaunt_rrr

  contains

    !> Returns the Wigner \(3j\)-symbol. There are many equivalent definitions for the 
    !> \(3j\)-symbols. The following provides high accuracy for \(j_i \le 50\):
    !> \[
    !>   \begin{align*}
    !>     \begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix} &=(-1)^{j1+j2+m3} \\
    !>     &\times\sqrt{\frac{(j_1+m_1)!(j_2+m_2)!(j_3+m_3)!(j_3-m_3)!(j_1-m_1)!
    !>     (j_2-m_2)!}{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!(1+j_1+j_2+j_3)!}} \\
    !>     &\times\sum_{\max(0,j_2-j_3-m_1,j_1-j_3+m_2)}^{\min(j_1+j_2-j_3,j_1-m_1,j_2+m_2)}
    !>     (-1)^k\frac{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!}
    !>     {(j_3-j_1-m_2+k)!(j_3-j_2+m_1+k)!(j_1+j_2-j_3-k)!k!(j_1-m_1-k)!(j_2+m_2-k)}\,.
    !>   \end{align*}
    !> \]
    real(dp) function wigner3j( j1, j2, j3, m1, m2, m3)
      !> angular momentum quantum numbers
      integer, intent(in) :: j1, j2, j3
      !> magnetic quantum numbers
      integer, intent(in) :: m1, m2, m3

      integer :: k, k1, k2, l1, l2, l3, n1, n2
      real(dp) :: sgn, add, t1

      real(dp), external :: factnm, factr

      call assert( j1 >= 0, 'wigner3j: j1 is negative.')
      call assert( j2 >= 0, 'wigner3j: j2 is negative.')
      call assert( j3 >= 0, 'wigner3j: j3 is negative.')
      call assert( j1 <= 50, 'wigner3j: j1 is too large (>50).')
      call assert( j2 <= 50, 'wigner3j: j2 is too large (>50).')
      call assert( j3 <= 50, 'wigner3j: j3 is too large (>50).')
      call assert( abs(m1) <= j1, 'wigner3j: m1 out of range -j1 <= m1 <= j1.')
      call assert( abs(m2) <= j2, 'wigner3j: m2 out of range -j2 <= m2 <= j2.')
      call assert( abs(m3) <= j3, 'wigner3j: m3 out of range -j3 <= m3 <= j3.')

      if( (j1 == 0) .and. (j2 == 0) .and. (j3 == 0)) then
        wigner3j = 1._dp
        return
      end if

      l1 = j2 - j1 + j3
      l2 = j1 - j2 + j3
      l3 = j1 + j2 - j3

      if( (m1+m2+m3 /= 0) .or. any( [l1,l2,l3] < 0)) then
        wigner3j = 0._dp
        return
      end if

      n1 = j1 - m1
      n2 = j2 + m2
      k1 = max( 0, j2-j3-m1, j1-j3+m2)
      k2 = min( l3, n1, n2)
      sgn = (-1._dp)**(k1 + j1 + j2 + m3)
      add = 0._dp
      do k = k1, k2
        t1 = sgn * factr( l1, l1-n2+k) * factr( l2, l2-n1+k) * factr( l3, l3-k)
        add = add + t1 / (factnm( k, 1) * factnm( n1-k, 1) * factnm( n2-k, 1))
        sgn = -sgn
      end do
      t1 = factr( j1+m1, l1) * factr( j2+m2, l2) * factr( j3+m3, l3) * factr( j3-m3, 1+j1+j2+j3) * &
           factnm( j1-m1, 1) * factnm( j2-m2, 1)
      wigner3j = add * sqrt(t1)
    end function wigner3j

    !> Returns the Clebsch-Gordan coefficient using the Wigner \(3j\)-symbol.
    !> \[
    !>   C(J_1\; J_2\; J_3 | m_1\; m_2\; m_3) = (-1)^{J_1-J_2+m_3}\sqrt{2J_3+1}
    !>   \begin{pmatrix} J_1 & J_2 & J_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}\,.
    !> \]
    !> Suitable for \(J_i \le 50\).
    !> See [[wigner3j(function)]] for details.
    real(dp) function clebsch_gordan( j1, j2, j3, m1, m2, m3)
      !> angular momentum quantum numbers
      integer, intent(in) :: j1, j2, j3
      !> magnetic quantum numbers
      integer, intent(in) :: m1, m2, m3

      if( (m1+m2-m3 /= 0) .or. (j1+j2 < j3) .or. (j2+j3 < j1) .or. (j3+j1 < j2)) then
        clebsch_gordan = 0._dp
      else if( mod(j1-j2+m3, 2) == 0) then
        clebsch_gordan = sqrt( dble( 2*j3+1)) * wigner3j( j1, j2, j3, m1, m2,-m3)
      else
        clebsch_gordan = -sqrt( dble( 2*j3+1)) * wigner3j( j1, j2, j3, m1, m2,-m3)
      end if
    end function clebsch_gordan

    !> Returns the Gaunt coefficient of three complex spherical harmonics given by
    !> \[\begin{align*}
    !>   \langle Y_{l_1 m_1} | Y_{l_2 m_2} | Y_{l_3 m_3} \rangle 
    !>   &= \int_0^{2\pi} {\rm d}\varphi \int_0^\pi \sin\theta \, {\rm d}\theta \,
    !>      Y_{l_1 m_1}^\ast(\theta,\varphi)\, Y_{l_2 m_2}(\theta,\varphi)\, Y_{l_3 m_3}(\theta,\varphi) \\
    !>   &= (-1)^{m_1} \left[ \frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi} \right]^{\frac{1}{2}}
    !>      \begin{pmatrix} l_1 & l_2 & l_3 \\  0   & 0   & 0   \end{pmatrix}
    !>      \begin{pmatrix} l_1 & l_2 & l_3 \\ -m_1 & m_2 & m_3 \end{pmatrix} \;.
    !> \end{align*}\]
    !> Suitable for \(l_i \le 50\).
    !> See [[wigner3j(function)]] for details.
    real(dp) function gaunt_yyy( l1, l2, l3, m1, m2, m3)
      !> angular momentum quantum numbers
      integer, intent(in) :: l1, l2, l3
      !> magnetic quantum numbers
      integer, intent(in) :: m1, m2, m3

      integer :: j, j1, j2, j3, jh
      real(dp) :: t1

      real(dp), external :: factr, factnm

      if( m1-m2-m3 /= 0) then
        gaunt_yyy = 0._dp
        return
      end if

      j1 = l2 - l1 + l3
      j2 = l1 - l2 + l3
      j3 = l1 + l2 - l3

      if( (j1 < 0) .or. (j2 < 0) .or. (j3 < 0)) then
        gaunt_yyy = 0._dp
        return
      end if

      j = l1 + l2 + l3
      if( mod( j, 2) /= 0) then
        gaunt_yyy = 0._dp
        return
      end if

      jh = j / 2
      t1 = sqrt( dble( (2*l1+1) * (2*l2+1) * (2*l3+1)) * factr( j1, j+1) * factnm( j2, 1) * factnm( j3, 1))
      t1 = t1 * factr( jh, jh-l1) / (factnm( jh-l2, 1) * factnm( jh-l3, 1))
      gaunt_yyy = t1 * isqrtfourpi * wigner3j( l1, l2, l3, -m1, m2, m3)
      if( mod( m1+jh, 2) /= 0) gaunt_yyy = -gaunt_yyy
      return
    end function gaunt_yyy

    !> Returns the Gaunt coefficient of a real spherical harmonic sandwiched between
    !> two complex spherical harmonics, i.e. \(\langle Y_{l_1 m_1} | R_{l_2 m_2} | Y_{l_3 m_3} \rangle\).
    !> Whereby the real spherical harmonics are defined as
    !> \[
    !>   R_{lm}(\theta,\phi) = \begin{cases} 
    !>     \sqrt(2) \Re \lbrace Y_{lm}(\theta,\phi) \rbrace & m > 0 \\
    !>     \sqrt(2) \Im \lbrace Y_{lm}(\theta,\phi) \rbrace & m < 0 \\
    !>     \Re \lbrace Y_{lm}(\theta,\phi) \rbrace & m = 0
    !>   \end{cases}\;.
    !> \]
    !> Suitable for \(l_i \le 50\).
    !> See [[gaunt_yyy(function)]] for details.
    complex(dp) function gaunt_yry( l1, l2, l3, m1, m2, m3)
      !> angular momentum quantum numbers
      integer, intent(in) :: l1, l2, l3
      !> magnetic quantum numbers
      integer, intent(in) :: m1, m2, m3

      real(dp) :: t1, t2

      if( m2 > 0) then
        if( mod( m2, 2) == 0) then
          t1 = 1._dp
        else
          t1 = -1._dp
        end if
        t2 = sqrthalf * (gaunt_yyy( l1, l2, l3, m1, m2, m3) + t1 * gaunt_yyy( l1, l2, l3, m1, -m2, m3))
        gaunt_yry = cmplx( t2, 0._dp, dp)
      else if (m2 < 0) then
        if( mod( m2, 2) == 0) then
          t1 = 1._dp
        else
          t1 = -1._dp
        end if
        t2 = sqrthalf * (gaunt_yyy( l1, l2, l3, m1, m2, m3) - t1 * gaunt_yyy( l1, l2, l3, m1, -m2, m3))
        gaunt_yry = cmplx( 0._dp, -t2, 8)
      else
        gaunt_yry = cmplx( gaunt_yyy( l1, l2, l3, m1, m2, m3), 0._dp, 8)
      end if
      return
    end function gaunt_yry

    !> Returns the Gaunt coefficient of a complex and two real spherical harmonics, 
    !> i.e. \(\langle Y_{l_1 m_1} | R_{l_2 m_2} | R_{l_3 m_3} \rangle\).
    !> Suitable for \(l_i \le 50\).
    !> See [[gaunt_yry(function)]] for details.
    complex(dp) function gaunt_yrr( l1, l2, l3, m1, m2, m3)
      !> angular momentum quantum numbers
      integer, intent(in) :: l1, l2, l3
      !> magnetic quantum numbers
      integer, intent(in) :: m1, m2, m3

      real(dp) :: t1

      if( m3 > 0) then
        if( mod( m3, 2) == 0) then
          t1 = 1._dp
        else
          t1 = -1._dp
        end if
        gaunt_yrr = sqrthalf * (gaunt_yry( l1, l2, l3, m1, m2, m3) + t1 * gaunt_yry( l1, l2, l3, m1, m2, -m3))
      else if( m3 < 0) then
        if( mod( m3, 2) == 0) then
          t1 = 1._dp
        else
          t1 = -1._dp
        end if
        gaunt_yrr = sqrthalf * (gaunt_yry( l1, l2, l3, m1, m2, m3) - t1 * gaunt_yry( l1, l2, l3, m1, m2, -m3))*cmplx( 0._dp, -1._dp, dp)
      else
        gaunt_yrr = gaunt_yry( l1, l2, l3, m1, m2, m3)
      end if
      return
    end function gaunt_yrr

    !> Returns the Gaunt coefficient of three real spherical harmonics, 
    !> i.e. \(\langle R_{l_1 m_1} | R_{l_2 m_2} | R_{l_3 m_3} \rangle\).
    !> Suitable for \(l_i \le 50\).
    !> See [[gaunt_yry(function)]] for details.
    real(dp) function gaunt_rrr( l1, l2, l3, m1, m2, m3)
      !> angular momentum quantum numbers
      integer, intent(in) :: l1, l2, l3
      !> magnetic quantum numbers
      integer, intent(in) :: m1, m2, m3

      real(dp) :: t1

      if( m1 > 0) then
        if( mod( m1, 2) == 0) then
          t1 = 1._dp
        else
          t1 = -1._dp
        end if
        gaunt_rrr = sqrthalf * dble( gaunt_yrr( l1, l2, l3, m1, m2, m3) + t1 * gaunt_yrr( l1, l2, l3, -m1, m2, m3))
      else if( m1 < 0) then
        if( mod( m1, 2) == 0) then
          t1 = 1._dp
        else
          t1 = -1._dp
        end if
        gaunt_rrr = -sqrthalf * aimag( gaunt_yrr( l1, l2, l3, m1, m2, m3) - t1 * gaunt_yrr( l1, l2, l3, -m1, m2, m3))
      else
        gaunt_rrr = dble( gaunt_yrr( l1, l2, l3, m1, m2, m3))
      end if
      return
    end function gaunt_rrr

end module wigner3j_symbol
