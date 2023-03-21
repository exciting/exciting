!> Compute the product of two real valued functions given as a
!> real spherical harmonics expansion and add the result to
!> a third function, i.e., compute
!> \[ f_3({\bf r}) \leftarrow a\, f_1({\bf r})\, f_2({\bf r}) + b\, f_3({\bf r}) \;, \]
!>
!> where 
!> \[ f_i({\bf r}) = \sum\limits_{l=0}^{l_{i,{\rm max}}} \sum\limits_{m=-l}^l
!>                   f_{i,lm}(r)\, R_{lm}(\hat{\bf r}) \;, \]
!> with \(R_{lm}\) being real spherical harmonics.
!>
!> For the Gaunt coefficients, the property
!> \[  \langle R_{l_3 m_3} | R_{l_1 m_1} | R_{l_2 m_2} \rangle 
!>   = \langle R_{l_1 m_1} | R_{l_3 m_3} | R_{l_2 m_2} \rangle \]
!> is used.
!>
!> @note
!> In general, the expansion of the third function contains non-zero
!> expansion coeffiencts up to order
!> \( l_{3,{\rm max}} = l_{1,{\rm max}} + l_{2,{\rm max}} \).
!> The expansion of the third function can be truncated to smaller order
!> by setting
!> \( l_{3,{\rm max}} < l_{1,{\rm max}} + l_{2,{\rm max}} \).
!> @endnote
subroutine dshmul( lmax1, lmax2, lmax3, nr, a, f1, ld1, f2, ld2, b, f3, ld3)
  use precision, only: dp
  use asserts, only: assert
  use gaunt
  !> maximum angular momentum \(l\) in expansion of first function \(f_1({\bf r})\)
  integer, intent(in) :: lmax1
  !> maximum angular momentum \(l\) in expansion of second function \(f_2({\bf r})\)
  integer, intent(in) :: lmax2
  !> maximum angular momentum \(l\) in expansion of third function \(f_3({\bf r})\)
  integer, intent(in) :: lmax3
  !> number of radial points
  integer, intent(in) :: nr
  !> coefficent \(a\)
  real(dp), intent(in) :: a
  !> radial functions \(f_{1,lm}(r)\)
  real(dp), intent(in) :: f1(ld1,*)
  !> leading dimension of `f1` as allocated in calling scope
  integer, intent(in) :: ld1
  !> radial functions \(f_{2,lm}(r)\)
  real(dp), intent(in) :: f2(ld2,*)
  !> leading dimension of `f2` as allocated in calling scope
  integer, intent(in) :: ld2
  !> coefficient \(b\)
  real(dp), intent(in) :: b
  !> radial functions \(f_{3,lm}(r)\)
  real(dp), intent(inout) :: f3(ld3,*)
  !> leading dimension of `f3` as allocated in calling scope
  integer, intent(in) :: ld3

  integer :: lmmax1, lmmax2, lmmax3
  integer :: i, lm1, lm2

  real(dp), allocatable :: tmp(:), res(:,:), f1c(:,:), f2c(:,:)

  call assert( lmax1 >= 0, 'dshmul: Invalid argumemt: lmax1 must not be negative.')
  call assert( lmax2 >= 0, 'dshmul: Invalid argumemt: lmax2 must not be negative.')
  call assert( lmax3 >= 0, 'dshmul: Invalid argumemt: lmax3 must not be negative.')
  call assert( nr > 0, 'dshmul: Invalid argument: nr must be positive.')

  lmmax1 = (lmax1 + 1)**2
  lmmax2 = (lmax2 + 1)**2
  lmmax3 = (lmax3 + 1)**2

  if( b == 0._dp) then
    f3( 1:lmmax3, 1:nr) = 0._dp
  else if( b /= 1._dp) then
    f3( 1:lmmax3, 1:nr) = b * f3( 1:lmmax3, 1:nr)
  end if
  if( a == 0._dp) return

  ! make copy of input function for quicker memory access
  allocate( f1c, source=transpose(f1(1:lmmax1,1:nr)))
  allocate( f2c, source=transpose(f2(1:lmmax2,1:nr)))

  ! ensure availability of Gaunt coefficients
  if( .not. gaunt_coeff_rrr%check_bounds( lmax1, lmax3, lmax2)) &
    gaunt_coeff_rrr = non_zero_gaunt_rrr( lmax1, lmax3, lmax2)

  ! allocate temporary array for products
  allocate( tmp(nr), res(nr,lmmax3))
  res = 0._dp

  ! compute product
  do lm2 = 1, lmmax2
    do lm1 = 1, lmmax1
      tmp = f1c(:,lm1) * f2c(:,lm2)
      do i = 1, gaunt_coeff_rrr%num(lm1, lm2)
        if( gaunt_coeff_rrr%lm2(i, lm1, lm2) > lmmax3 ) exit
        call daxpy( nr, gaunt_coeff_rrr%val(i, lm1, lm2), tmp, 1, res(1, gaunt_coeff_rrr%lm2(i, lm1, lm2)), 1 )
      end do
    end do
  end do

  f3(1:lmmax3,1:nr) = f3(1:lmmax3,1:nr) + a * transpose( res)
  deallocate( tmp, res, f1c, f2c)
end subroutine dshmul

!> Compute the product of two complex valued functions given as a
!> complex spherical harmonics expansion and add the result to
!> a third function, i.e., compute
!> \[ f_3({\bf r}) \leftarrow a\, f_1({\bf r})\, f_2({\bf r}) + b\, f_3({\bf r}) \;, \]
!>
!> where 
!> \[ f_i({\bf r}) = \sum\limits_{l=0}^{l_{i,{\rm max}}} \sum\limits_{m=-l}^l
!>                   f_{i,lm}(r)\, Y_{lm}(\hat{\bf r}) \;, \]
!> with \(Y_{lm}\) being real spherical harmonics.
!>
!> For the Gaunt coefficients, the property
!> \[  \langle Y_{l_3 m_3} | Y_{l_1 m_1} | Y_{l_2 m_2} \rangle 
!>   = (-1)^{m_2} \langle Y_{l_3 m_3} | Y_{l_1 m_1} | Y_{l_2 -m_2}^\ast \rangle
!>   = (-1)^{m_2} \langle Y_{l_1 m_1} | Y_{l_3 m_3} | Y_{l_2 -m_2} \rangle \]
!> is used.
!>
!> @note
!> In general, the expansion of the third function contains non-zero
!> expansion coeffiencts up to order
!> \( l_{3,{\rm max}} = l_{1,{\rm max}} + l_{2,{\rm max}} \).
!> The expansion of the third function can be truncated to smaller order
!> by setting
!> \( l_{3,{\rm max}} < l_{1,{\rm max}} + l_{2,{\rm max}} \).
!> @endnote
subroutine zshmul( lmax1, lmax2, lmax3, nr, a, f1, ld1, f2, ld2, b, f3, ld3)
  use precision, only: dp
  use asserts, only: assert
  use gaunt
  !> maximum angular momentum \(l\) in expansion of first function \(f_1({\bf r})\)
  integer, intent(in) :: lmax1
  !> maximum angular momentum \(l\) in expansion of second function \(f_2({\bf r})\)
  integer, intent(in) :: lmax2
  !> maximum angular momentum \(l\) in expansion of third function \(f_3({\bf r})\)
  integer, intent(in) :: lmax3
  !> number of radial points
  integer, intent(in) :: nr
  !> coefficent \(a\)
  complex(dp), intent(in) :: a
  !> radial functions \(f_{1,lm}(r)\)
  complex(dp), intent(in) :: f1(ld1,*)
  !> leading dimension of `f1` as allocated in calling scope
  integer, intent(in) :: ld1
  !> radial functions \(f_{2,lm}(r)\)
  complex(dp), intent(in) :: f2(ld2,*)
  !> leading dimension of `f2` as allocated in calling scope
  integer, intent(in) :: ld2
  !> coefficient \(b\)
  complex(dp), intent(in) :: b
  !> radial functions \(f_{3,lm}(r)\)
  complex(dp), intent(inout) :: f3(ld3,*)
  !> leading dimension of `f3` as allocated in calling scope
  integer, intent(in) :: ld3

  integer :: lmmax1, lmmax2, lmmax3
  integer :: i, lm1, lm2, l2, m2, lmm2, sgn

  complex(dp), allocatable :: tmp(:), res(:,:)

  call assert( lmax1 >= 0, 'zshmul: Invalid argumemt: lmax1 must not be negative.')
  call assert( lmax2 >= 0, 'zshmul: Invalid argumemt: lmax2 must not be negative.')
  call assert( lmax3 >= 0, 'zshmul: Invalid argumemt: lmax3 must not be negative.')
  call assert( nr > 0, 'zshmul: Invalid argument: nr must be positive.')

  lmmax1 = (lmax1 + 1)**2
  lmmax2 = (lmax2 + 1)**2
  lmmax3 = (lmax3 + 1)**2

  if( b == cmplx( 0.0, 0.0, dp)) then
    f3(1:lmmax3,1:nr) = 0._dp
  else if( b /= cmplx( 1.0, 0.0, dp)) then
    f3(1:lmmax3,1:nr) = b * f3(1:lmmax3,1:nr)
  end if
  if( a == cmplx( 0.0, 0.0, dp)) return

  ! ensure availability of Gaunt coefficients
  if( .not. gaunt_coeff_yyy%check_bounds( lmax1, lmax3, lmax2)) &
    gaunt_coeff_yyy = non_zero_gaunt_yyy( lmax1, lmax3, lmax2)

  ! allocate temporary array for products
  allocate( tmp(nr), res(nr,lmmax3))
  res = 0._dp

  ! compute product
  do lm2 = 1, lmmax2
    l2 = int( sqrt( dble(lm2 - 1))) ! get l of (l,m)
    m2 = lm2 - (l2 + 1)**2 + l2     ! get m of (l,m)
    lmm2 = (l2 + 1)**2 - l2 - m2    ! get (l,-m)
    sgn = (-1)**m2
    do lm1 = 1, lmmax1
      tmp = f1( lm1, 1:nr) * f2( lm2, 1:nr)
      do i = 1, gaunt_coeff_yyy%num(lm1, lmm2)
        if( gaunt_coeff_yyy%lm2(i, lm1, lmm2) > lmmax3 ) exit
        call zaxpy( nr, cmplx(sgn*gaunt_coeff_yyy%val(i, lm1, lmm2), 0, dp), tmp, 1, res(1, gaunt_coeff_yyy%lm2(i, lm1, lmm2)), 1 )
      end do
    end do
  end do

  f3(1:lmmax3,1:nr) = f3(1:lmmax3,1:nr) + a * transpose( res)
end subroutine zshmul

!> Compute the product of the complex conjugate of a first function
!> with a second complex valued function, both given as a
!> complex spherical harmonics expansion, and add the result to
!> a third function, i.e., compute
!> \[ f_3({\bf r}) \leftarrow a\, f_1^\ast({\bf r})\, f_2({\bf r}) + b\, f_3({\bf r}) \;, \]
!>
!> where 
!> \[ f_i({\bf r}) = \sum\limits_{l=0}^{l_{i,{\rm max}}} \sum\limits_{m=-l}^l
!>                   f_{i,lm}(r)\, Y_{lm}(\hat{\bf r}) \;, \]
!> with \(Y_{lm}\) being real spherical harmonics.
!>
!> For the Gaunt coefficients, the property
!> \[  \langle Y_{l_3 m_3} | Y_{l_1 m_1}^\ast | Y_{l_2 m_2} \rangle 
!>   = \langle Y_{l_2 m_2} | Y_{l_3 m_3} | Y_{l_1 m_1} \rangle \]
!> is used.
!>
!> @note
!> In general, the expansion of the third function contains non-zero
!> expansion coeffiencts up to order
!> \( l_{3,{\rm max}} = l_{1,{\rm max}} + l_{2,{\rm max}} \).
!> The expansion of the third function can be truncated to smaller order
!> by setting
!> \( l_{3,{\rm max}} < l_{1,{\rm max}} + l_{2,{\rm max}} \).
!> @endnote
subroutine zshmulc( lmax1, lmax2, lmax3, nr, a, f1, ld1, f2, ld2, b, f3, ld3)
  use precision, only: dp
  use asserts, only: assert
  use gaunt
  !> maximum angular momentum \(l\) in expansion of first function \(f_1({\bf r})\)
  integer, intent(in) :: lmax1
  !> maximum angular momentum \(l\) in expansion of second function \(f_2({\bf r})\)
  integer, intent(in) :: lmax2
  !> maximum angular momentum \(l\) in expansion of third function \(f_3({\bf r})\)
  integer, intent(in) :: lmax3
  !> number of radial points
  integer, intent(in) :: nr
  !> coefficent \(a\)
  complex(dp), intent(in) :: a
  !> radial functions \(f_{1,lm}(r)\)
  complex(dp), intent(in) :: f1(ld1,*)
  !> leading dimension of `f1` as allocated in calling scope
  integer, intent(in) :: ld1
  !> radial functions \(f_{2,lm}(r)\)
  complex(dp), intent(in) :: f2(ld2,*)
  !> leading dimension of `f2` as allocated in calling scope
  integer, intent(in) :: ld2
  !> coefficient \(b\)
  complex(dp), intent(in) :: b
  !> radial functions \(f_{3,lm}(r)\)
  complex(dp), intent(inout) :: f3(ld3,*)
  !> leading dimension of `f3` as allocated in calling scope
  integer, intent(in) :: ld3

  integer :: lmmax1, lmmax2, lmmax3
  integer :: i, lm1, lm2

  complex(dp), allocatable :: tmp(:), res(:,:)

  call assert( lmax1 >= 0, 'zshmulc: Invalid argumemt: lmax1 must not be negative.')
  call assert( lmax2 >= 0, 'zshmulc: Invalid argumemt: lmax2 must not be negative.')
  call assert( lmax3 >= 0, 'zshmulc: Invalid argumemt: lmax3 must not be negative.')
  call assert( nr > 0, 'zshmulc: Invalid argument: nr must be positive.')

  lmmax1 = (lmax1 + 1)**2
  lmmax2 = (lmax2 + 1)**2
  lmmax3 = (lmax3 + 1)**2

  if( b == cmplx( 0.0, 0.0, dp)) then
    f3( 1:lmmax3, 1:nr) = 0._dp
  else if( b /= cmplx( 1.0, 0.0, dp)) then
    f3( 1:lmmax3, 1:nr) = b * f3( 1:lmmax3, 1:nr)
  end if
  if( a == cmplx( 0.0, 0.0, dp)) return

  ! ensure availability of Gaunt coefficients
  if( .not. gaunt_coeff_yyy%check_bounds( lmax2, lmax3, lmax1)) &
    gaunt_coeff_yyy = non_zero_gaunt_yyy( lmax2, lmax3, lmax1)

  ! allocate temporary array for products
  allocate( tmp(nr), res(nr,lmmax3))
  res = 0._dp

  ! compute product
  do lm2 = 1, lmmax2
    do lm1 = 1, lmmax1
      tmp = conjg( f1( lm1, 1:nr)) * f2( lm2, 1:nr)
      do i = 1, gaunt_coeff_yyy%num(lm2, lm1)
        if( gaunt_coeff_yyy%lm2(i, lm2, lm1) > lmmax3 ) exit
        call zaxpy( nr, cmplx(gaunt_coeff_yyy%val(i, lm2, lm1), 0, dp), tmp, 1, res(1, gaunt_coeff_yyy%lm2(i, lm2, lm1)), 1 )
      end do
    end do
  end do

  f3(1:lmmax3,1:nr) = f3(1:lmmax3,1:nr) + a * transpose( res)
end subroutine zshmulc
