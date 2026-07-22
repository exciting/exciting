Module modmixer_special_functions
  Use modmpi, Only: terminate_if_false
  implicit none

Contains

  !> Estimate the order-dependent envelope used in Miller's algorithm.
  Real(8) Function envj(n, x)
    implicit none
    !> Bessel-function order.
    integer, intent(in) :: n
    !> Function argument.
    real(8), intent(in) :: x

    envj = 0.5d0 * log10(6.28d0 * n) - n * log10(1.36d0 * x / n)
  End Function

  !> Determine a stable starting order for downward recursion.
  Integer Function msta1(x, mp)
    implicit none
    !> Function argument.
    real(8), intent(in) :: x
    !> Requested decimal precision target.
    integer, intent(in) :: mp
    real(8) :: a0, f, f0, f1
    integer :: it, n0, n1, nn

    a0 = abs(x)
    n0 = int(1.1d0 * a0) + 1
    f0 = envj(n0, a0) - mp
    n1 = n0 + 5
    f1 = envj(n1, a0) - mp
    do it = 1, 20
      nn = n1 - (n1 - n0) / (1.0d0 - f0 / f1)
      f = envj(nn, a0) - mp
      if (abs(nn - n1) .lt. 1) exit
      n0 = n1
      f0 = f1
      n1 = nn
      f1 = f
    end do

    msta1 = nn
  End Function

  !> Refine the starting order for downward recursion at target order `n`.
  Integer Function msta2(x, n, mp)
    implicit none
    !> Function argument.
    real(8), intent(in) :: x
    !> Target order and requested decimal precision.
    integer, intent(in) :: n, mp
    real(8) :: a0, ejn, f, f0, f1, hmp, obj
    integer :: it, n0, n1, nn

    a0 = abs(x)
    hmp = 0.5d0 * mp
    ejn = envj(n, a0)

    if (ejn .le. hmp) then
      obj = mp
      n0 = int(1.1d0 * a0) + 1
    else
      obj = hmp + ejn
      n0 = n
    end if

    f0 = envj(n0, a0) - obj
    n1 = n0 + 5
    f1 = envj(n1, a0) - obj
    do it = 1, 20
      nn = n1 - (n1 - n0) / (1.0d0 - f0 / f1)
      f = envj(nn, a0) - obj
      if (abs(nn - n1) .lt. 1) exit
      n0 = n1
      f0 = f1
      n1 = nn
      f1 = f
    end do

    msta2 = nn + 10
  End Function

  !> Compute modified spherical Bessel functions \( i_l(x) \) for \( 0 \le l \le lmax \).
  Subroutine msbesseli(lmax, x, il)
    implicit none
    !> Maximum order to evaluate.
    integer, intent(in) :: lmax
    !> Function argument.
    real(8), intent(in) :: x
    !> Output values \( i_l(x) \) for all orders up to `lmax`.
    real(8), intent(out) :: il(0:lmax)
    integer :: l, m
    real(8) :: cs, f, f0, f1, i0, xi, xm

    call terminate_if_false((lmax .ge. 0) .and. (lmax .le. 50), &
 &   'Error(msbesseli): lmax out of range')
    call terminate_if_false((x .ge. 0.d0) .and. (x .le. 1.d5), &
 &   'Error(msbesseli): x out of range')

    xi = 1.d0 / x
    xm = 1.d-8

    if (x .lt. xm) then
      il = 0.0d0
      il(0) = 1.d0
    else
      il(0) = sinh(x) * xi
      if (lmax .ge. 1) il(1) = xi * (cosh(x) - sinh(x) * xi)
      if (lmax .ge. 2) then
        i0 = il(0)
        m = msta1(x, 200)
        if (m .lt. lmax) then
          m = lmax
        else
          m = msta2(x, lmax, 15)
        end if
        f0 = 0.d0
        f1 = -99.d0
        do l = m, 0, -1
          f = (2.d0 * l + 3.d0) * f1 / x + f0
          if (l .le. lmax) il(l) = f
          f0 = f1
          f1 = f
        end do
        cs = i0 / f
        do l = 0, lmax
          il(l) = cs * il(l)
        end do
      end if
    end if
  End Subroutine

  !> Compute modified spherical Bessel functions \( k_l(x) \) for \( 0 \le l \le lmax \).
  Subroutine msbesselk(lmax, x, kl)
    implicit none
    !> Maximum order to evaluate.
    integer, intent(in) :: lmax
    !> Function argument.
    real(8), intent(in) :: x
    !> Output values \( k_l(x) \) for all orders up to `lmax`.
    real(8), intent(out) :: kl(0:lmax)
    integer :: l
    real(8) :: k0, k1, kt, t1, t2, xi

    call terminate_if_false((lmax .ge. 0) .and. (lmax .le. 50), &
 &   'Error(msbesselk): lmax out of range')
    call terminate_if_false((x .ge. 0.d0) .and. (x .le. 1.d5), &
 &   'Error(msbesselk): x out of range')

    xi = 1.d0 / x
    if (x .lt. 1.d-8) then
      kl(0) = -xi
      t1 = 1.d0
      t2 = xi
      do l = 1, lmax
        t1 = t1 * dble(2 * l - 1)
        t2 = t2 * xi
        kl(l) = t2 * t1
      end do
    else
      kl(0) = exp(-x) * xi
      if (lmax .ge. 1) kl(1) = kl(0) * (1.d0 + xi)
      if (lmax .ge. 2) then
        k0 = kl(0)
        k1 = kl(1)
        do l = 2, lmax
          kt = k0 + dble(2 * l + 1) * xi * k1
          k0 = k1
          k1 = kt
          kl(l) = k1
        end do
      end if
    end if
  End Subroutine

End Module
