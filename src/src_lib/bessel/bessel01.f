! Bessel J_0(x),J_1(x),Y_0(x),Y_1(x),
!        I_0(x),I_1(x),K_0(x),K_1(x) functions 
!        in multiple precision
!
! usage1
!     call bessel01(x,1,eps,j_0,j_1,y_0,y_1)
! usage2
!     call bessel01(x,-1,eps,i_0,i_1,k_0,k_1)
! input parameters
!     x   : argument (x > 0)
!     eps : error requested
!
      subroutine bessel01(x, l, eps, s, t, u, v)
      implicit real*8 (a - h, o - z)
      pi = 4 * atan(1.0d0)
      euler = 0.577215664901532860606512090082402431042159335939
     &    9235988057672348848677267776646709369470632917467495d0
      prc = -log(eps)
      if (x .lt. prc * 0.6d0) then
          s = eps
          t = 0
          u = 0
          v = 0
          w = 2 / x
          y = 0
          z = euler - log(w)
          if (x .gt. eps) then
              k = int(prc / 15 + ((2 + l) / 6.0d0) * x + 
     &            ((prc * prc / 7) * x) ** (1.0d0 / 3))
          else
              k = 1
          end if
          if (l .gt. 0) then
              do n = 2 * k, 2, -2
                  y = s + y
                  u = s / n - u
                  v = t * (n + 1) / (n * (n + 2)) - v
                  t = s * n * w - t
                  s = t * (n - 1) * w - s
              end do
              y = 1 / (2 * y + s)
              s = s * y
              t = t * y
              u = (4 * y * u + s * z) * (2 / pi)
              v = (4 * y * v + t * z - s * w / 2 - t) * (2 / pi)
          else
              do n = 2 * k, 2, -2
                  y = s + y
                  u = s / n + u
                  t = s * n * w + t
                  s = t * (n - 1) * w + s
              end do
              y = cosh(x) / (2 * y + s)
              s = s * y
              t = t * y
              if (x .lt. 1.5d0) then
                  u = 4 * y * u - s * z
              else
                  v = -x / 2
                  w = v
                  u = exp(-x) / 2
                  h = 9.0d0 / (prc + x)
                  y = exp(-h)
                  k = int(log(2 + 2 * prc / x) / h)
                  do n = 1, k
                      v = v / y
                      w = w * y
                      u = u + exp(v + w)
                  end do
                  u = u * h
              end if
              v = (1 / x - t * u) / s
          end if
      else
          w = 1 / (4 * x)
          k = int(prc * 0.12d0)
          do m = -1, 3, 4
              v = sqrt(8 / pi * w)
              y = v
              z = 0
              do n = 2, 8 * k + 6, 4
                  v = v * (w * (m * 1.0d0 / n - n + 2))
                  z = v - z * l
                  v = v * (w * (m * 1.0d0 / (n + 2) - n))
                  y = v - y * l
              end do
              if (m .eq. -1) then
                  t = y
                  u = z
              end if
          end do
          if (l .gt. 0) then
              v = sin(x - pi / 4)
              w = cos(x - pi / 4)
              s = v * u + w * t
              u = v * t - w * u
              t = v * y - w * z
              v = -v * z - w * y
          else
              v = exp(x) / 2
              w = pi / (4 * v)
              s = v * (t - u)
              u = w * (t + u)
              t = v * (y - z)
              v = w * (y + z)
          end if
      end if
      end
!
