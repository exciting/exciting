
module mod_pade
    implicit none

contains    

    ! Recursion for Pade coefficient (J.Serene)
    subroutine padecof(m,z,g,n,p)
        implicit none
        integer,    intent(in)  :: m
        complex(8), intent(in)  :: z(m)
        complex(8), intent(in)  :: g(m)
        integer,    intent(in)  :: n
        complex(8), intent(out) :: p(n,n)
        ! local
        integer :: i, j
        do j = 1, n
          p(1,j) = g(j)
        end do
        do j = 2, n
          do i = 2, j
            p(i,j) = (p(i-1,i-1)-p(i-1,j)) / (z(j)-z(i-1)) / p(i-1,j)
          end do
        end do
        return
    end subroutine

    ! Calculation of a green's function for a given pade-coeff-p(i,j)
    ! on the real axis e=e+i0
    subroutine gpade(m,z,n,p,e,ge)
        implicit none
        integer,    intent(in)  :: m
        complex(8), intent(in)  :: z(m)
        integer,    intent(in)  :: n
        complex(8), intent(in)  :: p(n,n)
        complex(8), intent(in)  :: e
        complex(8), intent(out) :: ge
        ! local
        integer :: i
        complex(8) :: aw(0:n), b(0:n)
        aw(0) = (0.d0,0.d0)
        aw(1) = p(1,1)
        b(0) = (1.d0,0.d0)
        b(1) = (1.d0,0.d0)
        do i = 1, n-1
          aw(i+1) = aw(i)+(e-z(i))*p(i+1,i+1)*aw(i-1)
          b(i+1) = b(i)+(e-z(i))*p(i+1,i+1)*b(i-1)
        end do
        ge = aw(n)/b(n)
        return
    end subroutine

end module