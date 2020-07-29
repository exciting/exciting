
module mod_pade

    use modmain, only: zzero, zone
    implicit none

contains    

    !----------------------------------------------------------------------
    ! Calculate the pade approximant and the corresponding derivative in zz 
    ! the function f calculated at the n points z
    !----------------------------------------------------------------------
    subroutine pade_approximant(n,z,f,zz,pade,dpade)
        implicit none
        integer(4), intent(in)  :: n
        complex(8), intent(in)  :: z(n)
        complex(8), intent(in)  :: f(n)
        complex(8), intent(in)  :: zz
        complex(8), intent(out) :: pade
        complex(8), intent(out) :: dpade
        ! Local
        complex(8) :: a(n)
        complex(8) :: Az(0:n), Bz(0:n)
        complex(8) :: dAz(0:n), dBz(0:n)
        integer :: i
        
        call calculate_pade_a(n,z,f,a)

        Az(0)  = zzero
        Az(1)  = a(1)
        Bz(0)  = zone
        Bz(1)  = zone
        dAz(0) = zzero
        dAz(1) = zzero
        dBz(0) = zzero
        dBz(1) = zzero

        do i = 1, n-1
            Az(i+1)  = Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
            Bz(i+1)  = Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
            dAz(i+1) = dAz(i)+a(i+1)*Az(i-1)+(zz-z(i))*a(i+1)*dAz(i-1)
            dBz(i+1) = dBz(i)+a(i+1)*Bz(i-1)+(zz-z(i))*a(i+1)*dBz(i-1)
        end do

        pade  = Az(n)/Bz(n)
        dpade = dAz(n)/Bz(n) - Az(n)*dBz(n)/(Bz(n)*Bz(n))

    end subroutine
   
    !----------------------------------------------------------------------
    subroutine calculate_pade_a(n,z,f,a)
        implicit none
        integer(4), intent(in)  :: n
        complex(8), intent(in)  :: z(n)
        complex(8), intent(in)  :: f(n)
        complex(8), intent(out) :: a(n)
        ! Local
        integer(4) :: i, j
        complex(8) :: g(n,n)

        g(1,1:n) = f(1:n)
  
        do i = 2, n
            do j = i, n
                g(i,j) = (g(i-1,i-1)-g(i-1,j)) / ((z(j)-z(i-1))*g(i-1,j))
            end do
        end do

        do i = 1, n
            a(i) = g(i,i)
        end do

    end subroutine calculate_pade_a

end module
