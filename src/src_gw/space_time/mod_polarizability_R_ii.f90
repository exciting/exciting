module mod_polarizability_R_ii
    use precision, only: wp
    use constants, only: zzero, zone
    use modgw, only: kqset, Gset
    use mod_green_R_ii, only: green_R_ii

    implicit none
    private
    public :: polarizability_R_ii

    contains

    !!---------------------------------------------------------
    !! Polarizability in I-I region in R-space
    !> \begin{equation}
    !>   P^{\mathbf R}_{\mathbf r, \mathbf r'}(\tau) = - 
    !>   G^{\mathbf R}_{\mathbf r, \mathbf r'}(\tau) 
    !>   G^{\mathbf R}_{\mathbf r, \mathbf r'}(-\tau) 
    !>\end{equation}
    !!---------------------------------------------------------

    subroutine polarizability_R_ii(tau, pola_R_ii)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Polarizability in imaginary time and R representation in I-I
        complex(wp), allocatable, intent(out) :: pola_R_ii(:,:,:)
        ! complex(wp), allocatable :: pola_R_ii(:,:,:)

        ! Local variables
        !> Running index for the Bravias vector R 
        integer :: ir
        !> Number of Bravais lattice R
        integer :: nbigR
        !> Running index for r-mesh
        integer :: ir_grid, ir_grid1
        !> Green function in imaginary time and R representation occupied part
        complex(wp), allocatable :: green_R_occ(:,:,:)
        !> Green function in imaginary time and R representation unoccupied part
        complex(wp), allocatable :: green_R_uno(:,:,:)

        ! Test and delete
        real :: t_i, t_f

        !> External routines for matrix matrix multiplication
        external :: zgemm

        nbigR = kqset%nkpt
        !!print*, '*******', kqset%nkpt, nbigR

        allocate(pola_R_ii(Gset%ngrtot, Gset%ngrtot, nbigR))
        pola_R_ii = zzero

        call green_R_ii(tau, green_R_occ, green_R_uno)

        call CPU_TIME(t_i)

        ! do ir = 1, nbigR
        !   do ir_grid = 1, Gset%ngrtot
        !     do ir_grid1 = 1, Gset%ngrtot
        !       pola_R_ii(ir_grid1, ir_grid, ir) = - green_R_uno(ir_grid1, ir_grid, ir) * &
        !                                            conjg(green_R_occ(ir_grid1, ir_grid, ir))
        !     enddo 
        !   enddo
        ! enddo 

        !!----------------------------------------------------------------------------------------------- OLD
        ! OLD 
        !! TODO: do we need $G^*(-\tau)$ or $G(-\tau)$
        do ir = 1, nbigR
          call zgemm('n', 'c', Gset%ngrtot, Gset%ngrtot, Gset%ngrtot, zone, green_R_uno(1,1,ir), &
                      Gset%ngrtot, green_R_occ(1,1,ir), Gset%ngrtot, zzero, pola_R_ii(1,1,ir), Gset%ngrtot)
        enddo
        !!----------------------------------------------------------------------------------------------- OLD

        !! Spin
        pola_R_ii = 2.0_wp * pola_R_ii   !! spin included here but (TODO:) do it properly 
        ! pola_R_ii = - 2.0_wp * pola_R_ii   !! "-" sign present in the eqn of polarizability

        call CPU_TIME(t_f)
        print*, 'time, sum real and imag pola R:', t_f-t_i, sum(real(pola_R_ii)), sum(imag(pola_R_ii))
        print*, '"-" time, sum real and imag pola R:', t_f-t_i, sum(real(pola_R_ii)), sum(imag(pola_R_ii))
        print*, 'test pola_R_ii', minval(real(pola_R_ii)), maxval(imag(pola_R_ii))

    end subroutine polarizability_R_ii
end module mod_polarizability_R_ii