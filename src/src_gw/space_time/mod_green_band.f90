module mod_green_band
    use precision, only: wp
    use mod_bands, only: evalfv, nomax, numin, nstdf
    !! TEST and DELETE
    use mod_eigenvalue_occupancy, only: efermi,  nstfv 
    
    implicit none
    private
    public :: green_band

    contains

    !!-------------------------------------------------------------------------------
    !> Green' function (GF) in k-space and in the band-representation
    !> \begin{equation}
	  !>   G_{n \bf k}^{0}(\tau)= e^{-\left(\epsilon_{n \mathbf{k}}-\mu\right)\tau} 
    !>   \left [ \Theta(-\tau) \Theta(\mu-\epsilon_{n \mathbf{k}}) - 
    !>   \Theta(\tau) \Theta(\epsilon_{n \mathbf{k}}-\mu) \right ]
	  !> \end{equation}
    !> $\mu$ shifted to zero
    !!-------------------------------------------------------------------------------

    subroutine green_band(tau, ik, green_n)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> k index 
        integer, intent(in) :: ik
        !> Geen's function in band-representation
        real(wp), intent(out) :: green_n(:)
        
        !> Local variables
        integer :: i_band

        !> Occupied (nomax -> occupied band index) part with negative tau
        do i_band = 1, nomax
          green_n(i_band) = - exp(evalfv(i_band, ik) * tau)   ! old
          ! green_n(i_band) = exp(evalfv(i_band, ik) * tau)   !! according to Lucia's book
          write(2,*) tau, i_band, ik, evalfv(i_band, ik)
          ! if(abs(green_n(i_band))>150._wp) green_n(i_band) = 0.0_wp   ! *** this is done due to weird eigen-values from DFT
        enddo
            
        !> Unoccupied (numin -> first unoccupied band index) part with positive tau
        !> nstdf -> Total bands used to calculate epsilon (no. of states for dielectric function)
        do i_band = numin, nstdf
          green_n(i_band) = - exp(- evalfv(i_band, ik) * tau)
          write(2,*) tau, i_band, ik, evalfv(i_band, ik)
          ! if(abs(green_n(i_band))>150._wp) green_n(i_band) = 0.0_wp   ! *** this is done due to weird eigen-values from DFT
        enddo

        print*, 'sum G_nn', sum(green_n)
        write(3,*) sum(green_n)
        write(3,*) 'nomax, numin, nstdf', nomax, numin, nstdf
        write(3,*) '----------------'
        do i_band = 1, nstdf
          ! write(4, '(i4, 2f15.7, i4, f15.7)') i_band, evalfv(i_band, ik), tau, ik, green_n(i_band)
          write(4, '(i4, 3f15.7, i4, f15.7)') i_band, evalfv(i_band, ik), efermi, tau, ik, green_n(i_band)
        enddo
        write(4,*) '-----------------------------------', nomax, numin, nstdf
        print*, 'shape(evalfv)', shape(evalfv), 'k-point', ik 
    end subroutine green_band
end module mod_green_band