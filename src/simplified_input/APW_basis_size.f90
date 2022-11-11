module APW_basis_size
    use modmpi, only: mpiglobal
    use precision, only: dp
    use errors_warnings, only: terminate_if_false
    use predefined_rgkmax, only : get_predefined_rgkmax
    
    implicit none 
    private
    public :: determine_rgkmax, & 
              determine_APWprecision
            
    contains
        
    !> Given the \(\text{APWprecision}\) value and the atomic number of the species which has the smallest muffin-tin radius \(\text{Z}_{\text{RMTmin}}\),
    !> the \(\text{rgkmax}\) value for the calculation is determined.
    !>
    !> For the given atomic number \(\text{Z}_{\text{RMTmin}}\) a predefined rgkmax is extracted from an array (see module predefined_rgkmax), 
    !> which contains predefined rgkmax values: 
    !> \[ 
    !>     \text{predef_rgkmax}_{\text{RMTmin}} = \textbf{predef_rgkmax} \left[ \text{Z}_{\text{RMTmin}}\right]
    !> \]
    !> The rgkmax value for the calculation is then calculated using the following formula:
    !> \[ 
    !>     \text{rgkmax} = \text{APWprecision} \cdot \text{predef_rgkmax}_{\text{RMTmin}}
    !> \]
    function determine_rgkmax(APWprecision, z_rmt_min) result(rgkmax)
        !> APWprecision
        real(dp), intent(in):: APWprecision
        !> Atomic number for species with smallest muffin-tin radius
        real(dp), intent(in):: z_rmt_min
        
        !> Predefined rgkmax value for species with smallest muffin-tin radius 
        real(dp) :: rgkmax_rmt_min
        !> Calculated rgkmax value 
        real(dp):: rgkmax

        ! Label for atomic number
        character(10) :: z_label

        rgkmax_rmt_min = get_predefined_rgkmax(z_rmt_min)
        write(z_label,'(I3)') idnint(Abs(z_rmt_min))
        call terminate_if_false(mpiglobal, rgkmax_rmt_min /= -1.0_dp ,"Error (APW_basis_size):&
                                & rgkmax cannot be calculated using the APWprecision attribute because a predefined&
                                & rgkmax value for the given atomic number"// trim(z_label) // " does not exist.& 
                                & Remove the APWprecision attribute from your input and manually choose a rgkmax-value&
                                & for the calculation.")

        rgkmax = APWprecision * rgkmax_rmt_min

    end function determine_rgkmax
    
    !> Given the \(\text{rgkmax}\) value used for the calculation and the atomic number of the species which has the smallest muffin-tin radius
    !> \(\text{Z}_{\text{RMTmin}}\), the \(\text{APWprecision}\) value for the corresponding \(\text{rgkmax}\) is determined.
    !>
    !> For the given atomic number \(\text{Z}_{\text{RMTmin}}\) a predefined rgkmax is extracted from an array (see module predefined_rgkmax), 
    !> which contains predefined rgkmax values: 
    !> \[ 
    !>     \text{predef_rgkmax}_{\text{RMTmin}} = \textbf{predef_rgkmax} \left[ \text{Z}_{\text{RMTmin}}\right]
    !> \]
    !> The \(\text{APWprecision}\) value is then calculated using the following formula:
    !> \[ 
    !>     \text{APWprecision} = \frac{ \text{rgkmax}}{\text{predef_rgkmax}_{\text{RMTmin}}} 
    !> \]
    function determine_APWprecision(rgkmax, z_rmt_min) result(APWprecision)
        real(dp), intent(in):: rgkmax
        !> Atomic number for species with smallest muffin-tin radius 
        real(dp), intent(in):: z_rmt_min
        
        !> Predefined rgkmax value for species with smallest rmt 
        real(dp) :: rgkmax_rmt_min
        !> Calculated APWprecision value 
        real(dp):: APWprecision

        ! Label for atomic number
        character(10) :: z_label
        
        rgkmax_rmt_min = get_predefined_rgkmax(z_rmt_min)
        write(z_label,'(I3)') idnint(Abs(z_rmt_min))

        if (rgkmax_rmt_min == -1.0_dp) then
            call warning("Warning (APW_basis_size): Automatic computation of the attribute APWprecision&
            & is not possible, because a predefined rgkmax value for the given atomic number"// trim(z_label) // " does not exist.")
            APWprecision = -1.0_dp 
        else 
            APWprecision = rgkmax / rgkmax_rmt_min
        end if

    end function determine_APWprecision

end module 