module bse_utils
  
  implicit none

  contains 

  !> Return the coefficients for the type of the BSE calculation from a character string
  subroutine bse_type_to_bool(bse_type, calculate_vexc, calcualte_wscr)
    !> BSE type (IP, RPA, singlet, triplet)
    character(*), intent(in) :: bse_type

    logical, intent(out) :: calculate_vexc, calcualte_wscr

    select case(bse_type)
      case("singlet")
        calculate_vexc = .true.
        calcualte_wscr = .true.
      
      case("triplet")
        calculate_vexc = .false.
        calcualte_wscr = .true.

      case("RPA")
        calculate_vexc = .true.
        calcualte_wscr = .false.

      case("IP")
        calculate_vexc = .false.
        calcualte_wscr = .false.
      
    end select

  end subroutine bse_type_to_bool

end module bse_utils