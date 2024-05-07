!> This module is supposed to be home of everything needed for a RT-TDDFT calculation with Hybrid functionals
module rttddft_hybrids
   use vx_enums, only: HYB_PBE0, HYB_HSE
   use mod_potential_and_density, only: xctype

   implicit none

   integer, save :: dimension_mixed_product_basis !dimension The matrix sizes of mixed product basis quantities depend on if the core electrons are treated as valence

   private :: dimension_mixed_product_basis
contains

!> Function identifying the use of a Hybrid functional
   logical function hybrids_used()
      if (xctype(1) == HYB_PBE0 .OR. xctype(1) == HYB_HSE) then
         hybrids_used = .true.
      else
         hybrids_used = .false.
      end if

   end function hybrids_used

   subroutine set_barecoul_basis()
   use modinput, only: input
      input%gw%barecoul%basis = input%groundstate%Hybrid%BasisBareCoulomb
   end subroutine
 

   subroutine Set_Dimension_mixed_product_basis( number_to_set )
      integer, intent(in) :: number_to_set
      dimension_mixed_product_basis = number_to_set
   end subroutine
 
   integer function Get_Dimension_mixed_product_basis()
      Get_Dimension_mixed_product_basis = dimension_mixed_product_basis
   end function

end module