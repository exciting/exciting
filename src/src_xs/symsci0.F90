subroutine symsci0(flag,scrnh0,scrnih0,scrnisym)
  implicit none
  ! arguments
  integer, intent(in) :: flag
  real(8), intent(in) :: scrnh0(3)
  real(8), intent(in) :: scrnih0(3)
  real(8), intent(out) :: scrnisym
  select case(flag)
  case(0)
     ! Peter's original choice in his BSE implementation:
     ! average the screening diagonal tensor components and take
     ! inverse of this value
     scrnisym=1.d0/((sum(scrnh0))/3.d0)
  case(1)
     ! Treatment found in the BSE implementation of R. Laskowski, also
     ! mentioned in [M. Hybertsen, PRB 34, 5390 (1986), p5411]
     scrnisym=sum(scrnih0)/3.d0
  case default
     write(*,*)
     write(*,'("Error(symsci0): not a valid choice for symmetrizing")')
     write(*,'(" the screened Coulomb interaction for q=0:",i8)') flag
     write(*,*)
     call terminate
  end select
end subroutine symsci0
