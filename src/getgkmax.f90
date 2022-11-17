subroutine Get_gkmax( input, gkmax )
  use modinput, only: input_type
  use mod_atoms, only: nspecies
  use mod_muffin_tin, only: rmt
  implicit none
  type(input_type), intent(in) :: input
  real(8), intent(out) :: gkmax

  if ( (input%groundstate%isgkmax >= 1) .and. &
     & (input%groundstate%isgkmax <= nspecies) ) then
    gkmax = input%groundstate%rgkmax / rmt(input%groundstate%isgkmax)
  else
    gkmax = input%groundstate%rgkmax / 2.d0
  end if
  if ( 2.d0 * gkmax >= input%groundstate%gmaxvr + input%structure%epslat ) then
    write (*,*)
    write (*, '("Error(Get_gkmax): 2*gkmax > gmaxvr  ", 2G18.10)') &
         & 2.d0 * gkmax, input%groundstate%gmaxvr
    write (*,*)
    stop
  end if
end subroutine Get_gkmax
