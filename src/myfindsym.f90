subroutine myfindsym( vpl, nsym, isym)
  use modmain
  use modinput
  use mod_kpointset
  implicit none

  real(8), intent( in) :: vpl(3)
  integer, intent( out) :: nsym
  integer, intent( out) :: isym( 100)

  integer :: i, lspl, iv(3)
  real(8) :: s(3,3), v1(3), v2(3), t1

  nsym = 0
  isym = 0
  do i = 1, nsymcrys
    lspl = lsplsymc( i)
    s = dble( symlat( :, :, lspl))
    call r3mtv( s, vpl, v1)
    call r3frac( input%structure%epslat, v1, iv)
    t1 = sum( abs( v1 - vpl))
    if( t1 .lt. input%structure%epslat) then
      nsym = nsym + 1
      isym( nsym) = lspl
    end if
  end do
  return
End Subroutine
