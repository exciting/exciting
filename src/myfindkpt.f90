subroutine myfindkpt( vpl, kset, isym, ik)
  use modmain
  use modinput
  use mod_kpointset
  implicit none

  real(8), intent( in) :: vpl(3)
  type( k_set), intent( in) :: kset
  integer, intent( out) :: isym
  integer, intent( out) :: ik

  integer :: lspl, iv(3)
  real(8) :: s(3,3), v1(3), v2(3), t1

  do isym = 1, nsymcrys
    lspl = lsplsymc( isym)
    s = dble( symlat( :, :, lspl))
    call r3mtv( s, vpl, v1)
    call r3frac( input%structure%epslat, v1, iv)
    do ik = 1, kset%nkpt
      v2(:) = kset%vkl( :, ik)
      call r3frac( input%structure%epslat, v2, iv)
      t1 = sum( abs( v1-v2))
      if( t1 .lt. input%structure%epslat) return
    end do
  end do
  write (*, '("Error (myfindkpt): equivalent k-point not in set")')
  write (*, '("Requested k-point : ", 3G18.10)') vpl
  call terminate
  return
End Subroutine
