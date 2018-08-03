subroutine findequivkpt( vpl, kset, nk, isym, ik)
  use modmain
  use modinput
  use mod_kpointset
  implicit none

  real(8), intent( in) :: vpl(3)
  type( k_set), intent( in) :: kset
  integer, intent( out) :: nk
  integer, intent( out) :: isym( kset%nkpt)
  integer, intent( out) :: ik( kset%nkpt)

  integer :: iq, is, lspl, iv(3)
  real(8) :: s(3,3), v1(3), v2(3), t1

  nk = 0
  isym(:) = 0
  ik(:) = 0
  do is = 1, nsymcrys
    lspl = lsplsymc( is)
    s = dble( symlat( :, :, lspl))
    call r3mtv( s, vpl, v1)
    call r3frac( input%structure%epslat, v1, iv)
    do iq = 1, kset%nkpt
      v2(:) = kset%vkl( :, iq)
      call r3frac( input%structure%epslat, v2, iv)
      t1 = sum( abs( v1-v2))
      if( (t1 .lt. input%structure%epslat) .and. .not. (any( ik .eq. iq))) then
        nk = nk + 1
        ik( nk) = iq
        isym( nk) = is
        exit
      end if
    end do
  end do
  if( nk .eq. 0) then
    write (*, '("Error (findequivkpt): equivalent k-point not in set")')
    write (*, '("Requested k-point : ", 3G18.10)') vpl
    stop
  end if
  return
end subroutine findequivkpt
