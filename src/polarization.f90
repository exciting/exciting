subroutine polarization
  real(8) :: pol(3,3)

  call macro_polarization( pol, 'calc')
  call macro_polarization( pol, 'write')

  return
end subroutine polarization
