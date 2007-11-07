
module a2str
  implicit none
contains

  character(80) function r2str(r,fmt)
    ! arguments
    real(8), intent(in) :: r
    character(*), intent(in) :: fmt
    write(r2str,fmt=fmt) r
    r2str=adjustl(r2str)
  end function r2str

  character(80) function i2str(i,fmt)
    ! arguments
    integer, intent(in) :: i
    character(*), intent(in) :: fmt
    write(i2str,fmt=fmt) i
    i2str=adjustl(i2str)
  end function i2str

end module a2str
