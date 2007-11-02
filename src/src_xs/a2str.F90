
module a2str
  implicit none
contains

  character(256) function r2str(r,fmt)
    ! arguments
    real(8), intent(in) :: r
    character(*), intent(in) :: fmt
    write(r2str,fmt=fmt) r
    r2str=adjustl(r2str)
  end function r2str

end module a2str
